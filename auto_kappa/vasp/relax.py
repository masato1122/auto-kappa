# -*- coding: utf-8 -*-
#
# relax.py
#
# StrincRelaxation is a class for a strict structure relaxation with
# Birch-Murnaghan equation of state
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
#
import sys
import os.path
import math
import numpy as np
import pandas as pd
import ase
import yaml

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun

# from auto_kappa.io.vasp import wasfinished
from auto_kappa.calculators.vasp import get_vasp_calculator, run_vasp
from auto_kappa.structure import change_structure_format
from auto_kappa.structure.two import get_thickness, get_normal_index

## Birch-Murnaghan equation of state
from pymatgen.analysis.eos import BirchMurnaghan

import matplotlib
matplotlib.use('Agg')

import logging
logger = logging.getLogger(__name__)

class StrictRelaxation():

    def __init__(self, initial_structure, dim=3, outdir="./volume"):
        """
        Args
        ----
        initial_structure : pymatgen Structure or ASE Atoms object
            Initial structure to be relaxed
        """
        if isinstance(initial_structure, Structure) == False:
            initial_structure = change_structure_format(
                initial_structure, format='pmg-Structure')
        
        self.initial_structure = initial_structure
        self.outdir = outdir
        self.outfile_yaml = outdir + "/result.yaml"
        self._optimal_volume = None
        self._dim = dim
        
        self._volumes = None
        self._energies = None
    
    @property
    def optimal_volume(self):
        if self._optimal_volume is not None:
            return self._optimal_volume
        else:
            msg = "\n Caution: Optimal volume is not yet obtained.\n"
            logger.warning(msg)
    
    @property
    def volumes(self):
        if self._volumes is not None:
            return self._volumes
        else:
            msg = "\n Caution: volumes were not yet set.\n"
            logger.warning(msg)

    @property
    def energies(self):
        if self._energies is not None:
            return self._energies
        else:
            msg = "\n Caution: energies were not yet set.\n"
            logger.warning(msg)
    
    @property
    def dim(self):
        return self._dim
    
    def with_different_volumes(
        self, 
        initial_strain_range=[-0.03, 0.05], nstrains=15, tol_strain=None, 
        kpts=[2,2,2],
        command={'mpirun': 'mpirun', 'nprocs': 1, 'vasp': 'vasp'},
        encut_factor=1.3, 
        verbosity=1,
        params_mod=None
        ):
        """ Run VASP jobs for different volumes and calculate the energy.
        
        Args
        ----
        strain_range : list of float
            The range of strain to explore (min, max)
        nstrains : int
            The number of strain points to calculate
        max_strain : float
            The maximum strain from the minimum energy point
        kpts : list of int
            The k-point grid for VASP calculations
        command : dict
            The command to run VASP
        encut_factor : float
            The factor to scale the plane-wave cutoff energy
        verbosity : int
            The verbosity level for logging
        params_mod : dict
            A dictionary of parameters to modify in the VASP input
        """
        if tol_strain is None:
            tol_strain = (initial_strain_range[1] - initial_strain_range[0]) / nstrains / 2.
        
        ### Run VASP jobs for different volumes
        line = " Strict structure optimization"
        msg = "\n" + line
        msg += "\n " + "=" * len(line)
        if verbosity > 0:
            logger.info(msg)
        
        ## volume of the initial structure
        init_vol = get_volume(self.initial_structure, dim=self.dim)
        
        ##
        max_iterations = 3
        flag = False
        count = 0
        while flag == False:
            
            ## Check the number of iterations
            if count == max_iterations:
                msg = "\n Number of iterations reached the limit (%d)." % count
                logger.warning(msg)
                break
            
            if count == 0:
                verbosity = 2
                strains = np.linspace(initial_strain_range[0], initial_strain_range[1], nstrains)
            else:
                verbosity = 1
            
            df_results = relaxation_with_different_volumes(
                self.initial_structure, strains,
                base_directory=self.outdir,
                tol_strain=tol_strain,
                #
                kpts=kpts, encut_factor=encut_factor,
                command=command,
                params_mod=params_mod,
                verbosity=verbosity,
                dim=self.dim
            )
            df_results = df_results.sort_values(by='strain')
            
            if len(df_results) < nstrains:
                msg = "\n Caution: Not all the strains were calculated."
                logger.warning(msg)
                count += 1
                continue
            
            ## Get the volume minimizing the energy
            if len(df_results) != 0:
                imin = np.argmin(df_results['energy'].values)
                opt_vol = df_results['volume'].values[imin]
            else:
                opt_vol = init_vol
            
            ## Set the next strain list
            min_vol = opt_vol * (1. + initial_strain_range[0])
            max_vol = opt_vol * (1. + initial_strain_range[1])
            s0 = (min_vol - init_vol) / init_vol
            s1 = (max_vol - init_vol) / init_vol
            strains = []
            for new_str in np.linspace(s0, s1, nstrains):
                if np.min(abs(df_results['strain'].values - new_str)) > tol_strain:
                    strains.append(float(new_str))
            
            if len(strains) == 0:
                flag = True
                logger.info("\n Optimal volume was found properly.")
                break
                
            count += 1
        
        ##
        Vs, Es = self.get_calculated_volumes_and_energies()
        self._volumes = Vs
        self._energies = Es
        
        self._fit()
        vol_strains = (Vs - self.optimal_volume) / self.optimal_volume
        min_strain = np.min(vol_strains)
        max_strain = np.max(vol_strains)
        if abs(min_strain) < 0.005 or abs(max_strain) < 0.005:
            msg = "\n Error: The applied strains are too small."
            msg += "\n Please check the result of the volume relaxation."
            logger.error(msg)
            sys.exit()
        
        self.output_structures()
        return Vs, Es
    
    def output_structures(self):
        
        logger.info("")

        outfile = self.outdir + "/POSCAR.init"
        self.initial_structure.to(filename=outfile)
        msg = " Output %s" % outfile.replace(os.getcwd(), ".")
        logger.info(msg)

        outfile = self.outdir + "/POSCAR.opt"
        self.get_optimal_structure().to(filename=outfile)
        msg = " Output %s" % outfile.replace(os.getcwd(), ".")
        logger.info(msg)
    
    def get_calculated_volumes_and_energies(self):
        init_cell = self.initial_structure.lattice.matrix
        df = _get_calculated_results(self.outdir, cell_pristine=init_cell, dim=self.dim)
        df_sort = df.sort_values(by='volume')
        self._volumes = df_sort['volume'].values
        self._energies = df_sort['energy'].values
        return self._volumes, self._energies
     
    def _fit(self):
        Vs, Es = self.get_calculated_volumes_and_energies()
        bm = BirchMurnaghan(Vs, Es)
        bm.fit()
        self._optimal_volume = bm.v0
        
    def get_fitting_error(self, type='mae'):
        bm = BirchMurnaghan(self.volumes, self.energies)
        bm.fit()
        Es_fit = bm.func(self.volumes)
        if type == 'mae':
            mae = np.sum(abs(self.energies - Es_fit)) / len(self.energies)
        else:
            logger.error(" Not yet supported.")
            sys.exit()
        return mae
    
    def plot_bm(self, figname='fig_bm.png'):
        """ Plot a result of fitting with Birch-Murnaghan EOS
        """
        Vs, Es = self.get_calculated_volumes_and_energies()
        bm = BirchMurnaghan(Vs, Es)
        bm.fit()
        
        if self.dim == 3:
            thickness = None
            modulus = bm.b0_GPa # GPa
            modulus_unit = "GPa"
        elif self.dim == 2:
            opt_struct = self.get_optimal_structure(format='pmg')
            thickness = get_thickness(opt_struct) # Angstrom
            modulus = bm.b0_GPa * thickness * 0.1 # N/m
            modulus_unit = "N/m"
        modulus_info = {'value': float(modulus), 'unit': modulus_unit, 'thickness': thickness}
        
        from auto_kappa.plot.fitting import plot_bm_result
        plot_bm_result(bm, figname=figname, dim=self.dim, modulus_info=modulus_info)
        
    def _strain_volume2length(self, s_vol):
        """ Convert volume strain to length strain
        """
        return math.pow((s_vol + 1.), 1 / self.dim) - 1
    
    def get_optimal_structure(self, format='pmg'):
        vinit = get_volume(self.initial_structure, dim=self.dim)
        s_vol = (self.optimal_volume - vinit) / vinit
        s_len = self._strain_volume2length(s_vol)
        struct_opt = get_strained_structure(
                self.initial_structure, s_len, format=format, dim=self.dim)
        
        vol_opt = get_volume(struct_opt, dim=self.dim)
        if abs(vol_opt - self.optimal_volume) > 1e-4:
            msg = "\n Error: There is a discrepancy in the optimal volumes."
            msg += "\n Volume 1: %.4f, volume 2: %.4f" % (self.optimal_volume, vol_opt)
            logger.error(msg)
            sys.exit()
        
        return struct_opt
        
    def print_results(self):
        
        Vs, Es = self.get_calculated_volumes_and_energies()
        bm = BirchMurnaghan(Vs, Es)
        bm.fit()
        
        # vinit= self.initial_structure.volume
        vinit = get_volume(self.initial_structure, dim=self.dim)
        s_vol = (vinit - self.optimal_volume) / self.optimal_volume
        s_len = self._strain_volume2length(s_vol)
        mae = self.get_fitting_error()
        
        msg = "\n"
        msg += " Minimum energy : %8.3f eV\n" % bm.e0
        if self.dim == 3:
            modulus = bm.b0_GPa # GPa
            modulus_unit = "GPa"
        elif self.dim == 2:
            opt_struct = self.get_optimal_structure(format='pmg')
            thickness = get_thickness(opt_struct)
            modulus = bm.b0_GPa * thickness * 0.1 # N/m
            modulus_unit = "N/m"
        msg += " Bulk modulus   : %8.3f %s\n" % (modulus, modulus_unit)
        msg += " Optimal volume : %8.3f A^3\n" % self.optimal_volume
        msg += " Initial volume : %8.3f A^3\n" % vinit
        msg += " Initial strain : %8.3f\n" % (s_len)
        msg += " Error (MAE)    : %8.5f eV" % (mae)
        if self.dim == 2:
            thickness = get_thickness(self.initial_structure)
            msg += "\n"
            msg += " 2D thickness   : %8.3f A" % (thickness)
        logger.info(msg)
        
        out = {}
        out['minimum_energy'] = [float(bm.e0), 'eV']
        out['bulk_modulus'] = [float(modulus), modulus_unit]
        out['optimal_volume'] = [float(self.optimal_volume), 'A^3']
        out['initial_volume'] = [float(vinit), 'A^3']
        out['initial_strain'] = [float(s_len), '-']
        out['mae'] = [float(mae), 'eV']
        if self.dim == 2:
            out['thickness'] = [float(thickness), 'A']
        
        with open(self.outfile_yaml, 'w') as f:
            yaml.dump(out, f)
            msg = "\n Output %s" % self.outfile_yaml.replace(os.getcwd(), ".")
            logger.info(msg)
        
def _get_calculated_results(base_dir, cell_pristine=None, num_max=200, dim=3):
    """ Get all the results which have been already obtained.
    
    Args
    ----
    base_dir : string
        VASP calculations were performed in base_dir+"/*".
    
    cell_pristine : ndarray, shape=(3,3)
        cell size w/o strain
    
    Return
    -------
    all_results : list of dictionary
        each element contains the info in strain.yaml file
    """
    ### get possible directories 
    dirs = []
    for ii in range(num_max):
        dir_each = base_dir + "/%d" % ii
        file_xml = dir_each + "/vasprun.xml"
        if os.path.exists(file_xml):
            dirs.append(dir_each)
    
    ### check each directory
    all_results = []
    for diri in dirs:
        
        ### check result (vasprun.xml)
        ### ver.1
        #if wasfinished(diri, filename='vasprun.xml') == False:
        #    continue
        
        ### ver.2
        file_xml = diri + "/vasprun.xml"
        try:
            vasprun = Vasprun(file_xml, parse_potcar_file=False)
        except Exception:
            continue
        
        ### Make strain.yaml if it doesn't exist and possible.
        file_yaml = diri + "/strain.yaml"
        if os.path.exists(file_yaml) == False or dim == 2:
            try:
                _make_strain_yaml(diri, cell_pristine=cell_pristine, dim=dim, verbosity=0)
            except Exception:
                continue
        
        with open(file_yaml, 'r') as yml:
            data = yaml.safe_load(yml)
            data['filename'] = "./" + os.path.relpath(file_xml, base_dir)
            all_results.append(data)
    
    if len(all_results) == 0:
        return None
    else:
        return pd.DataFrame(all_results)

def _make_strain_yaml(directory, cell_pristine=None, dim=3, verbosity=1):
    """ Make strain.yaml file """
    
    ### read file 
    file_xml = directory + "/vasprun.xml"
    vasprun = Vasprun(file_xml, parse_potcar_file=False)
    structure = vasprun.final_structure
    
    ### get strain
    l1 = np.linalg.norm(structure.lattice.matrix[0])
    l0 = np.linalg.norm(cell_pristine[0])
    strain = (l1 - l0) / l0
    
    ### get parameters
    out_data = {}
    out_data['strain'] = float(strain)
    out_data['energy'] = float(vasprun.final_energy.real)
    out_data['energy_unit'] = str(vasprun.final_energy.unit)
    out_data['volume'] = float(get_volume(structure, dim=dim))
    out_data['volume_unit'] = 'A^3'
    out_data['dim'] = int(dim)
    if dim == 2:
        thickness = get_thickness(structure)
        out_data['thickness'] = float(thickness)
        out_data['thick_unit'] = 'A'
    
    ### output file
    outfile = f"{directory}/strain.yaml"
    with open(outfile, "w") as f:
        yaml.dump(out_data, f)
        if verbosity > 0:
            msg = "\n Output ./%s" % os.path.relpath(outfile, os.getcwd())
            logger.info(msg)
    return out_data

# def _get_results_all(base_dir, keys=['strain'], cell_pristine=None, dim=3):
#     all_results = _get_calculated_results(base_dir, cell_pristine=cell_pristine, dim=dim)
#     extracted_data = []
#     for each in all_results:
#         try:
#             values = []
#             for key in keys:
#                 values.append(each.get(key, None))
#             extracted_data.append(values)
#         except Exception:
#             pass
#     return extracted_data

# def _get_calculated_strains(base_dir, cell_pristine=None, dim=3):
#     """ Get list of strains for which the energy has been already calculated. 
#     """
#     df_results = _get_calculated_results(base_dir, cell_pristine=cell_pristine, dim=dim)
    
#     if len(df_results) == 0:
#         return []
#     else:
#         strains_all = np.sort(df_results['strain'].values)
#         return strains_all
        
#         # ### remove duplicative data
#         # strains = []
#         # for i in range(len(strains_all)):
#         #     if i == 0:
#         #         strains.append(strains_all[0])
#         #     else:
#         #         if strains_all[i] - strains_all[i-1] > tol_strain:
#         #             strains.append(strains_all[i])
#         # strains = np.asarray(strains)
#         # return strains

# def _get_calculated_volumes_and_energies(
#     base_dir, keys=['volume', 'energy', 'filename'], cell_pristine=None, dim=3):
#     """ Read info from strain.yaml files and return a DataFrame. """
#     value_tmp = _get_calculated_results(base_dir, cell_pristine=cell_pristine, dim=dim)
    
#     dump = []
#     for each in value_tmp:
#         dump.append([])
#         for i, name in enumerate(keys):
#             if name == 'energy':
#                 val = each[i][0]
#             else:
#                 val = each[i]
#             dump[-1].append(val)
    
#     df = pd.DataFrame(dump, columns=keys)
#     df = df.sort_values(by='volume')
    
#     return df

def get_strained_structure(structure0, strain, format='pmg', dim=3):
    """
    Args
    ----
    structure0 : pymatgen structure obj.
        reference structure
    
    strain : float, unit=[-]
        applied strain
    
    Return
    --------
    structure : 
        strained structure with the assigened format
    """
    if isinstance(structure0, Structure) == False:
        structure0 = change_structure_format(structure0, format='pmg-Structure')
    
    cell0 = structure0.lattice.matrix.copy()
    positions = structure0.cart_coords
    symbols = []
    for el in structure0.species:
        symbols.append(el.name)
    
    ### prepare a strained structure
    if dim == 3:
        s = (1. + np.array(strain)) * np.eye(3)
    elif dim == 2:
        normal_idx = get_normal_index(structure0, base='xyz')
        s = (1. + np.array(strain)) * np.eye(3)
        s[normal_idx, normal_idx] = 1.0
    else:
        logger.error(f" dim == {dim} is not supported.")
        sys.exit()
    
    cell = np.dot(cell0.T, s).T
    positions = np.dot(positions, s.T)
    
    if format.lower().startswith('pmg') or format.lower() == 'pymatgen':
        ## pymatgen
        structure = Structure(
                lattice=cell,
                species=symbols,
                coords=positions,
                coords_are_cartesian=True
                )
    elif format.lower() == 'ase':
        ## ASE
        structure = ase.Atoms(
                cell=cell, pbc=True,
                positions=positions,
                symbols=symbols,
                )
    else:
        logger.error("\n Error: Structure format {} is not supported".format(format))
        sys.exit()
    
    return structure

def get_optimal_volume(volumes, energies, min_num_data=5):
    if len(volumes) >= min_num_data:
        bm = BirchMurnaghan(volumes, energies)
        bm.fit()
        return bm.v0
    else:
        return None

def relaxation_with_different_volumes(
        struct_init, strains,
        base_directory='./volume',
        tol_strain=1e-5,
        kpts=[2,2,2], encut_factor=1.3, 
        command={'mpirun': 'mpirun', 'nprocs': 1, 'vasp': 'vasp'},
        params_mod=None,
        verbosity=1, dim=3,
        ):
    """ Structure relaxation with the Birch-Murnaghan equation of state
    Args
    ----
    strains : array of float, shape=(nstrains), unit=[-]
        strain values to be applied to the initial structure
    
    nstrains : integer
        number of strains to be applied
    """
    ### get results in previous calculations
    df_results = _get_calculated_results(
        base_directory, cell_pristine=struct_init.lattice.matrix, dim=dim)
    if df_results is None:
        calculated_strains = []
    else:
        calculated_strains = np.sort(df_results['strain'].values)
    
    ### print already-analyzed-strains
    if len(calculated_strains) > 0:
        line = "Strains already calcuclated (%):"
        msg = "\n %s" % line
        msg += "\n " + "-" * len(line)
        msg += "\n "
        for each in calculated_strains:
            msg += "%.3f, " % (each*100.)
        logger.info(msg)
    
    ### calculate the potential energy for each strain value
    for ist, strain in enumerate(strains):
        
        if len(calculated_strains) > 0:
            if np.min(abs(calculated_strains - strain)) < tol_strain:
                continue
        
        ### print strain value
        line = "Strain: %f" % strain
        msg = "\n " + line
        msg += "\n " + "-" * len(line)
        if verbosity > 0:
            logger.info(msg)
        
        ### set the output directory
        flag = False
        outdir = None
        count = 1
        while flag == False:
            outdir = base_directory + "/%d" % int(count)
            if os.path.exists(outdir) == False:
                flag = True
            count += 1
        
        if verbosity > 0:
            logger.info("\n Output directory: %s" % (
                outdir.replace(os.getcwd(), ".")))
        
        #### prepare a strained structure
        atoms = get_strained_structure(struct_init, strain, format='ase', dim=dim)
        
        ### set calculator object
        calc = get_vasp_calculator(
                'relax-freeze', 
                directory=outdir,
                atoms=atoms,
                kpts=kpts,
                encut_scale_factor=encut_factor,
                **params_mod
                )
        
        mpirun = command.get('mpirun', 'mpirun')
        nprocs = command.get('nprocs', 1)
        if list(kpts) == [1, 1, 1]:
            vasp = command.get('vasp_gam')
        else:
            vasp = command.get('vasp', 'vasp')
        calc.command = f"{mpirun} -n {nprocs} {vasp}"
        
        ### run VASP
        run_vasp(calc, atoms, method='custodian')
        
        ### Get each result
        out_data = _make_strain_yaml(outdir, cell_pristine=struct_init.lattice.matrix, dim=dim)
        volume = out_data['volume']
        vol_unit = out_data['volume_unit']
        energy = out_data['energy']
        ene_unit = out_data['energy_unit']
        
        if verbosity > 0:
            msg = (
                f" strain(%): {strain*100:.3f}   "
                f"volume({vol_unit}): {volume:.3f}   "
                f"energy({ene_unit}): {energy:.3f}")
            logger.info(msg)
    
    ###
    init_cell = struct_init.lattice.matrix
    df_results = _get_calculated_results(base_directory, cell_pristine=init_cell, dim=dim)
    
    outfile = base_directory + "/volume_energy.csv"
    df_results.to_csv(outfile, index=False)
    msg = "\n Output %s" % outfile.replace(os.getcwd(), ".")
    logger.info(msg)
    
    return df_results

def get_volume(struct, dim=3):
    
    struct_pmg = change_structure_format(struct, format='pmg-Structure')
    vol = struct_pmg.volume
    cell = struct_pmg.lattice.matrix
    
    if dim == 2:
        norm_idx_abc = get_normal_index(struct_pmg, base='abc')
        # norm_idx_xyz = get_normal_index(struct_pmg, base='xyz')
        thickness = get_thickness(struct_pmg)
        cell_height = np.linalg.norm(cell[norm_idx_abc])
        
        ## Check cell angles
        angles = struct_pmg.lattice.angles
        angles = np.delete(angles, norm_idx_abc)
        if all(abs(angle - 90.0) < 1e-5 for angle in angles) == False:
            msg = "\n Error: The angles of the cell are not 90 degrees."
            msg += "\n Please check the structure."
            logger.error(msg)
            sys.exit()
        
        vol = vol * thickness / cell_height
    
    return vol

#from pymatgen.io.vasp import Poscar
#
#struct_init = Poscar.from_file("POSCAR.init").structure
#
#relax = StrictRelaxation(struct_init)
#Vs, Es = relax.with_different_volumes()
#relax.plot_bm(figname='fig_bm.png')
#relax.print_results()
#
#struct_opt = relax.get_optimal_structure()
#struct_opt.to(filename="POSCAR.opt")
#
