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
import ase
import glob
import yaml

import ase.io
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun

from auto_kappa.io.vasp import wasfinished
from auto_kappa.calculators.vasp import get_vasp_calculator, run_vasp

## Birch-Murnaghan equation of state
from pymatgen.analysis.eos import BirchMurnaghan

import matplotlib
matplotlib.use('Agg')

import logging
logger = logging.getLogger(__name__)

class StrictRelaxation():

    def __init__(self, initial_structure, outdir="./volume"):
        
        self.initial_structure = initial_structure
        self.outdir = outdir
        self._optimal_volume = None
        self.outfile_yaml = outdir + "/result.yaml"

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
    
    def with_different_volumes(self, 
            initial_strain_range=[-0.01, 0.05], nstrains=11,
            tol_frac_ediff=0.5,
            kpts=[2,2,2],
            command={'mpirun': 'mpirun', 'nprocs': 1, 'vasp': 'vasp'},
            encut_factor=1.3, 
            verbosity=1,
            params_mod=None
            ):
        
        ### Run VASP jobs for different volumes
        line = " Strict structure optimization"
        msg = "\n" + line
        msg += "\n " + "=" * len(line)
        if verbosity > 0:
            logger.info(msg)
         
        ### get the strain interval
        strain_interval = (
                (initial_strain_range[1] - initial_strain_range[0]) /
                (nstrains - 1)
                )
        
        ### calculate energies with the initial setting
        Vs0, Es0 = relaxation_with_different_volumes(
                self.initial_structure,
                base_directory=self.outdir,
                strain_range=initial_strain_range,
                nstrains=nstrains,
                #
                kpts=kpts, encut_factor=encut_factor,
                command=command,
                verbosity=2,
                params_mod=params_mod
                )
        
        def _get_energy_info(energies):
            emin = np.min(energies)
            emax = np.max(energies)
            de_tot = emax - emin
            return emin, emax, de_tot
        
        ### analyze smaller or larger strains
        emin, emax, de_tot = _get_energy_info(Es0)
        flag = False
        if Es0[0] - emin < de_tot * tol_frac_ediff:
            direction = 'smaller'
        elif Es0[-1] - emin < de_tot * tol_frac_ediff:
            direction = 'larger'
        else:
            flag = True
        
        strain_range = initial_strain_range
        while flag == False:
            
            n = 2
            if direction == 'smaller':
                s0 = strain_range[0] - strain_interval * n
                s1 = strain_range[0] - strain_interval
            else:
                s0 = strain_range[1] + strain_interval
                s1 = strain_range[1] + strain_interval * n
            
            strain_range = [s0, s1]
            Vs, Es = relaxation_with_different_volumes(
                    self.initial_structure,
                    base_directory=self.outdir,
                    strain_range=strain_range,
                    nstrains=n,
                    #
                    kpts=kpts, encut_factor=encut_factor,
                    command=command,
                    verbosity=1,
                    params_mod=params_mod
                    )
            
            ### check energies
            emin, emax, de_tot = _get_energy_info(Es)
            if direction == 'smaller':
                if Es[0] - emin > tol_frac_ediff * de_tot:
                    flag = True
            elif direction == 'larger':
                if Es[1] - emin > tol_frac_ediff * de_tot:
                    flag = True
        
        ###
        Vs, Es = self.get_calculated_volumes_and_energies()
        self._fit()
        self._volumes = Vs
        self._energies = Es
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
        Vs, Es = _get_calculated_volumes_and_energies(self.outdir)
        self._volumes = Vs
        self._energies = Es
        return Vs, Es
     
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
        
        from auto_kappa.plot.fitting import plot_fitting_result
        fig = plot_fitting_result(bm, figname=figname)
    
    def _strain_volume2length(self, s_vol):
        """ Convert volume strain to length strain
        """
        return math.pow((s_vol + 1.), 1/3) - 1
    
    def get_optimal_structure(self, format='pmg'):
        
        vinit = self.initial_structure.volume
        s_vol = (self.optimal_volume - vinit) / vinit
        s_len = self._strain_volume2length(s_vol)
        
        struct_opt = get_strained_structure(
                self.initial_structure, s_len, format=format)
        
        return struct_opt 
        
    def print_results(self):
        
        Vs, Es = self.get_calculated_volumes_and_energies()
        bm = BirchMurnaghan(Vs, Es)
        bm.fit()
        
        vinit= self.initial_structure.volume
        s_vol = (vinit - self.optimal_volume) / self.optimal_volume
        s_len = self._strain_volume2length(s_vol)
        mae = self.get_fitting_error()
        
        msg = "\n"
        msg += " Minimum energy : %8.3f eV\n" % bm.e0
        msg += " Bulk modulus   : %8.3f GPa\n" % bm.b0_GPa
        msg += " Optimal volume : %8.3f A^3\n" % self.optimal_volume
        msg += " Initial volume : %8.3f A^3\n" % vinit
        msg += " Initial strain : %8.3f\n" % (s_len)
        msg += " Error (MAE)    : %8.5f eV" % (mae)
        logger.info(msg)
        
        out = {}
        out['minimum_energy'] = [float(bm.e0), 'eV']
        out['bulk_modulus'] = [float(bm.b0_GPa), 'GPa']
        out['optimal_volume'] = [float(self.optimal_volume), 'A^3']
        out['initial_volume'] = [float(vinit), 'A^3']
        out['initial_strain'] = [float(s_len), '-']
        out['mae'] = [float(mae), 'eV']
         
        with open(self.outfile_yaml, 'w') as f:
            yaml.dump(out, f)
            msg = "\n Output %s" % self.outfile_yaml.replace(os.getcwd(), ".")
            logger.info(msg)
    
def _get_calculated_results(base_dir, cell_pristine=None, num_max=1000):
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
        if os.path.exists(file_yaml) == False:
            try:
                _make_strain_yaml(diri, cell_pristine=cell_pristine)
            except Exception:
                continue
        
        with open(file_yaml, 'r') as yml:
            data = yaml.safe_load(yml)
            all_results.append(data)
    
    return all_results

def _make_strain_yaml(directory, cell_pristine=None):
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
    out_data['energy'] = [
            float(vasprun.final_energy.real),
            str(vasprun.final_energy.unit)]
    out_data['volume'] = float(structure.volume)
    
    ### output file
    outfile = directory + "/strain.yaml"
    with open(outfile, "w") as f:
        yaml.dump(out_data, f)
        msg = "\n Output %s" % outfile
        logger.info(msg)
    
    return 0

def _get_results_all(base_dir, keys=['strain'], cell_pristine=None):
    
    all_results = _get_calculated_results(base_dir, cell_pristine=cell_pristine)
    extracted_data = []
    for each in all_results:
        try:
            values = []
            for key in keys:
                values.append(each[key])
            extracted_data.append(values)
        except Exception:
            pass
    return extracted_data

def _get_calculated_strains(
        base_dir, key='strain', cell_pristine=None, tol_strain=1e-9):
    """ Get list of strains for which the energy has been already calculated. 
    """
    value_tmp = _get_results_all(
            base_dir, keys=[key], 
            cell_pristine=cell_pristine)
    
    if len(value_tmp) == 0:
        return []
    else:
        strains_all = np.sort(np.asarray(value_tmp)[:,0])
        
        ### remove duplicative data
        strains = []
        for i in range(len(strains_all)):
            if i == 0:
                strains.append(strains_all[0])
            else:
                if strains_all[i] - strains_all[i-1] > tol_strain:
                    strains.append(strains_all[i])
        strains = np.asarray(strains)
        return strains

def _get_calculated_volumes_and_energies(
        base_dir, keys=['volume', 'energy'], cell_pristine=None):
    
    value_tmp = _get_results_all(
            base_dir, keys=keys, 
            cell_pristine=cell_pristine)
    
    Vs = []
    Es = []
    for each in value_tmp:
        Vs.append(each[0])
        Es.append(each[1][0])
    
    isort = np.argsort(Vs)
    return (
            np.asarray(Vs)[isort], 
            np.asarray(Es)[isort], 
            )

def get_strained_structure(structure0, strain, format='pmg'):
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
    cell0 = structure0.lattice.matrix.copy()
    scaled_positions = structure0.frac_coords
    symbols = []
    for el in structure0.species:
        symbols.append(el.name)
    
    ### prepare a strained structure
    s = (1. + np.array(strain)) * np.eye(3)
    cell = np.dot(cell0.T, s).T
    
    if (format.lower() == 'pmg' or
            format.lower() == 'pymatgen'
            ):
        ## pymatgen
        structure = Structure(
                lattice=cell,
                species=symbols,
                coords=scaled_positions,
                )
    elif format.lower() == 'ase':
        ## ASE
        structure = ase.Atoms(
                cell=cell, pbc=True,
                scaled_positions=scaled_positions,
                symbols=symbols,
                )
    
    return structure

def relaxation_with_different_volumes(
        struct_init, strain_range=[-0.01, 0.05], nstrains=11, kpts=[2,2,2],
        base_directory='./volume', encut_factor=1.3, 
        command={'mpirun': 'mpirun', 'nprocs': 1, 'vasp': 'vasp'},
        tol_strain=1e-5,
        verbosity=1,
        params_mod=None
        ):
    """ Structure relaxation with the Birch-Murnaghan equation of state
    Args
    ----
    struct_init : pymatgen Structure

    strain_range : array of float, shape=(2), unit=[-]
    
    nstrains : integer
        number of strains to be applied
    """
    ### list of strains to be analyzed
    strains = np.linspace(strain_range[0], strain_range[1], nstrains)
    
    ### get the calculated strains
    calculated_strains = _get_calculated_strains(
            base_directory, 
            cell_pristine=struct_init.lattice.matrix,
            tol_strain=tol_strain*0.1)
    
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
        #structure = get_strained_structure(struct_init, strain, format='pmg')
        atoms = get_strained_structure(struct_init, strain, format='ase')
        
        ### set calculator object
        calc = get_vasp_calculator(
                'relax-freeze', 
                directory=outdir,
                atoms=atoms,
                kpts=kpts,
                encut_scale_factor=encut_factor,
                **params_mod
                )
        calc.command = "%s -n %d %s" % (
                command['mpirun'],
                command['nprocs'],
                command['vasp'],
                )
        
        ### run VASP
        run_vasp(calc, atoms, method='custodian')
        
        ### output data
        vasprun = Vasprun(outdir + "/vasprun.xml", parse_potcar_file=False)
        
        volume = float(atoms.get_volume())
        energy = float(vasprun.final_energy.real)
        out_data = {}
        out_data['strain'] = float(strain)
        out_data['energy'] = [energy, str(vasprun.final_energy.unit)]
        out_data['volume'] = float(atoms.get_volume())
        
        ### output yamle file
        outfile = outdir + "/strain.yaml"
        with open(outfile, "w") as f:
            yaml.dump(out_data, f)
        
        if verbosity > 0:
            msg = "\n strain: %f  volume: %f  energy: %f" % (
                    strain, volume, energy)
            logger.info(msg)
    
    ###
    Vs, Es = _get_calculated_volumes_and_energies(base_directory)
    return Vs, Es

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

