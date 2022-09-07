# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

from auto_alamode2.io.phonondb import Phonondb
from auto_alamode2.apdb import ApdbVasp
from auto_alamode2.almcalc import AlmCalc
from auto_alamode2 import output_directories

parser = OptionParser()

parser.add_option("-d", "--directory", dest="directory", type="string",
        default="../mp-149", help="directory of phonondb")

parser.add_option("--mpid", dest="mpid", type="string",
        default="mp-149", help="material ID [mp-149]")

parser.add_option("--cutoff3", dest="cutoff3", type="float",
        default=4.3, help="cutoff3, unit=Ang [4.3]")

parser.add_option("-n", "--nprocs", dest="nprocs", type="int",
        default=2, help="nprocs [2]")

(options, args) = parser.parse_args()

### Read data of phonondb
phdb = Phonondb(options.directory, mpid=options.mpid)

pmat = phdb.primitive_matrix
smat = phdb.scell_matrix

### Set output directoreis
out_dirs = {}
for k1 in output_directories.keys():
    values1 = output_directories[k1]
    if type(values1) == str:
        out_dirs[k1] = './' + options.mpid + '/' + values1
    else:
        out_dirs[k1] = {}
        for k2 in values1.keys():
            values2 = values1[k2]
            out_dirs[k1][k2] = './' + options.mpid + '/' + values2

### command to run a VASP job
vasp_cmd = {'mpirun': 'mpirun', 'nprocs': options.nprocs, 'vasp': 'vasp'}

### Set ApdbVasp object
apdb = ApdbVasp(
        phdb.get_unitcell(format='ase'), 
        primitive_matrix=pmat,
        scell_matrix=smat,
        command=vasp_cmd)

### Relaxation calculation
mode = 'relax'
apdb.run_vasp(mode, 
        out_dirs[mode],
        phdb.get_kpoints(mode=mode).kpts[0],
        print_params=True
        )

### Born effective charge
if phdb.nac == 1:
    mode = 'nac'
    apdb.run_vasp(mode,
            out_dirs[mode],
            phdb.get_kpoints(mode=mode).kpts[0],
            print_params=True
            )

### Set AlmCalc
almcalc = AlmCalc(
        apdb.primitive,
        mpid=options.mpid,
        primitive_matrix=pmat,
        scell_matrix=smat,
        cutoff2=-1, cutoff3=options.cutoff3, 
        nbody=[2,3], mag=0.01,
        nac=phdb.nac
        )

### ASE calculator for forces
mode = 'force'
calc_force = apdb.get_calculator(mode,
        out_dirs['harm'][mode],
        phdb.get_kpoints(mode=mode).kpts[0]
        )

### calculate forces for harmonic IFCs
almcalc.calc_forces(1, calc_force, directory=out_dirs['harm'][mode])
almcalc.calc_harmonic_force_constants()

#### calculate band and DOS
almcalc.write_anphon_input(propt='band')
almcalc.run_anphon(propt='band', force=True)
#
almcalc.write_anphon_input(propt='dos')
almcalc.run_anphon(propt='dos')
almcalc.plot_bandos()

### Check negative frequency
if almcalc.frequency_range[0] < -0.001:
    print("")
    print(" Eigenvalues are negative. Stop the calculation.")
    print("")
    exit()

### calculate forces for cubic IFCs
mode = 'force'
almcalc.calc_forces(2, calc_force, out_dirs['cube'][mode],
        suggest_number_limit=1000, npattern_random=160)
almcalc.calc_cubic_force_constants()

### calculate kappa
almcalc.write_anphon_input(propt='kappa')
almcalc.run_anphon(propt='kappa')
almcalc.plot_kappa()

