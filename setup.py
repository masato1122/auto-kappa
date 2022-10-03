#from distutils.core import setup

import pathlib
from setuptools import setup

def _get_version():
    
    filename = "./auto_kappa/version.py"
    
    try:
        lines = open(filename, 'r').readlines()
        for line in lines:
            if '__version__' in line:
                line = line.translate(str.maketrans(
                    {"\"": " ", "\'": " ", "=": " "}
                    ))
                break
        data = line.split()
        return data[-1]
    
    except Exception:
        return "x.x"

def main(build_dir):
    
    packages_autokappa = [
            'auto_kappa',
            'auto_kappa.alamode',
            'auto_kappa.cui',
            'auto_kappa.io',
            'auto_kappa.math',
            'auto_kappa.plot',
            'auto_kappa.structure',
            ]

    scripts_autokappa = [
            'scripts/akrun'
            ]
    
    version = _get_version()

    
    setup(
            name='auto_kappa',
            version=version,
            description='auto-kappa module',
            author='Masato Ohnishi',
            author_email='ohnishi@photon.t.u-tokyo.ac.jp',
            packages=packages_autokappa,
            install_requires=[
              'numpy', 'phonopy', 'spglib', 'ase', 'pymatgen', 'custodian',
              'xmltodict', 'mkl', 'f90nml', 'PyYAML',
              ],
            scripts=scripts_autokappa,
            url='https://github.com/masato1122/auto_kappa.git',
            license='MIT',
            provides=['auto_kappa'],
            )
 

if __name__ == "__main__":
    
    build_dir = pathlib.Path.cwd() / "_build"
    
    main(build_dir)

