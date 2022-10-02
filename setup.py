from distutils.core import setup

packages_pyscat = [
        'auto_kappa',
        'auto_kappa.alamode',
        'auto_kappa.cui',
        'auto_kappa.io',
        'auto_kappa.math',
        'auto_kappa.plot',
        'auto_kappa.structure',
        ]

setup(name='auto_kappa',
      version='0.1',
      description='Auto-KAPPA is a software to calculate thermal conductivity '\
        'automatically with VASP and Alamode',
      author='Masato Ohnishi',
      author_email='ohnishi@photon.t.u-tokyo.ac.jp',
      packages=packages_pyscat,
      setup_requires=[
          'numpy', 'ALM', 'phonopy', 'spglib', 
          'ase', 'pymatgen', 'custodian',
          'xmltodict', 'mkl', 'f90nml'
          ],
      url='https://github.com/masato1122/auto_kappa.git',
      license='MIT',
      provides=['auto_kappa']
      )

