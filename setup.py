from distutils.core import setup

packages_pyscat = [
        'auto_kappa',
        'auto_kappa.cui',
        'auto_kappa.io',
        'auto_kappa.plot',
        'auto_kappa.structure',
        ]

setup(name='auto_kappa',
      version='0.0',
      description='Auto-Alamode is a software to calculate thermal conductivity '\n
        'automatically with VASP and Alamode',
      author='Masato Ohnishi',
      author_email='ohnishi@photon.t.u-tokyo.ac.jp',
      packages=packages_pyscat,
      requires=[
          'numpy', 'ALM', 'phonopy', 'spglib', 
          'ase', 'pymatgen', 'custodian',
          'xmltodict', 'mkl', 'f90nml'
          ],
      url='https://github.com/masato1122/auto_kappa.git',
      license='MIT',
      provides=['auto_kappa']
      )

