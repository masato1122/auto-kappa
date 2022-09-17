from distutils.core import setup

packages_pyscat = [
        'auto_alamode2',
        'auto_alamode2.cui',
        'auto_alamode2.io',
        'auto_alamode2.plot',
        'auto_alamode2.structure',
        ]

setup(name='auto_alamode2',
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
      url='https://github.com/masato1122/auto_alamode2.git',
      license='MIT',
      provides=['auto_alamode2']
      )

