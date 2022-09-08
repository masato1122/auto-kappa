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
      description='Automation software to calculate thermal conductivity with'\
              ' VASP and Alamode',
      author='Masato Ohnishi',
      author_email='ohnishi@photon.t.u-tokyo.ac.jp',
      packages=packages_pyscat,
      requires=['numpy', 'ALM', 'phonopy', 'ase', 'pymatgen', 'custodian'],
      url='https://github.com/masato1122/auto_alamode2.git',
      license='MIT',
      provides=['auto_alamode2']
      )

