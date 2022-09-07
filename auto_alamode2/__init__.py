# -*- coding: utf-8 -*-

output_directories = {
        'relax': 'relax',
        'nac': 'nac',
        'harm':{
            'force': 'harm/force',
            'bandos': 'harm/bandos',
            },
        'cube':{
            'force': 'cube/force',
            'kappa': 'cube/kappa',
            },
        'result': 'result'
        }

calculator_parameters = {
        'setups': 'recommended',
        'mpirun': 'mpirun',
        'vasp': 'vasp',
        'xc': 'pbesol',
        }

default_vasp_parameters = {
        'relax':{
            'prec': 'Accurate',
            'ibrion': 2,
            'nsw': 100,
            'isif': 3,
            'ediffg': -0.0001,
            'lreal': False,
            'addgrid': True,
            },
        'force':{
            'ibrion': -1,
            'prec': 'Accurate',
            },
        'nac':{
            'ibrion': -1,
            'lepsilon': True,
            'prec': 'Accurate',
            },
        'md':{
            'ibrion': 0,
            'potim': 1.0,
            'tebeg': 500.,
            'teend': 500.,
            'nsw':  3000,
            },
        'shared':{
            'nelmin': 5,
            'sigma': 0.02,
            'ediff': 1e-8,
            'ismear': 0,
            'ialgo': 38,
            'addgrid': True,
            'lwave': False,
            'lcharg': False,
            }
        }

output_files = {
        'harm_dfset': 'DFSET.harm',
        'cube_dfset': 'DFSET.cube',
        'harm_xml': "FCs_harm.xml",
        'cube_xml': "FCs_cube.xml",
        }

