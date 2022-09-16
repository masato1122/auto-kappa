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
        'lasso':{    
            'evec' : 'lasso/evec',
            'force': 'lasso/force',
            'cv'   : 'lasso/cv',
            'lasso': 'lasso/lasso',
            'kappa': 'lasso/kappa',
            },
        'result': 'result'
        }

#default_vaspcalc_parameters = {
#        'setups': 'recommended',
#        'mpirun': 'mpirun',
#        'vasp': 'vasp',
#        'xc': 'pbesol',
#        }

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
        'harm_xml': "FCs_harm.xml",
        'cube_dfset': 'DFSET.cube',
        'cube_xml': "FCs_cube.xml",
        'lasso_dfset': 'DFSET.lasso',
        'lasso_xml': "FCs_lasso.xml",
        }

default_lasso_parameters = {
        'linear_model': 2,
        'cross_validation': 5,
        'num_l1_alpha': 50
        }

#default_lassobyscikit_parameters = {
#        'n_alphas': 100,
#        'eps': 1.e-7,
#        'n_splits': 5,
#        'standardize': True,
#        }

