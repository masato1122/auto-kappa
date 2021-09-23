
__all__ = ['alm_variables', 'str_vars', 'int_vars', 'double_vars',
        'array_str_vars', 'array_int_vars', 'array_double_vars', 
        'hyphen_int_vars']

alm_variables = {
        'general':{
            'hessian': 0,
            'fcsym_basis': None,
            'kd': None,
            'magmom': None,
            'mode': None,
            'nat': None,
            'nkd': None,
            'nmaxsave': None,
            'noncollinear': None,
            'periodic': None,
            'prefix': None,
            'printsym': None,
            'tolerance': 1e-3,
            },
        'interaction':{
            'nbody': None,
            'norder': None,
            },
        'cutoff':{
            None       ## {'Si-Si': [None, 10.0]}
            },
        'cell':{
            None
            },
        'position':{
            None
            },
        'optimize':{
            'lmodel': None,
            'dfset': None,
            'ndata': None,
            'nstart': None,
            'nend': None,
            'skip': None,
            'iconst': None,
            'rotaxis': None,
            'fc2xml': None,
            'fc3xml': None,
            'sparse': None,
            'sparsesolver': None,
            'maxiter': None,
            'conv_tol': None,
            'l1_ratio': None,
            'l1_alpha': None,
            'cv': None,
            'dfset_cv': None,
            'ndata_cv': None,
            'nstart_cv': None,
            'nend_cv': None,
            'cv_minalpha': None,
            'cv_maxalpha': None,
            'cv_nalpha': None,
            'debias_ols': None,
            'enet_dnorm': None,
            'solution_path': None,
            'standardize': None,
            }
        }

str_vars = [
        'prefix', 'mode', 'fcsym_basis',  'lmodel', 'dfset', 'rotazis', 
        'fc2xml', 'fc3xml', 'sparsesolver'
        ]

int_vars = [
        'nat', 'nkd', 'printsym', 'noncollinear', 'nmaxsave', 'hessian',
        'norder', 'ndata', 'nstart', 'nend', 'iconst', 'sparse', 'maxiter',
        'cv', 'ndata_cv', 'nstart_cv', 'nend_cv', 'cv_nalpha', 'standardize', 
        'solution_path', 'debias_ols']

double_vars = [
        'tolerance', 'conv_tol', 'l1_ratio', 'l1_alpha', 'cv_minalpha',
        'cv_maxalpha', 'enet_dnorm'
        ]

array_str_vars = ['kd']

array_int_vars = ['periodic', 'nbody']

array_double_vars = ['magmom']

hyphen_int_vars = ['skip']

