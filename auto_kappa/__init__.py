from auto_kappa.version import __version__
from auto_kappa.cui.ak_log import set_logging

output_directories = {
        'relax': 'relax',
        'nac': 'nac',
        ### Harmonic FCs with finite displacement method
        'harm':{
            'suggest': 'harm/suggest',
            'force':   'harm/force',     ## FC2 with finite displacement
            'bandos':  'harm/bandos',
            'evec' :   'harm/evec',
            'pes'  :   'harm/pes',
            },
        ### cubic FCs
        'cube':{
            'suggest':  'cube/suggest',
            'force_fd': 'cube/force_fd',         ## FC3 with finite displacement
            'kappa_fd': 'cube/kappa_fd',
            'force_lasso': 'cube/force_lasso',   ## FC3 with random displacement
            'cv'   :       'cube/lasso',
            'lasso':       'cube/lasso',
            'kappa_lasso': 'cube/kappa_lasso',
            'gruneisen': 'cube/gruneisen',
            },
        ### high-order FCs using LASSO and fixed harmonic and cubic FCs
        'higher':{
            'suggest': 'higher/suggest',
            'force'  : 'higher/force',
            'cv'     : 'higher/lasso',
            'lasso'  : 'higher/lasso',
            'scph'   : 'higher/scph',
            'kappa_scph'    : 'higher/kappa_scph',
            'kappa_4ph'     : 'higher/kappa_4ph',
            'kappa_scph_4ph': 'higher/kappa_scph_4ph',
            },
        'result': 'result'
        }

output_files = {
        'harm_dfset' : 'DFSET.harm',
        'harm_xml'   : "FC2.xml",
        #
        'cube_fd_dfset'   : 'DFSET.cube_fd',
        'cube_fd_xml'     : "FC3_fd.xml",
        'cube_lasso_dfset': 'DFSET.cube_lasso',
        'cube_lasso_xml'  : "FC3_lasso.xml",
        #
        'higher_dfset': 'DFSET.high_lasso',
        'higher_xml'  : "FCs_high.xml",
        ##'lasso_dfset': 'DFSET.lasso',
        ##'lasso_xml'  : "FCs_lasso.xml",
        }

#
# According the manual for VASP:
# "If you have no a priori knowledge of your system, for instance, if you do not know whether your system is an insulator, 
# semiconductor or metal then always use Gaussian smearing ISMEAR=0 in combination with a small SIGMA=0.03-0.05."
# ISMEAR = 0 : Gaussian smearing
# SIGMA : width of smearing
#
default_vasp_parameters = {
        ### for relaxation
        'relax':{
            'prec': 'Accurate',
            'ibrion': 2,
            #'nsw': 500,
            'ediffg': -1e-6,
            'addgrid': True,
            },
        'relax-full':{    
            'isif': 3,
            },
        'relax-freeze':{
            'isif': 2,
            },
        ### force calculation
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
            'ismear': 0,      # Gaussian smearing
            'ialgo': 38,
            'lreal': False,
            'addgrid': True,
            'lwave': False,
            'lcharg': False,
            }
        }

default_amin_parameters = {"value": 0.01, "tol_length": 50, "num_of_errors": 1}

#default_lasso_parameters = {
#        'lmodel': 'enet',
#        'cv': 5,
#        'num_l1_alpha': 50
#        }

#default_lassobyscikit_parameters = {
#        'n_alphas': 100,
#        'eps': 1.e-7,
#        'n_splits': 5,
#        'standardize': True,
#        }

