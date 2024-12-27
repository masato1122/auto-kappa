#
# hubbard.py
#
# This script helps to make VASP parameters for Hubbard U correction.
#
# Copyright (c) 2024 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
from auto_kappa.calculators import transition_metals, rare_earths

def use_hubbard(symbols):
    """
    Args
    ------
    symbols : list of str
        chemical symbols
    """
    # fluoride or oxide
    if 'F' not in symbols and 'O' not in symbols:
        return False
    
    # transition metals or rare earths
    flag_hubbard = False
    for metals in [transition_metals, rare_earths]:
        if bool(set(symbols) & set(metals)):
            flag_hubbard = True
    
    return flag_hubbard


# ind = np.argsort(atoms.symbols)
# symbols = atoms.symbols[ind]
# coord = coord[ind]
# if constraints_present:
#     sflags = sflags[ind]


def get_hubbard_params(
    atoms, default_params=None, ldautype=2, default_magmom=0.0, sort=True):
    """ 
    Generate Hubbard U parameters for VASP from an ASE Atoms object.
    
    Args
    -----
    atoms : Atoms object
    
    default_params : dict
        default Hubbard U parameters
    
    default_magmom : float
        Default magnetic moment to assign to atoms without specific values (float).
    
    Returns
    --------
    dict : Dictionary with LDAUL, LDAUU, LDAUJ, and MAGMOM parameters.
    """
    for key in ['LDAUU', 'LDAUJ', 'MAGMOM']:
        if key not in default_params:
            raise ValueError(f' "default_params" is missing required key: {key}')
    
    ## sort with alphabetical order
    if sort:
        ind = np.argsort(atoms.symbols)
        symbols = atoms.symbols[ind]
    else:
        symbols = atoms.symbols
    
    unique_elements = sorted(set(symbols), key=symbols.index)
    
    ## Get LDAUL, LDAUU, LDAUJ
    ldaul = []
    ldauu = []
    ldauj = []
    
    for element in unique_elements:
        
        ## LDAUL
        ldaul.append(get_suggested_ldaul(element))
        
        ldauu_each = 0.
        ldauj_each = 0.
        for fo_element in ['F', 'O']:
            if element in fo_element:
                
                ## LDAUU
                if element in default_params['LDAUU'][fo_element]:
                    ldauu_each = default_params['LDAUU'][fo_element][element]
                
                ## LDAUJ
                if element in default_params['LDAUJ'][fo_element]:
                    ldauj_each = default_params['LDAUJ'][fo_element][element]
        
        ldauu.append(ldauu_each)
        ldauj.append(ldauj_each)
    
    ## MAGMOM
    magmom = []
    for element in symbols:
        magmom.append(default_params.get("MAGMOM", {}).get(element, default_magmom))
    
    atoms.set_initial_magnetic_moments(magmom)
    
    hubbard_params = {
            "ldau": True,
            "ldautype": ldautype,
            "ldaul": " ".join(map(str, ldaul)),
            "ldauu": " ".join(map(str, ldauu)),
            "ldauj": " ".join(map(str, ldauj)),
            "ispin": 2,
            "magmom": magmom,
            }
    
    return hubbard_params

def get_suggested_ldaul(element):
    """
    Get suggestd value for LDAUL
    
    Args
    -----
    element : str
        The chemical symbol of the element (e.g., "Fe", "Ce", "O").
    
    Returns
    --------
    int : suggestd LDAUL parameter
    """
    # Define the periodic table groups for transition metals and rare earths
    if element in transition_metals:
        return 2
    elif element in rare_earths:
        return 3
    else:
        return -1
