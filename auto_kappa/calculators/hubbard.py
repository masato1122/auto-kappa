# -*- coding: utf-8 -*-
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

def get_hubbard_params(atoms, ldautype=2):
    """
    Generate Hubbard U parameters for VASP from an ASE Atoms object.
    
    Parameters:
    - atoms: ASE Atoms object representing the structure.
    - hubbard_u_values: Dictionary mapping chemical symbols to (U, L, J) values. Example:
        {
            "Fe": {"U": 4.0, "L": 2, "J": 0.0},
            "O": {"U": 0.0, "L": -1, "J": 0.0}
        }
    - default_magmom: Default magnetic moment to assign to atoms without specific values (float).
    
    Returns:
    - dict: Dictionary with LDAUL, LDAUU, LDAUJ, and MAGMOM parameters.
    """
    from auto_kappa.calculators.vasp import get_default_params_pymatgen
    params_pmg = get_default_params_pymatgen()
    
    #unique_elements = sorted(set(atoms.get_chemical_symbols()))
    
    ldaul = []
    ldauu = []
    ldauj = []
    magmom = [0.0] * len(atoms)
    
    symbols = atoms.get_chemical_symbols()
    for i, element in enumerate(symbols):
        
        print(element)
        
        ## LDAUL
        ldaul.append(get_suggested_ldaul(element))
        
        ldauu_each = 0.
        ldauj_each = 0.
        for fo_element in ['F', 'O']:
            if element in fo_element:
                
                ## LDAUU
                if element in params_pmg['LDAUU'][fo_element]:
                    ldauu_each = params_pmg['LDAUU'][fo_element][element]
                
                ## LDAUJ
                if element in params_pmg['LDAUJ'][fo_element]:
                    ldauj_each = params_pmg['LDAUJ'][fo_element][element]
        
        ## MAGMOM
        if element in params_pmg['MAGMOM']:
            magmom[i] = params_pmg['MAGMOM'][element]
        
        ldauu.append(ldauu_each)
        ldauj.append(ldauj_each)
        
    atoms.set_initial_magnetic_moments(magmom)
    
    hubbard_params = {
            "ldau": True,
            "ldautype": ldautype,
            "ldaul": " ".join(map(str, ldaul)),
            "ldauu": " ".join(map(str, ldauu)),
            "ldauj": " ".join(map(str, ldauj)),
            "ispin": 2,
            "magmom": " ".join(map(str, magmom)),
            }
    
    print(magmom)
    exit()
    return hubbard_params

def get_suggested_ldaul(element):
    """
    Get suggestd value for LDAUL
    
    Parameters:
    - element (str): The chemical symbol of the element (e.g., "Fe", "Ce", "O").
    
    Returns:
    - str: "Transition Metal", "Rare Earth", or "Other".
    """
    # Define the periodic table groups for transition metals and rare earths
    if element in transition_metals:
        return 2
    elif element in rare_earths:
        return 3
    else:
        return -1
