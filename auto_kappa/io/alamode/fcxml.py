# -*- coding: utf-8 -*-
import math
import numpy as np
import pandas as pd
from optparse import OptionParser
import xml.etree.cElementTree as ET
import ase

from auto_kappa.structure.crystal import get_primitive_structure_spglib
import auto_kappa.units as unit

import logging
logger = logging.getLogger(__name__)

class FCXML():
    """ Parse XML file for force constants created by ALAMODE.
    Note that this class uses SI units (Angstrom, eV, ...) while atomic units 
    (Bohr, Ry, ...) are used in XML file of Alamode.
    
    How to Use
    =============
    >>> 
    >>> auto_kappa.io.alamode.fcxml import FCXML
    >>> filename = "Si_ANHARM.xml"
    >>> xml = FCXML(filename)
    >>> print(xml.get_fc_info(order=2, which="unique")) 
    >>> print(xml.get_force_constants(order=2))
    >>> print(xml.p2s_map)
    >>> print(xml.s2p_map)
    >>> print(xml.get_translations())
    >>> 
    """
    def __init__(self, filename):
        
        self._filename = filename
        self._root = None
        self._contents = None

        ###
        self._structure = None
        self._supercell = None
        self._primitive = None
        
        self._translations = None
        self._p2s_map = None
        self._s2p_map = None
        
        self._force_constants = {}   ## dictionary of ndarray
        self._fc_info = {"all": {}, "unique": {}}   ## dictionary of DataFrame
        
        if self.filename is not None:
            self._get_contents()
    
    @property
    def filename(self):
        return self._filename
    
    @property
    def root(self):
        if self._root is None:
            tree = ET.parse(self.filename)
            self._root = tree.getroot()
        return self._root
    
    def _get_contents(self):
        self._contents = _parse_alm_xml(self.root)
        return self._contents
    
    @property
    def contents(self):
        if self._contents is None:
            self._contents = self._get_contents()
        return self._contents
    
    @property
    def structure(self):
        if self._structure is None:
            info_struct = self.contents["Structure"]
            
            cell = info_struct["LatticeVector"]
            scaled_pos = info_struct["Position"]['position']
            pos = np.dot(scaled_pos, cell)
            
            self._structure = ase.Atoms(
                    cell=cell, pbc=True, 
                    symbols=info_struct["Position"]["element"],
                    positions=pos)
            
        return self._structure
    
    @property
    def supercell(self):
        if self._supercell is None:
            self._supercell = self.structure
        return self._supercell
    
    @property
    def primitive(self):
        if self._primitive is None:
            self._primitive = get_primitive_structure_spglib(self.supercell)
        return self._primitive
    
    @property
    def p2s_map(self):
        """ Return the atomic indices in the supercell (>= 1) corresponding to 
        the atomic indices in the primitive cell.
        
        Return
        --------
        
        array of integer : shape=(natoms_prim,)
        
        """
        if self._p2s_map is None:
            mapping = self.get_translations(inverse=True)
            self._p2s_map = mapping[0,:]
        return self._p2s_map
    
    @property
    def s2p_map(self):
        """ Return the atomic (site) indices (>= 1) in the primitive cell 
        corresponding to the atoms in the supercell.
        
        Return
        ---------
        
        array of integer : shape=(natoms_sc)
        
        """
        if self._s2p_map is None:
            mapping = self.get_translations(inverse=False)
            self._s2p_map = [each['atom'] for each in mapping]
        return self._s2p_map
     
    def get_translations(self, inverse=False):
        """ Get info of the supercell, primitive cell, and translational vector.
        The ``iat_sc``-th atom in the supercell corresponds to ``isite``-th 
        site in the ``icell``-th cell.
        
        Return
        ------
        
        info[iat_sc][key] : array of dictionary, shape=(natoms_sc)
            When inverse == False, return the cell index ("trans") and/or the 
            site index ("atom") of the given atom in the supercell ("iat_sc").
            
            ``iat_sc`` : atomic index in the supercell (0 to natoms_sc-1)
            ``key`` : "trans", "atom", or "translation"
            "trans" (>=1) is the cell index of primitive cell in the supercell,
            "atom" (>=1) in the site index, and
            "translation" (>=1) is the atomic index, which should be arranged in
            increasing order.
        
        info[icell][isite] : array, shape=(num_cells, natoms_prim)
            When inverse == True, return the atom index in the supercell (1 to
            natoms_sc) of the ``isite``-th site in the ``icell``-th cell.
            
            Atomic index in the supercell (1 to # of atoms in the supercell)
            ``icell`` : index of primitive cell in the supercell
            ``isite`` : site index in a primitive cell
        
        How to Use
        ------------
        >>> translations = self.get_translations()
        >>> iat_sc = 0
        >>> icell = translations[iat_sc]['tran'] - 1
        >>> isite = translations[iat_sc]['atom'] - 1
        
        >>> map_inv = self.get_translations(inverse=True)
        >>> iat = map_inv[icell,isite]
        """
        dump = self.contents["Symmetry"]["Translations"]
        
        natoms_sc = len(self.structure)
        
        ###
        if inverse == False:
            
            info = []
            for iat_sc in range(natoms_sc):
                info.append({})
                for key in dump:
                    info[iat_sc][key] = None
            
            for iat_sc in range(natoms_sc):
                icell = dump["translation"][iat_sc] - 1
                for key in dump:
                    info[icell][key] = int(dump[key][iat_sc])
        else:
            ncells = self.contents["Symmetry"]["NumberOfTranslations"]
            natoms_prim = len(self.primitive)
            info = np.zeros((ncells, natoms_prim), dtype=int)
            for iat_sc in range(natoms_sc):
                icell = dump["tran"][iat_sc] - 1
                isite = dump["atom"][iat_sc] - 1
                info[icell,isite] = int(dump["translation"][iat_sc])
        
        return info
    
    def get_force_constants(self, order=2):
        """ Get and return force constants with an ndarray 
        
        Returns
        ---------
        
        fcs : ndarray, shape=(natoms_primt, natoms_sc, ..., 3, 3, ...)
        
        To Do
        --------
        Check with NaCl, including two kinds of elements.

        Output force constants with Phonopy format and calcualte phonon
        properties with Phonopy to check the force constants.
        
        """
        ### Return values if they have been already extracted.
        if order in self._force_constants:
            if len(self._force_constants[order]) is not None:
                return self._force_constants[order]
        
        msg = "\n WARNING: get_force_constants is still under testing."
        msg += "\n Please use with carefull attention!!"
        logger.warning(msg)
        
        natoms_prim = len(self.primitive)
        natoms_sc = len(self.supercell)
        
        idx_dirs = {"x": 0, "y": 1, "z": 2}
        
        fc_info = self.get_fc_info(order=order, which="all")
        
        ### fc_info => matrix
        if order == 2:
            """ order 2 """        
            fcs = np.zeros((natoms_prim, natoms_sc, 3, 3))
            
            for ifc in range(len(fc_info)):
                
                iat1 = fc_info["atom1"].values[ifc] - 1
                dir1 = fc_info["direction1"].values[ifc] 
                idir1 = idx_dirs[dir1]
                
                iat2 = fc_info["atom2"].values[ifc] - 1
                dir2 = fc_info["direction2"].values[ifc]
                idir2 = idx_dirs[dir2]
                
                val = fc_info["FC"+str(order)].values[ifc]
                fcs[iat1, iat2, idir1, idir2] = val
            
            self._force_constants[order] = fcs
        
        elif order == 3:
            """ order 3 """ 
            fcs = np.zeros((natoms_prim, natoms_sc, natoms_sc, 3, 3, 3))
            
            for ifc in range(len(fc_info)):
                
                iat1 = fc_info["atom1"].values[ifc] - 1
                dir1 = fc_info["direction1"].values[ifc] 
                idir1 = idx_dirs[dir1]
                
                iat2 = fc_info["atom2"].values[ifc] - 1
                dir2 = fc_info["direction2"].values[ifc]
                idir2 = idx_dirs[dir2]
                
                iat3 = fc_info["atom3"].values[ifc] - 1
                dir3 = fc_info["direction3"].values[ifc]
                idir3 = idx_dirs[dir3]
                
                val = fc_info["FC"+str(order)].values[ifc]
                fcs[iat1, iat2, iat3, idir1, idir2, idir3] = val
            
            self._force_constants[order] = fcs
        
        elif order == 4:
            """ order 4 """
            fcs = np.zeros((
                natoms_prim, natoms_sc, natoms_sc, natoms_sc, 
                3, 3, 3, 3))
            
            for ifc in range(len(fc_info)):
                
                iat1 = fc_info["atom1"].values[ifc] - 1
                dir1 = fc_info["direction1"].values[ifc] 
                idir1 = idx_dirs[dir1]
                
                iat2 = fc_info["atom2"].values[ifc] - 1
                dir2 = fc_info["direction2"].values[ifc]
                idir2 = idx_dirs[dir2]
                
                iat3 = fc_info["atom3"].values[ifc] - 1
                dir3 = fc_info["direction3"].values[ifc]
                idir3 = idx_dirs[dir3]
                
                iat4 = fc_info["atom4"].values[ifc] - 1
                dir4 = fc_info["direction4"].values[ifc]
                idir4 = idx_dirs[dir4]
                
                val = fc_info["FC"+str(order)].values[ifc]
                fcs[iat1, iat2, iat3, iat4, idir1, idir2, idir3, idir4] = val
            
            self._force_constants[order] = fcs
        
        else:
            msg = "\n Error: order = %d is not yet supported." % order
            logger.error(msg)
            return None
        
        return self._force_constants[order]
    
    def get_fc_info(self, order=2, which="all"):
        """ Get and return info of force constants with DataFrame. 
        If ``which`` == "all", a matrix of FCs with the shape of 
        (natoms_prim, natoms_sc, ..., 3, 3, ...) is returned. 
        If ``which`` == "unique", an array of dictionaries is
        returned. Each element contains the indices of atoms and directions, 
        multiplicity, and the magnitude of a FC.
        
        Return (``which`` == "all")
        ----------------------------
        
        fcs : ndarray, float, shape=(natoms_prim, natoms_sc, ..., 3, 3, ...)
            if which == "all".
        
        Return (``which`` == "unique")
        --------------------------------
        
        fcs : array of dictionary, shape=(num_fcs,)
            Each dictionary contains "atom_info", "multiplicity", and "FC*".
            "multiplicity" and "FC*" contains the same info as that of input.
            
        fcs["atom_info"] : array of dictionary, shape=(order,)
            Each dictionary contains "atom" and "direction" keys, which denote the
            atom (1, 2, ...) and direction (1, 2, 3) indices.
        
        Format
        -------
        
        each_fc = self.contents["ForceConstants"]["HARMONIC"]
        each_fc["pair1"][0] : atom in the primitive cell
        each_fc["pair1"][1] : direction (1.x, 2.y, 3.z)
        each_fc["pair2"][0] : atom in the supercell 
        each_fc["pair2"][1] : direction (1.x, 2.y, 3.z)
        each_fc["pair2"][2] : cell index of the primitive cell in the supercell
        
        >>> print(each_fc)
        {'pair1': [1, 1], 'pair2': [5, 1, 6], 'FC2': -0.002526518981468988}
        
        """
        ### Return values if they have been already extracted.
        if which in self._fc_info:
            if order in self._fc_info[which]:
                if len(self._fc_info[which][order]) is not None:
                    return self._fc_info[which][order]
        
        ### get key for different orders
        if which == "all":
            if order == 2:
                key_order = "HARMONIC"
            else:
                key_order = "ANHARM" + str(order)
        elif which == "unique":
            if order == 2:
                key_order = "HarmonicUnique"
            elif order == 3:
                key_order = "CubicUnique"
            else:
                msg = "\n Warning: no %s FC data for order = %d" % (which, order)
                logger.warning(msg)
                return None
        
        ### get FC info including atomic and site indices, value, etc.
        dump = self.contents["ForceConstants"][key_order]
        
        key_fc = "FC%ds" % (order)
        
        nfcs = len(dump[key_fc])
        natoms_prim = len(self.primitive)
        natoms_sc = len(self.supercell)
        
        ### Get FCs as a matrix with the shape of (natoms_prim, natoms_sc)
        if which == "all":
            df = _extract_all_fc_info(
                    dump[key_fc], order=order,
                    natoms_prim=len(self.primitive), 
                    natoms_sc=len(self.supercell))
        elif which == "unique":
            df = _extract_unique_fc_info(dump[key_fc], order=order)
        else:
            msg = "\n Error: %s is not supported." % which
            logger.error(msg)
            return None
        
        ### add distance info
        all_distances = self.structure.get_all_distances(mic=True)
        distances = []
        for ifc in range(len(df)):
            idx_atoms = []
            
            """ 
            CAUTION!!
            ==========
            For iorder == 0, atom index may denote the index in the primitive 
            cell, but not in the supercell. In this case, this part should be 
            corrected!!
            """
            for iorder in range(order):
                idx_atoms.append(df["atom"+str(iorder+1)].values[ifc] - 1)
            
            ds_each = []
            for i1 in range(order):
                iat1 = idx_atoms[i1]
                for i2 in range(i1+1, order):
                    iat2 = idx_atoms[i2]
                ds_each.append(all_distances[iat1, iat2])
            
            ### get the maximum value
            distances.append(np.max(ds_each))
        
        df["distance[A]"] = distances
        
        ###
        self._fc_info[which][order] = df
        return self._fc_info[which][order]
    
    #def get_fcs_and_distances(self, order=None):
    #    """ Get FCs of the given order and return it. If order is 2 or 3, the
    #    unique FCs are used. 
    #    
    #    Return
    #    --------
    #    
    #    this DataFrame contains the similar info 
    #
    #    
    #    """
    #    
    #    if order is None:
    #        _order_is_not_given()
    #        return None 
    #    
    #    ### decide extracted data type
    #    if order <= 3:
    #        which = "unique"
    #    else:
    #        which = "all"
    #    
    #    ### key for each FC
    #    key_fc = "FC" + str(order)
    #    
    #    all_distances = self.structure.get_all_distances(mic=True)
    #    directions = ["x", "y", "z"]
    #    
    #    ### get force constants
    #    fcs = self.get_force_constants(order=order, which=which)
    #    
    #    if which == "unique":
    #        """ Get unique FCs and corresponding atomic distances """
    #        
    #        nfcs = len(fcs)
    #        
    #        ### prepare a dictionary
    #        info = {}
    #        for iorder in range(order):
    #            info["atom" + str(iorder+1)] = []
    #            info["direction" + str(iorder+1)] = []
    #        info["multiplicity"] = []
    #        info["distance[A]"] = []
    #        info[key_fc] = []
    #        
    #        for ifc in range(nfcs):
    #            
    #            ### atomic and direction indices
    #            idx_atoms = []
    #            for iorder in range(order):
    #                iatom = fcs[ifc]["atom_info"][iorder]["atom"]
    #                idir = fcs[ifc]["atom_info"][iorder]["direction"]
    #                idx_atoms.append(iatom - 1)
    #
    #                info["atom" + str(iorder+1)].append(iatom)
    #                info["direction" + str(iorder+1)].append(directions[idir-1])
    #            
    #            ### multiplicity
    #            info["multiplicity"].append(fcs[ifc]["multiplicity"])
    #
    #            ### force constant
    #            val = fcs[ifc][key_fc]
    #            info[key_fc].append(val)
    #            
    #            ###
    #            ds_each = []
    #            for i1 in range(order):
    #                iat1 = idx_atoms[i1]
    #                for i2 in range(i1+1, order):
    #                    iat2 = idx_atoms[i2]
    #                    ds_each.append(all_distances[iat1, iat2])
    #            
    #            ### compare with .fcs file
    #            #alpha_fc = unit.EvToRy / math.pow(unit.AToBohr, order)
    #            #print(np.max(ds_each) * unit.AToBohr, val * alpha_fc)
    #            
    #            info["distance[A]"].append(np.max(ds_each))
    #
    #        df = pd.DataFrame(info)
    #        return df
    #    
    #    elif which == "all":
    #
    #        ds_each = []
    #        for i1 in range(len(self.primitive)):
    #            iat1 = self.p2s_map[i1] - 1
    #            for i2 in range(len(self.supercell)):
    #                iat2 = i2
    #                d12 = all_distances[iat1, iat2]
    #                ds_each.append(d12)
    #                
    #                #print(iat1, iat2, d12, "\n", fcs[i1, i2])
    #        
    #        msg = "\n Error: order >= 4 is not supported yet."
    #        logger.error(msg)
    #        return None


def _prepare_dict_for_FCs(order, which):
    """ prepare a dictionary """
    info = {}
    for iorder in range(order):
        info["atom" + str(iorder+1)] = []
        info["direction" + str(iorder+1)] = []
        if iorder > 0 and which == "all":
            info["cell" + str(iorder+1)] = []
    if which == "unique":
        info["multiplicity"] = []
    info["FC" + str(order)] = []
    return info

def _extract_unique_fc_info(dump, order=None):
    """ Extract unique force constants 
    
    Args
    -------
    
    dump : array of dictionary, shape=(num_fcs,)
        Each dictionary contains "pairs", "multiplicity", and "FC*".
        
    dump["pairs"] : array, integer, shape=(order)
        The i-th element provides the info of atom and direction.
    
    dump["multiplicity"] : integer
    
    dump["FC*"] : float
        force constant (FC)
    
    Returns
    --------
    
    df : pandas.DataFrame
        contents are similar to the info in .fcs file
        
    """
    key_fc = "FC" + str(order)
    info = _prepare_dict_for_FCs(order, "unique")
    directions = ["x", "y", "z"]
    for each in dump:
        
        ### get index of atom and direction
        idx_atoms = []
        for iorder in range(order):
            ipair = each["pairs"][iorder]

            ### atomic index
            iatom = int(ipair / 3) + 1
            info["atom" + str(iorder+1)].append(iatom)
            
            ### direction index
            idir = ipair % 3 + 1
            info["direction" + str(iorder+1)].append(directions[idir-1])
        
        ### multiplicity and FC data
        info["multiplicity"].append(each["multiplicity"])
        info[key_fc].append(each[key_fc])
    
    ###
    df = pd.DataFrame(info)
    return df

def _extract_all_fc_info(dump, order=None, natoms_prim=None, natoms_sc=None):
    """ Make and return the matrix of force constants (FCs).
    
    Args
    ------
    dump : array of dictionary, shape=(num_fcs,)
        Each dictionary contains "pair1", "pair2", ..., "pair*" and "FC2", 
        which are info of a FC.
    
    dump["pair*"] : array, integer, shape=(2 or 3)
        Shape is 2 for ``order`` = 2 and 3 for ``order`` >= 3.
        dump["pair*"][0] (1 to ``natoms_sc``) and [1] (1, 2, 3) are the atom 
        and direction (1, 2, 3) indices, and dump["pair*"][2] (>=1) is the cell 
        index of the primitive cell in the supercell.
        
    Returns
    --------
    
    df : DataFrame
    
    """
    key_value = "FC" + str(order)
    nfcs = len(dump)
    directions = ["x", "y", "z"]
    info = _prepare_dict_for_FCs(order, "all")
    for each in dump:
        for iorder in range(order):
            
            key = "pair" + str(iorder+1)
            iatom = each[key][0]
            idir = each[key][1]

            info["atom"+str(iorder+1)].append(iatom)
            info["direction"+str(iorder+1)].append(directions[idir-1])
            
            if iorder > 0:
                icell = each[key][2]
                info["cell"+str(iorder+1)].append(icell)

        ###
        key = "FC" + str(order)
        info[key].append(each[key])
    
    df = pd.DataFrame(info)
    return df

def _order_is_not_given():
    """ """
    msg = "\n Error: order must be given."
    logger.error(msg)

def _parse_alm_xml(root):
    """ Parse alamode XML file and return info.
    In xml file, Bohr and Ry are used.
    """
    info = {}
    for child in root:
        tag = child.tag
        if tag == "ALM_version":
            contents = child.text
        elif tag == "Optimize":
            contents = _parse_optimize(child)
        elif tag == "Structure":
            contents = _parse_structure(child)
        elif tag == "Symmetry":
            contents = _parse_symmetry(child)
        elif tag == "ForceConstants":
            contents = _parse_forceconstants(child)
        else:
            _unknwon_tag(tag)
            continue
        info[tag] = contents
    return info

def _unknwon_tag(tag):
    msg = "\n Warning: found unknown tag in XML file, %s" % tag
    logger.warning(msg)

def _parse_optimize(root):
    """ Read "Optimize" part """
    info = {}
    for child in root:
        tag = child.tag
        if tag == "Constraint":
            info[tag] = int(child.text)
        elif "DFSET":
            info[tag] = child.text
        else:
            _unknown_tag(tag)
            continue
    return info

def _parse_structure(root):
    """ Read "Structure"part and make a structure """
    info = {}
    for child in root:
        tag = child.tag
        if "NumberOf" in tag:
            ### Number of atoms and elements
            contents = int(child.text)
        elif tag == "AtomicElements":
            ### chemical symbols
            contents = []
            for grand in child:
                content = grand.attrib
                content["symbol"] = grand.text
                contents.append(content)
        elif tag == "LatticeVector":
            ### lattice vector
            matrix = {}
            for i, grand in enumerate(child):
                matrix[grand.tag] = np.asarray(
                        [float(dd) * unit.BohrToA for dd in grand.text.split()])
            contents = np.zeros((3,3))
            for j in range(3):
                contents[j] = matrix["a"+str(j+1)]
        elif tag == "Periodicity":
            ### periodicity
            contents = [int(ii) for ii in child.text.split()]
        elif tag == "Position":
            ### scaled positions
            count = 0
            dump = {"index": [], "element": [], "position": []}
            for grand in child:
                key = "index"
                dump[key].append(int(grand.attrib[key]))
                key = "element"
                dump[key].append(grand.attrib[key])
                
                ### scaled position
                key = "position"
                dump["position"].append([
                    float(dd) for dd in grand.text.split()])
                count += 1
            
            idx_order = [ii - 1 for ii in dump['index']]
            index = [i for i in idx_order]
            elements = [dump["element"][i] for i in idx_order]
            positions = np.asarray([dump["position"][i] for i in idx_order])
            contents = {
                    "index": index, "element": elements, "position": positions}
        else:
            _unknown_tag(tag)
            continue
        info[tag] = contents 
    return info

def _parse_symmetry(root):
    """ Read "Symmetry" """
    info = {}
    for child in root:
        tag = child.tag
        if tag == "NumberOfTranslations":
            info[tag] = int(child.text)
        elif tag == "Translations":
            contents = {"tran": [], "atom": [], "translation": []}
            count = 0
            for grand in child:
                key = "tran"
                contents[key].append(int(grand.attrib[key]))
                key = "atom"
                contents[key].append(int(grand.attrib[key]))
                key = "translation"
                contents[key].append(int(grand.text))
                count += 1
            info[tag] = contents
        else:
            _unknown_tag(tag)
            continue
    return info

def _parse_forceconstants(root):
    """ Read "ForceConstant". Tag should be "HarmonicUnique", "CubicUnique",
    "HARMONIC", "ANHARMONIC3", etc. """
    info = {}
    for child in root:
        tag = child.tag
        info[tag] = _parse_fcs_each_order(child)
    return info

def _parse_fcs_each_order(root):
    """ Read "HarmonicUniqe", "CubicUniue", "HARMONIC", etc. 
    
    Variables
    ============
    
    fcs : array of dictionary
        number of elements is equal to the number of FCs
    
    """
    order = None
    fc_tag = None
    info = {}
    fcs = []
    for child in root:
        tag = child.tag
        if "NFC" in tag:
            info[tag] = int(child.text)
            order = int(tag.replace("NFC", ""))
            info["order"] = order
        elif tag == "Basis":
            info[tag] = child.text
        elif tag.startswith("FC"):
             
            attrib = child.attrib

            fcs.append({})
            for key in attrib:
                if key == "multiplicity":
                    fcs[-1][key] = int(attrib[key])
                elif "pair" in key:
                    fcs[-1][key] = [int(i) for i in attrib[key].split()]
                else:
                    #msg = "\n Warning: unknown attribute (%s)" % key
                    #logger.warning(msg)
                    fcs[-1][key] = attrib[key]
            
            ### read FC value
            fcs[-1][tag] = float(child.text)
            fc_tag = tag
    
    ### get order
    if order is None:
        order = int(fc_tag[-1])
    
    ### Ry/Bohr^{order} => eV/A^{order}
    frac = unit.RyToEv / math.pow(unit.BohrToA, order)
    
    for i in range(len(fcs)):
        fcs[i][fc_tag] *= frac
    
    ### add FCs to info dictionary
    info[fc_tag+"s"] = fcs
    return info

