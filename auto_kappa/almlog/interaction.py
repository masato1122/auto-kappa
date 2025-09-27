# 
# interaction.py
# 
# This script parses "interaction" section in ALAMODE log files.
# 
# Created on September 16, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import re
from dataclasses import dataclass
from typing import List

@dataclass
class Atom:
    index: int
    element: str

@dataclass
class NeighborAtoms:
    serial: int
    distance: float
    num_neighbors: int
    neighbors: List[Atom]

@dataclass
class PairAtoms:
    center: Atom
    neighbors: List[NeighborAtoms]

def _parse_neighbor_line(line):
    
    blank = r"\s*"
    v_int = fr"{blank}(\d+){blank}"
    v_str = fr"{blank}(\S+){blank}"
    v_number = fr"{blank}([-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?){blank}"
    
    ##
    data_pattern = re.compile(
        fr"""
        {v_int}                   # index
        \s+{v_number}             # distance
        {blank}\({v_int}\){blank}  # number of neighbors
        -{blank}(.*)              # neighbor section
        """, re.VERBOSE)
    header_pattern = re.compile(
        fr"""
        {v_int}                 # serial
        \({v_str}\):            # element
        {data_pattern.pattern}  # data pattern
        """, re.VERBOSE)
    neighbor_pattern = re.compile(fr"{v_int}\({blank}{v_str}{blank}\)")
    
    match_head = header_pattern.match(line)
    if match_head:
        index1 = int(match_head.group(1))
        element1 = match_head.group(2)
        serial = int(match_head.group(3))
        distance = float(match_head.group(4))
        num_neighbors = int(match_head.group(5))
        neighbor_section = match_head.group(6)
        neighbors = neighbor_pattern.findall(neighbor_section)
        # print(f"{index1} ({element1:2s})  {distance:9.5f} ({num_neighbors:2d})", end=" ")
        # print(neighbors)
        
        center = PairAtoms(
            center=Atom(index=index1, element=element1),
            neighbors = [
                NeighborAtoms(
                    serial=serial,
                    distance=distance,
                    num_neighbors=num_neighbors,
                    neighbors=[Atom(index=int(n[0]), element=n[1]) for n in neighbors]
                    )
                ]
            )
        return center
    
    match_data = data_pattern.match(line)
    if match_data:
        serial = int(match_data.group(1))
        distance = float(match_data.group(2))
        num_neighbors = int(match_data.group(3))
        neighbor_section = match_data.group(4)
        neighbors = neighbor_pattern.findall(neighbor_section)
        # print(f"          {distance:9.5f} ({num_neighbors:2d})", end=" ")
        # print(neighbors)
        neighbor_atoms = NeighborAtoms(
            serial=serial,
            distance=distance,
            num_neighbors=num_neighbors,
            neighbors=[Atom(index=int(n[0]), element=n[1]) for n in neighbors]
        )
        return neighbor_atoms

def _get_list_of_neighbors(lines):
    
    info = []
    
    title = "List of neighboring atoms"
    reading_data = False
    for line in lines:
        
        if not line:
            continue
        line = line.strip()
        
        ## Title line
        if line.lower().startswith(title.lower()):
            reading_data = True
            continue
        
        if not reading_data:
            continue
        
        out = _parse_neighbor_line(line)
        if out is None:
            continue
        
        if hasattr(out, 'center'):
            info.append(out)
        else:
            info[-1].neighbors.append(out)
    
    return info
        
def _get_cutoff(lines):
    info = {}    
    pattern_title = r"\+\+\+ Cutoff Radii Matrix.*\+\+\+"
    print(" Not implemented yet.")


def read_interaction(lines):
    """ Read interaction information from log file lines (Interaction section)
    
    Args:
        lines (list): list of lines in the log file
    """
    info = {}
    
    # _get_cutoff(lines)
    
    info['pair_list'] = _get_list_of_neighbors(lines)
    
    return info

