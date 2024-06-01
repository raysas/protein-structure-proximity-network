#!/bin/env/python

import os, sys
from Bio import PDB
import networx as nx

os.chdir(os.path.expanduser('..'))

def get_PDB_structure(pdb_id:str)->PDB.Structure.Structure:
    '''
    Takes a pdb id and returns its sructure

    Starts by retriving the pdb file, saving it in data/{pdb_id}/{pdb_id}.pdb
    then reads the file and returns the structure

    param:
    ------
    - pdb_id: str, pdb id of the structure to be retrived

    return:
    -------
    - structure: Bio.PDB.Structure.Structure, structure of the pdb file
    '''
    pdb_list = PDB.PDBList()
    pdb_list.retrieve_pdb_file(pdb_id, pdir=f'data/{pdb_id}', file_format='pdb')
    os.rename(f'data/{pdb_id}/pdb{pdb_id.lower()}.ent', f'data/{pdb_id}/{pdb_id}.pdb')
        
    parser = PDB.PDBParser()
    structure = parser.get_structure(pdb_id, f'data/{pdb_id}/{pdb_id}.pdb')
    return structure

