import sys
import numpy as np
from scipy.spatial import distance
from Bio.PDB import PDBParser


def pdb_parsing(pdb):
    """Parser of pdb file.

    Searches for the atoms of the residues composing the targeted chain
    and store them.

    Parameter
    ---------
    pdb : str
        name of pdb file

    Return
    ------
    chain : Bio.PDB.chain.chain
        protein chain
    """
    parser = PDBParser()
    prot_id = pdb[0:4]
    prot_file = pdb
    structure = parser.get_structure(prot_id, prot_file)
    model = structure[0]
    chain = model[sys.argv[2]]

    return(chain)


def calc_hbond(chain):
    """calculates the energies of hydrogen bonds.

    calculates for each residue i its binding energy with the residue i+n.
    With n from 3 to length of the chain.
    The hydrogen bonds take place between the C=0 of residue i 
    and the N-H of residue j.
    F : dimensional factor = 332 kcal/mol
    q1 = partial charge of the oxygen of the carbonyl group = -0.42e
    q2 = partial charge of the hydrogen of the amide group = +0.20e

    Parameter
    ---------
    chain : Bio.PDB.chain.chain
        protein chain

    Return
    ------
    hbond : list
        list of hydrogen bonds
    """

    hbond = []
    q1q2 = 0.084
    F = 332

    for i in range(1, len(chain)-4):
        for j in range(i+3, len(chain)):
            try:
                # The amide group of proline cannot form a hydrogen bond.
                if (chain[j].resname == 'PRO'):
                    continue
                
                d = (q1q2*((1/distance.euclidean(chain[i]["O"],chain[j]["N"]))
                        +(1/distance.euclidean(chain[i]["C"],chain[j]["H"]))
                        -(1/distance.euclidean(chain[i]["O"],chain[j]["H"]))
                        -(1/distance.euclidean(chain[i]["C"],chain[j]["N"])))*F)
                # Threshold from which it's considered that there is a
                # hydrogen bond.
                if (d < (-0.5)):
                    hbond.append([i, j, d, 1])
            # If the chain in the pdb file does not start at residue 1.  
            except KeyError:
                pass

    return(hbond)


def hbond_cleaner(hbond):
    """Correct the list of hydrogen bonds.

    Removes possible duplicates.

    Paramater
    ---------
    hbond :  list
        list of hydrogen bonds

    Return
    ------
    hbond_clean : list
        cleaned list of hydrogen bonds
    """
    hbond_clean = []

    for i in range(0, len(hbond)):
        if (hbond[i][0] == hbond[i-1][0]):
            if (hbond[i][2] < hbond[i-1][2]):
                hbond[i-1][3] = 0
            else:
                hbond[i][3] = 0

    for bond in hbond:
        if bond[3] == 1:
            hbond_clean.append(bond)

    return(hbond_clean)


def turn_assignement(hbond_clean):
    """Turn assignement.

    Recognition of turn patterns.

    Paramter
    --------
    hbond_clean : list
        cleaned list of hydrogen bonds

    Return
    ------
    hbond_clean : list
        cleaned list with turn assignement
    """

    for i in range(0, len(hbond_clean)):
        # in [3,4,5] for 3-turn, 4-turn or 5-turn pattern.
        if hbond_clean[i][1] - hbond_clean[i][0] in [3, 4, 5]:
            hbond_clean[i][3] = "T"
        else:
            hbond_clean[i][3] = "C"

    return(hbond_clean)


def helix_assignement(hbond_turn):
    """Helix assignement.

    Recognition of alpha Helix, 3,10 Helix, pi Helix.

    Parameter
    ---------
    hbond_turn : list
        cleaned list with turn assignement

    Return
    ------
    hbond_turn : list
        cleaned list with turn and helix assignement
    """

    for i in range(0, len(hbond_turn)-1):
        if ((hbond_turn[i+1][0] - hbond_turn[i][0] == 1 and
             hbond_turn[i][3] == hbond_turn[i+1][3]) or
           (hbond_turn[i][0] - hbond_turn[i-1][0] == 1 and
            hbond_turn[i][3] == hbond_turn[i-1][3])):
            if (hbond_turn[i][1] - hbond_turn[i][0] == 4):
                        hbond_turn[i][3] = "H"
                        hbond_turn[i+1][3] = "H"
            elif (hbond_turn[i][1] - hbond_turn[i][0] == 3):
                        hbond_turn[i][3] = "G"
                        hbond_turn[i+1][3] = "G"
            elif (hbond_turn[i][1] - hbond_turn[i][0] == 5):
                        hbond_turn[i][3] = "I"
                        hbond_turn[i][3] = "I"

    return(hbond_turn)


def sheet_assignement(hbond_helix):
    """Beta sheet assignement.

    Recognition of sense and anti-sense sheet.

    Parameter
    ---------
    hbond_helix : list
        cleaned list with turn and helix assignement

    Return
    ------
    hbond_helix : list
        cleaned list with turn, helix and beta sheet assignement
    """

    for i in range(0, len(hbond_helix)-1):
        if (hbond_helix[i+1][0] - hbond_helix[i][0] == 2 and
            hbond_helix[i+1][1] - hbond_helix[i][1] == 2):
            hbond_helix[i][3] = "B"
            hbond_helix[i+1][3] = "B"

        elif (hbond_helix[i+1][0] - hbond_helix[i][0] == 2 and
              hbond_helix[i+1][1] - hbond_helix[i][1] == -2):
            hbond_helix[i][3] = "b"
            hbond_helix[i+1][3] = "b"

    return(hbond_helix)


def build_structure(hbond_sheet, chain):
    """Building of a list of residues and structures.

    Copy the structure assignment of the residues present in
    the hydrogen bonds list.

    Parameters
    ----------
    hbond_sheet : list
        cleaned list with turn, helix and beta sheet assignement
    chain : Bio.PDB.chain.chain
        protein chain 

    Return
    ------
    partial_structure : list
        list containing a partial assignment of the structure

    """

    partial_structure = []

    for i in range(1, len(chain)+1):
        try:
            partial_structure.append([i, chain[i].resname, "C"])
        except KeyError:
            pass

    for ligne in hbond_sheet:
        if ligne[3] == "T":

            for res in partial_structure:
                if ligne[0] == res[0]:
                    res[2] = "T"
                if ligne[1] == res[0]:
                    res[2] = "T"

    for ligne in hbond_sheet:
        if ligne[3] == "H":
            for res in partial_structure:
                if ligne[0] == res[0]:
                    res[2] = "H"
                if ligne[1] == res[0]:
                    res[2] = "H"

    for ligne in hbond_sheet:
        if ligne[3] == "G":
            for res in partial_structure:
                if ligne[0] == res[0]:
                    res[2] = "G"
                if ligne[1] == res[0]:
                    res[2] = "G"

    for ligne in hbond_sheet:
        if ligne[3] == "I":
            for res in partial_structure:
                if ligne[0] == res[0]:
                    res[2] = "I"
                if ligne[1] == res[0]:
                    res[2] = "I"

    for ligne in hbond_sheet:
        if ligne[3] == "B":
            for res in partial_structure:
                if ligne[0] == res[0]:
                    res[2] = "B"
                if ligne[1] == res[0]:
                    res[2] = "B"

    for ligne in hbond_sheet:
        if ligne[3] == "b":
            for res in partial_structure:
                if ligne[0] == res[0]:
                    res[2] = "b"
                if ligne[1] == res[0]:
                    res[2] = "b"

    return(partial_structure)


def structure_adjustement(partial_structure):
    """Adjustment of the structure assignment.

    Structure assignment for residues not in the hydrogen bonding list.

    Parameter
    ---------
    partial_structure : list
        list containing a partial assignment of the structure

    Return
    ------
    partial_structure : list
        updated list with structure assignment for all residues
    """

    for i in range(1, len(partial_structure)-1):
        if (partial_structure[i-1][2] == "B" and
            partial_structure[i+1][2] == "B"):
            partial_structure[i][2] = "B"

    for i in range(1, len(partial_structure)-1):
        if (partial_structure[i-1][2] == "b" and
            partial_structure[i+1][2] == "b"):
            partial_structure[i][2] = "b"

    for i in range(1, len(partial_structure)-1):
        if (partial_structure[i-1][2] == "H" and
            partial_structure[i+1][2] == "H"):
            partial_structure[i][2] = "H"

    for i in range(1, len(partial_structure)-1):
        if (partial_structure[i-1][2] == "G" and
            partial_structure[i+1][2] == "G"):
            partial_structure[i][2] = "G"

    for i in range(1, len(partial_structure)-1):
        if (partial_structure[i-1][2] == "I" and
            partial_structure[i+1][2] == "I"):
            partial_structure[i][2] = "I"

    for i in range(1, len(partial_structure)):
        if (partial_structure[i][2] == "T" and
            partial_structure[i+3][2] == "T"):
            partial_structure[i+1][2] = "T"
            partial_structure[i+2][2] = "T"

        if (partial_structure[i][2] == "T" and
            partial_structure[i+4][2] == "T"):
            partial_structure[i+1][2] = "T"
            partial_structure[i+2][2] = "T"
            partial_structure[i+3][2] = "T"

        if (partial_structure[i][2] == "T" and
            partial_structure[i+5][2] == "T"):
            partial_structure[i+1][2] = "T"
            partial_structure[i+2][2] = "T"
            partial_structure[i+3][2] = "T"
            partial_structure[i+4][2] = "T"

    return(partial_structure)
