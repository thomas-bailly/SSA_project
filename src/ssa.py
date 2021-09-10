import my_dssp as md

import sys
from datetime import datetime


def display_help():
    print("Usage: script.py pdb chain\n"
    "\nSecondary Structure Assignement\n"
    "\nArguments:\n  pdb\t\tpdb file\n  chain\t\ttargeted chain in pdb file")


def display_result(complete_structure):
    """Display the result.

    Parameter
    ---------
    complete_structure : list
        list with structure assignment for all residues
    """
    now = datetime.now()
    today = now.strftime("%Y/%m/%d %H:%M:%S")
    program = "Secondary Strcuture Assignement (ssa)"
    version = "1.0.0"
    author = "Thomas Bailly"
    texte = ("B : sense beta sheet\nb : anti-sense beta sheet\nC : random coil"
    "\nH : alpha helix\nG : 3,10 helix\nI : pi helix\nT : turn")
    delimiter = "="*50
    separetor = "-"*50
    input_pdb = sys.argv[1]
    chain = sys.argv[2]

    print("{}\nDate : {}\nProgram : {}\nVersion : {}\nAuthor : {}"
          "\nInput file : {}\nChaine : {}\n{}\n{}\n{}"
          .format(delimiter, today, program, version,author, input_pdb,
                  chain, separetor, texte, delimiter))

    for ligne in complete_structure:
        print("{:<6} {:<6} {:<6}".format(ligne[0], ligne[1], ligne[2]))

# Command line
try:
    if (sys.argv[1] == "help" or
    sys.argv[1]  == "-h"):
        display_help()
        sys.exit()
    elif sys.argv[1][-3:] != "pdb":
        print("Error: the file is not a pdb")
        display_help()
        sys.exit()
    elif sys.argv[1][-3:] == "pdb" and len(sys.argv) == 2:
        print("Error: You did not specify the chain")
        display_help()
        sys.exit()
    else:
        pass
except IndexError:
    display_help()
    sys.exit()

# Main Program

protein = md.pdb_parsing(sys.argv[1])
hbond = md.calc_hbond(protein)
hbond_clean = md.hbond_cleaner(hbond)
hbond_turn = md.turn_assignement(hbond_clean)
hbond_helix = md.helix_assignement(hbond_turn)
hbond_sheet = md.sheet_assignement(hbond_helix)
partial_structure = md.build_structure(hbond_sheet, protein)
complete_structure = md.structure_adjustement(partial_structure)
display_result(complete_structure)

