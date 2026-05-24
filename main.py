#!/usr/bin/env python3
import sys
from rna_classes import Atom, Nucleotide

def convert_2d_rna_structure(file_name):
#lecture du fichier pdb
    dict_nucleotide = {}
    with open(file_name, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                column = line.strip().split()
                atom = Atom(column[2],float(column[6]),float(column[7]),float(column[8]))
                if column[5] in dict_nucleotide:
                    pass
                else:
                    dict_nucleotide[column[5]] = Nucleotide(column[3],column[5])
                dict_nucleotide[column[5]].add_atom(atom)

    #on trie les keys cad les numéros des nucléotides dans l'ordre
    keys = list(dict_nucleotide.keys())
    bracket_list = ["."] * len(keys)

    for i in range(len(keys)):
        for j in range(i + 4, len(keys)):
            nucleotide_A = dict_nucleotide[keys[i]]
            nucleotide_B = dict_nucleotide[keys[j]]
            number_bonds = nucleotide_A.count_hydrogen_bonds(nucleotide_B)
            #pour la paire G-C
            if number_bonds >= 3:
                if (nucleotide_A.name == "G" and nucleotide_B.name == "C") or (nucleotide_A.name == "C" and nucleotide_B.name == "G"):
                    bracket_list[i] = "("
                    bracket_list[j] = ")"
            #pour les paires à deux liaisons
            elif number_bonds >= 2:
                #pour la paire A-U
                if (nucleotide_A.name == "A" and nucleotide_B.name == "U") or (nucleotide_A.name == "U" and nucleotide_B.name == "A"):
                    bracket_list[i] = "("
                    bracket_list[j] = ")"
                #pour la paire G-U
                elif (nucleotide_A.name == "G" and nucleotide_B.name == "U") or (nucleotide_A.name == "U" and nucleotide_B.name == "G"):
                    bracket_list[i] = "("
                    bracket_list[j] = ")"

    sequence_list = []
    for i in range(len(bracket_list)):
        sequence_list.append(dict_nucleotide[keys[i]].name)

    print("".join(sequence_list))
    print("".join(bracket_list))

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Erreur : Veuillez fournir un fichier PDB en argument. Exemple : python projet_1.py 1EHZ.pdb")
        sys.exit(1)
    file_name = sys.argv[1]
    convert_2d_rna_structure(file_name)



