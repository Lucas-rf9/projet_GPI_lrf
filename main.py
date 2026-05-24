#!/usr/bin/env python3

#création des différentes classes
class Atom:
    def __init__(self,name,coord_x,coord_y,coord_z):
        self.name = name
        self.x = coord_x
        self.y = coord_y
        self.z = coord_z

    def calculate_distance(self, other_atom):
        result = ((other_atom.x - self.x)**2 + (other_atom.y - self.y)**2 + (other_atom.z - self.z)**2)**0.5
        return result
        
class Nucleotide:
    def __init__(self,name,position):
        self.name = name
        self.position = position
        self.list_atom = []
    
    def add_atom(self,atoms):
        self.list_atom.append(atoms)

    def count_hydrogen_bonds(self, other_nucleotide):
        number_hydrogene_bonds = 0
        for atom_a in self.list_atom:
            for atom_b in other_nucleotide.list_atom:
                result = atom_a.calculate_distance(atom_b)
                if 2.5 <= result <= 3.5:
                    number_hydrogene_bonds += 1
        return number_hydrogene_bonds
            
#lecture du fichier pdb
dict_nucleotide = {}

with open("data/8D28.pdb", "r") as pdb_file:
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

for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
        nucleotide_A = dict_nucleotide[keys[i]]
        nucleotide_B = dict_nucleotide[keys[j]]
        number_bonds = nucleotide_A.count_hydrogen_bonds(nucleotide_B)
        #pour la paire G-C
        if number_bonds == 3:
            if (nucleotide_A.name == "G" and nucleotide_B.name == "C") or (nucleotide_A.name == "C" and nucleotide_B.name == "G"):
                print(nucleotide_A.position,"pb trouvée :)")




