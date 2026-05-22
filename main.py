#!/usr/bin/env python3

#création des différentes classes
class Atom:
    def __init__(self,name,coord_x,coord_y,coord_z):
        self.name = name
        self.x = coord_x
        self.y = coord_y
        self.z = coord_z

class Nucleotide:
    def __init__(self,name,position):
        self.name = name
        self.position = position
        self.list_atom = []
    
    def add_atom(self,atoms):
        self.list_atom.append(atoms)

#lecture du fichier pdb
dict_nucleotide = []

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

 
