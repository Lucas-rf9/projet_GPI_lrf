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
    
    def add_atoms(self,atoms):
        self.list_atom.append(atoms)

#lecture du fichier pdb
with open("data/8D28.pdb", "r") as pdb_file:
    for line in pdb_file:
        if line.startswith("ATOM"):
            print(line.strip())

 
