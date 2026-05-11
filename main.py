#!/usr/bin/env python3


class Protein:
    def __init__(self,name):
        self.nom = name

    pass 

protein1 = Protein('hemoglobin')
protein2 = Protein('cytochrom c')
print(protein2.nom)