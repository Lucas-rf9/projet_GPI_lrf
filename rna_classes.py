#création des différentes classes
class Atom:
    """
    Représente un atome avec ses coordonnées spatiales 3D issues du fichier PDB.
    """
    def __init__(self,name,coord_x,coord_y,coord_z):
        self.name = name
        self.x = coord_x
        self.y = coord_y
        self.z = coord_z

    def calculate_distance(self, other_atom):
        """Calcule et retourne la distance euclidienne entre cet atome et un autre."""
        result = ((other_atom.x - self.x)**2 + (other_atom.y - self.y)**2 + (other_atom.z - self.z)**2)**0.5
        return result
        
class Nucleotide:
    """
    Représente un nucléotide, contenant une liste de ses atomes constitutifs.
    """
    def __init__(self,name,position):
        self.name = name
        self.position = position
        self.list_atom = []
    
    def add_atom(self,atoms):
        """Ajoute un objet Atom à la liste du nucléotide."""
        self.list_atom.append(atoms)

    def count_hydrogen_bonds(self, other_nucleotide):
        """
        Parcourt les atomes de deux nucléotides pour compter les liaisons hydrogènes
        potentielles (distance comprise dans la fenêtre spatiale de 2.5 à 3.5 Å).
        """
        number_hydrogene_bonds = 0
        for atom_a in self.list_atom:
            for atom_b in other_nucleotide.list_atom:
                result = atom_a.calculate_distance(atom_b)
                if 2.5 <= result <= 3.5:
                    number_hydrogene_bonds += 1
        return number_hydrogene_bonds
            