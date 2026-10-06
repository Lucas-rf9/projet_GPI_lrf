"""Modélisation des atomes et nucléotides d'un ARN issus d'un fichier PDB."""

# Fenêtre de distance (en Å) entre un atome donneur et un atome accepteur
# pour considérer qu'une liaison hydrogène est formée.
HBOND_MIN_DISTANCE = 2.5
HBOND_MAX_DISTANCE = 3.5

# Couples d'atomes (donneur/accepteur) impliqués dans les liaisons hydrogènes
# de chaque paire de bases canonique, du côté Watson-Crick.
CANONICAL_HBONDS = {
    ("G", "C"): [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],
    ("A", "U"): [("N1", "N3"), ("N6", "O4")],
    ("G", "U"): [("N1", "O2"), ("O6", "N3")],
}

# Nombre minimal de liaisons hydrogènes pour valider chaque type de paire.
MIN_HBONDS = {
    ("G", "C"): 3,
    ("A", "U"): 2,
    ("G", "U"): 2,
}

# Nucléotides modifiés fréquents (notamment dans les ARNt) et leur base parente.
MODIFIED_NUCLEOTIDES = {
    "1MA": "A", "2MA": "A", "6IA": "A", "MIA": "A", "T6A": "A",
    "5MC": "C", "OMC": "C", "4OC": "C",
    "1MG": "G", "2MG": "G", "7MG": "G", "M2G": "G", "OMG": "G", "YYG": "G", "QUO": "G",
    "5MU": "U", "H2U": "U", "PSU": "U", "4SU": "U", "OMU": "U",
}

STANDARD_NUCLEOTIDES = {"A": "A", "C": "C", "G": "G", "U": "U",
                        "RA": "A", "RC": "C", "RG": "G", "RU": "U"}


def nucleotide_base(residue_name):
    """Retourne la base parente (A, C, G, U) d'un résidu PDB, ou None si ce n'est pas un nucléotide d'ARN."""
    residue_name = residue_name.strip().upper()
    return STANDARD_NUCLEOTIDES.get(residue_name) or MODIFIED_NUCLEOTIDES.get(residue_name)


class Atom:
    """
    Représente un atome avec ses coordonnées spatiales 3D issues du fichier PDB.
    """
    def __init__(self, name, coord_x, coord_y, coord_z):
        self.name = name
        self.x = coord_x
        self.y = coord_y
        self.z = coord_z

    def calculate_distance(self, other_atom):
        """Calcule et retourne la distance euclidienne entre cet atome et un autre."""
        return ((other_atom.x - self.x)**2 + (other_atom.y - self.y)**2 + (other_atom.z - self.z)**2)**0.5


class Nucleotide:
    """
    Représente un nucléotide, identifié par sa chaîne et son numéro de résidu,
    et contenant ses atomes constitutifs indexés par nom.
    """
    def __init__(self, name, position, chain="", residue_name=None):
        self.name = name                        # base parente : A, C, G ou U
        self.position = position                # numéro de résidu (+ code d'insertion)
        self.chain = chain
        self.residue_name = residue_name or name  # nom PDB d'origine (ex. PSU)
        self.atoms = {}

    def add_atom(self, atom):
        """Ajoute un objet Atom au nucléotide (le premier atome d'un nom donné est conservé)."""
        self.atoms.setdefault(atom.name, atom)

    def count_hydrogen_bonds(self, other_nucleotide):
        """
        Compte les liaisons hydrogènes Watson-Crick (ou Wobble G-U) formées avec
        un autre nucléotide, en ne regardant que les couples donneur/accepteur
        attendus pour ce type de paire. Retourne 0 si la paire n'est pas canonique.
        """
        bonds = CANONICAL_HBONDS.get((self.name, other_nucleotide.name))
        first, second = self, other_nucleotide
        if bonds is None:
            bonds = CANONICAL_HBONDS.get((other_nucleotide.name, self.name))
            first, second = other_nucleotide, self
        if bonds is None:
            return 0

        number_hydrogen_bonds = 0
        for name_a, name_b in bonds:
            atom_a = first.atoms.get(name_a)
            atom_b = second.atoms.get(name_b)
            if atom_a is None or atom_b is None:
                continue
            if HBOND_MIN_DISTANCE <= atom_a.calculate_distance(atom_b) <= HBOND_MAX_DISTANCE:
                number_hydrogen_bonds += 1
        return number_hydrogen_bonds

    def min_hydrogen_bonds(self, other_nucleotide):
        """Nombre minimal de liaisons requises pour apparier ces deux bases, ou None si la paire n'est pas canonique."""
        return MIN_HBONDS.get((self.name, other_nucleotide.name)) or MIN_HBONDS.get((other_nucleotide.name, self.name))
