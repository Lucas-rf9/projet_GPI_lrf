"""Tests unitaires (lancer avec : python -m unittest)."""
import os
import tempfile
import unittest

from main import convert_2d_rna_structure, find_base_pairs, parse_pdb, to_dot_bracket
from rna_classes import CANONICAL_HBONDS, Atom, Nucleotide

HBOND_LENGTH = 2.9


def pdb_line(serial, atom_name, residue_name, chain, residue_number, x, y, z,
             record="ATOM", alt_loc=" ", insertion_code=" "):
    """Formate une ligne ATOM/HETATM en respectant les colonnes du format PDB."""
    return (f"{record:<6}{serial:>5} {atom_name:<4}{alt_loc}{residue_name:>3} {chain}"
            f"{residue_number:>4}{insertion_code}   {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00\n")


def build_structure(sequence, pairs):
    """
    Construit des nucléotides dont seuls les atomes impliqués dans les paires
    données sont placés à distance de liaison hydrogène ; tous les autres sont éloignés.
    """
    coordinates = {index: [] for index in range(len(sequence))}
    for pair_index, (i, j) in enumerate(pairs):
        bonds = CANONICAL_HBONDS.get((sequence[i], sequence[j]))
        first, second = i, j
        if bonds is None:
            bonds = CANONICAL_HBONDS[(sequence[j], sequence[i])]
            first, second = j, i
        origin_z = 100.0 * (pair_index + 1)
        for bond_index, (name_a, name_b) in enumerate(bonds):
            y = 4.0 * bond_index
            coordinates[first].append((name_a, 0.0, y, origin_z))
            coordinates[second].append((name_b, HBOND_LENGTH, y, origin_z))
    for index in range(len(sequence)):
        # Un atome du squelette, loin de tout le reste.
        coordinates[index].append(("P", -900.0 + 10.0 * index, 0.0, 0.0))
    return coordinates


def write_pdb(sequence, pairs, extra_lines=""):
    coordinates = build_structure(sequence, pairs)
    lines = []
    serial = 1
    for index, base in enumerate(sequence):
        for name, x, y, z in coordinates[index]:
            lines.append(pdb_line(serial, name, base, "A", index + 1, x, y, z))
            serial += 1
    handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    handle.write("".join(lines) + extra_lines + "END\n")
    handle.close()
    return handle.name


class ParsingTest(unittest.TestCase):
    def setUp(self):
        self.files = []

    def tearDown(self):
        for file_name in self.files:
            os.remove(file_name)

    def write(self, content):
        handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
        handle.write(content)
        handle.close()
        self.files.append(handle.name)
        return handle.name

    def test_fixed_columns_with_touching_negative_coordinates(self):
        line = pdb_line(1, "N1", "G", "A", 1, -100.123, -200.456, -300.789)
        self.assertEqual(len(line[30:54].split()), 1)  # un split() naïf échouerait
        nucleotides = parse_pdb(self.write(line))
        atom = nucleotides[0].atoms["N1"]
        self.assertEqual((atom.x, atom.y, atom.z), (-100.123, -200.456, -300.789))

    def test_modified_residues_kept_and_non_nucleotides_skipped(self):
        content = (pdb_line(1, "N1", "G", "A", 1, 0, 0, 0)
                   + pdb_line(2, "N1", "PSU", "A", 2, 5, 0, 0, record="HETATM")
                   + pdb_line(3, "O", "HOH", "A", 101, 9, 0, 0, record="HETATM")
                   + pdb_line(4, "MG", "MG", "A", 102, 9, 9, 0, record="HETATM"))
        nucleotides = parse_pdb(self.write(content))
        self.assertEqual([n.name for n in nucleotides], ["G", "U"])
        self.assertEqual(nucleotides[1].residue_name, "PSU")

    def test_chains_insertion_codes_altlocs_and_models(self):
        content = (pdb_line(1, "N1", "G", "A", 1, 0, 0, 0)
                   + pdb_line(2, "N1", "C", "B", 1, 5, 0, 0)
                   + pdb_line(3, "N1", "A", "B", 1, 9, 0, 0, insertion_code="A")
                   + pdb_line(4, "N3", "A", "B", 1, 1, 1, 1, alt_loc="A", insertion_code="A")
                   + pdb_line(5, "N3", "A", "B", 1, 7, 7, 7, alt_loc="B", insertion_code="A")
                   + "ENDMDL\n"
                   + pdb_line(6, "N1", "U", "C", 1, 0, 0, 0))
        nucleotides = parse_pdb(self.write(content))
        self.assertEqual([(n.chain, n.position, n.name) for n in nucleotides],
                         [("A", "1", "G"), ("B", "1", "C"), ("B", "1A", "A")])
        self.assertEqual(nucleotides[2].atoms["N3"].x, 1.0)


class HydrogenBondTest(unittest.TestCase):
    def test_only_donor_acceptor_atoms_are_counted(self):
        guanine, cytosine = Nucleotide("G", "1"), Nucleotide("C", "2")
        guanine.add_atom(Atom("C8", 0, 0, 0))
        cytosine.add_atom(Atom("C5", 3.0, 0, 0))
        self.assertEqual(guanine.count_hydrogen_bonds(cytosine), 0)
        guanine.add_atom(Atom("N1", 0, 10, 0))
        cytosine.add_atom(Atom("N3", 3.0, 10, 0))
        self.assertEqual(guanine.count_hydrogen_bonds(cytosine), 1)
        self.assertEqual(cytosine.count_hydrogen_bonds(guanine), 1)

    def test_non_canonical_pair_has_no_bonds(self):
        self.assertEqual(Nucleotide("A", "1").count_hydrogen_bonds(Nucleotide("C", "2")), 0)
        self.assertIsNone(Nucleotide("A", "1").min_hydrogen_bonds(Nucleotide("C", "2")))


class StructureTest(unittest.TestCase):
    def test_hairpin(self):
        sequence = "GGAUAAAAUCCU"
        file_name = write_pdb(sequence, [(0, 10), (1, 9), (2, 8)])
        try:
            self.assertEqual(convert_2d_rna_structure(file_name), (sequence, "(((.....)))."))
        finally:
            os.remove(file_name)

    def test_wobble_and_au_pairs(self):
        sequence = "GAAAAAAU" + "AAAAAAU"
        nucleotides = parse_pdb_from(sequence, [(0, 7), (8, 14)])
        self.assertEqual(find_base_pairs(nucleotides), [(0, 7), (8, 14)])

    def test_each_nucleotide_paired_once(self):
        # On déplace C10 au contact de G0, déjà apparié à C5 : G0 a alors deux
        # partenaires possibles mais ne doit apparaître que dans une seule paire.
        sequence = "GAAAACAAAAC"
        nucleotides = parse_pdb_from(sequence, [(0, 5)])
        for name, atom in nucleotides[5].atoms.items():
            if name != "P":
                nucleotides[10].atoms[name] = atom
        pairs = find_base_pairs(nucleotides)
        used = [index for pair in pairs for index in pair]
        self.assertEqual(len(used), len(set(used)))
        self.assertIn((0, 5), pairs)

    def test_minimum_loop_length(self):
        nucleotides = parse_pdb_from("GAAC", [(0, 3)])
        self.assertEqual(find_base_pairs(nucleotides), [])

    def test_pseudoknot_brackets(self):
        self.assertEqual(to_dot_bracket(12, [(0, 6), (1, 5), (3, 10), (4, 9)]),
                         "((.[[))..]].")
        self.assertEqual(to_dot_bracket(8, [(0, 3), (1, 5), (2, 7)]), "([{).].}")


def parse_pdb_from(sequence, pairs):
    file_name = write_pdb(sequence, pairs)
    try:
        return parse_pdb(file_name)
    finally:
        os.remove(file_name)


if __name__ == "__main__":
    unittest.main()
