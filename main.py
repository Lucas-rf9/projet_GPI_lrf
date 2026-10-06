#!/usr/bin/env python3
"""Prédiction de la structure secondaire d'un ARN à partir de ses coordonnées 3D (fichier PDB)."""
import argparse
import sys

from rna_classes import Atom, Nucleotide, nucleotide_base

# Nombre minimal de nucléotides non appariés dans une boucle en épingle :
# deux bases i et j ne peuvent s'apparier que si j - i > MIN_LOOP_LENGTH.
MIN_LOOP_LENGTH = 3

# Symboles utilisés pour les niveaux successifs de pseudo-nœuds.
BRACKETS = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]


def parse_pdb(file_name):
    """
    Lit un fichier PDB et retourne la liste ordonnée des nucléotides d'ARN.
    Le parsing respecte les colonnes à largeur fixe du format PDB, ne lit que
    le premier modèle (structures RMN) et la première conformation alternative.
    """
    nucleotides = {}
    with open(file_name, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue

            alt_loc = line[16]
            if alt_loc not in (" ", "A"):
                continue
            residue_name = line[17:20].strip()
            base = nucleotide_base(residue_name)
            if base is None:  # eau, ions, ligands...
                continue

            atom_name = line[12:16].strip()
            chain = line[21]
            position = line[22:27].strip()  # numéro de résidu + code d'insertion
            atom = Atom(atom_name, float(line[30:38]), float(line[38:46]), float(line[46:54]))

            key = (chain, position)
            if key not in nucleotides:
                nucleotides[key] = Nucleotide(base, position, chain, residue_name)
            nucleotides[key].add_atom(atom)

    # Les dictionnaires conservent l'ordre d'insertion, c'est-à-dire l'ordre du fichier.
    return list(nucleotides.values())


def find_base_pairs(nucleotides):
    """
    Identifie les paires de bases canoniques (G-C, A-U, G-U).
    Chaque nucléotide n'est apparié qu'une seule fois : en cas de conflit,
    la paire présentant le plus de liaisons hydrogènes est retenue.
    Retourne une liste triée de tuples (i, j) avec i < j.
    """
    candidates = []
    for i in range(len(nucleotides)):
        for j in range(i + MIN_LOOP_LENGTH + 1, len(nucleotides)):
            nucleotide_a, nucleotide_b = nucleotides[i], nucleotides[j]
            required = nucleotide_a.min_hydrogen_bonds(nucleotide_b)
            if required is None:
                continue
            number_bonds = nucleotide_a.count_hydrogen_bonds(nucleotide_b)
            if number_bonds >= required:
                candidates.append((number_bonds, i, j))

    # Sélection gloutonne : les paires les plus liées d'abord.
    candidates.sort(key=lambda candidate: (-candidate[0], candidate[1], candidate[2]))
    paired = set()
    pairs = []
    for _, i, j in candidates:
        if i not in paired and j not in paired:
            paired.update((i, j))
            pairs.append((i, j))
    return sorted(pairs)


def pairs_cross(pair_a, pair_b):
    """Indique si deux paires (i, j) se croisent, c'est-à-dire forment un pseudo-nœud."""
    (i, j), (k, l) = pair_a, pair_b
    return i < k < j < l or k < i < l < j


def to_dot_bracket(length, pairs):
    """
    Construit la bracket notation. Les paires qui croisent des paires déjà
    placées (pseudo-nœuds) sont notées avec [], puis {}, puis <>.
    """
    structure = ["."] * length
    levels = [[] for _ in BRACKETS]
    for pair in pairs:
        for level, level_pairs in enumerate(levels):
            if not any(pairs_cross(pair, other) for other in level_pairs):
                level_pairs.append(pair)
                opening, closing = BRACKETS[level]
                structure[pair[0]] = opening
                structure[pair[1]] = closing
                break
    return "".join(structure)


def convert_2d_rna_structure(file_name):
    """
    Lit un fichier PDB d'ARN et retourne sa séquence et sa structure 2D
    au format bracket notation.
    """
    nucleotides = parse_pdb(file_name)
    pairs = find_base_pairs(nucleotides)
    sequence = "".join(nucleotide.name for nucleotide in nucleotides)
    return sequence, to_dot_bracket(len(nucleotides), pairs)


def main():
    parser = argparse.ArgumentParser(description="Prédit la structure secondaire d'un ARN à partir d'un fichier PDB.")
    parser.add_argument("pdb_file", help="fichier PDB contenant la structure 3D de l'ARN (ex. 1EHZ.pdb)")
    args = parser.parse_args()

    try:
        sequence, structure = convert_2d_rna_structure(args.pdb_file)
    except FileNotFoundError:
        sys.exit(f"Erreur : Le fichier '{args.pdb_file}' est introuvable. Vérifiez le chemin d'accès.")
    except ValueError as error:
        sys.exit(f"Erreur : Le fichier '{args.pdb_file}' n'est pas un fichier PDB valide ({error}).")

    if not sequence:
        sys.exit(f"Erreur : Aucun nucléotide d'ARN trouvé dans '{args.pdb_file}'.")

    print(sequence)
    print(structure)


if __name__ == '__main__':
    main()
