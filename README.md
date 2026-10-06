# projet_GPI_lrf
# Projet GPI #1 : Prédiction de structure secondaire d'ARN à partir de coordonnées 3D

## 1. Objectif du projet
Ce projet vise à analyser la structure tridimensionnelle d'un ARN à partir d'un fichier PDB afin d'en déduire sa structure secondaire (2D). L'outil extrait les coordonnées spatiales des atomes, calcule les distances interatomiques pour identifier les liaisons hydrogènes, et détermine les paires de bases (Watson-Crick G-C, A-U et Wobble G-U). Le résultat final affiche la séquence nucléotidique superposée à sa structure 2D au format *bracket notation*.

## 2. Installation et Utilisation
Ce projet a été développé en Python 3 et s'appuie uniquement sur la bibliothèque standard. Aucun module externe n'est requis.

Pour récupérer un fichier PDB d'exemple (ici l'ARN de transfert de la levure `1EHZ`), vous pouvez le télécharger directement depuis votre terminal avec la commande suivante :

```bash
curl -O https://files.rcsb.org/download/1EHZ.pdb
```

Pour exécuter le programme, placez-vous dans le répertoire du projet et lancez le script principal en lui passant le fichier PDB récupéré en argument :

```bash
python main.py 1EHZ.pdb
```
*(Note : Le programme inclut une gestion d'erreur via l'interface en ligne de commande si l'argument est manquant).*

## 3. Architecture Logicielle et Modélisation
Le code a été modularisé pour séparer la théorie (logique métier) de l'exécution (scripting) :

* **`rna_classes.py` (Modélisation Orientée Objet) :** Contient les classes `Atom` et `Nucleotide`. Un nucléotide est identifié par sa chaîne et son numéro de résidu, et indexe ses atomes par nom. Les liaisons hydrogènes sont détectées uniquement entre les couples donneur/accepteur attendus pour chaque paire canonique (ex. G:N1–C:N3, G:N2–C:O2, G:O6–C:N4), dans une fenêtre de 2.5 à 3.5 Å. Les nucléotides modifiés courants (PSU, 5MC, 2MG, H2U...) sont ramenés à leur base parente.
* **`main.py` (Exécution et Algorithmique) :**
    1. **Parsing PDB** sur les colonnes à largeur fixe du format (robuste aux coordonnées négatives collées), lecture des lignes `ATOM` et `HETATM`, du premier modèle uniquement et de la première conformation alternative ; l'eau, les ions et les ligands sont ignorés.
    2. **Détection des paires** G-C (3 liaisons), A-U et G-U (2 liaisons), avec une boucle minimale de 3 nucléotides entre deux bases appariées.
    3. **Appariement unique** : chaque nucléotide n'appartient qu'à une seule paire ; en cas de conflit, la paire la plus liée est retenue.
    4. **Bracket notation** avec gestion des pseudo-nœuds : les paires qui se croisent sont notées `[]`, puis `{}`, puis `<>`.
* **`test_main.py` :** tests unitaires sur des structures synthétiques, à lancer avec `python -m unittest`.
