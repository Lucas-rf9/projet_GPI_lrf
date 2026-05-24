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

* **`rna_classes.py` (Modélisation Orientée Objet) :** Contient les classes `Atom` et `Nucleotide`. Cette structure encapsule les coordonnées 3D (x, y, z) et intègre les méthodes de calcul de distance euclidienne. La détection des liaisons hydrogènes s'effectue via une fenêtre de tolérance spatiale stricte (entre 2.5 Å et 3.5 Å).
* **`main.py` (Exécution et Algorithmique) :** Gère le parsing (lecture) du fichier PDB et l'instanciation des objets. L'algorithme de prédiction applique ensuite un double filtre :
    1. Un seuil quantitatif de liaisons hydrogènes (>= 3 pour G-C, >= 2 pour A-U / G-U).
    2. Un filtre de décalage positionnel (i+4) dans la séquence afin d'ignorer les contacts diagonaux liés à l'empilement (*stacking*) des bases.
* **Génération du rendu :** Les paires identifiées modifient dynamiquement une liste d'états pour construire la séquence et sa *bracket notation* finale.