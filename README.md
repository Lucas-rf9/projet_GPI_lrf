# projet_GPI_lrf
# RNA PDB to Dot-Bracket Converter

> De la 3D à la 2D : extraction de la structure secondaire d'un ARN par calcul de distances atomiques.

## 📝 Description
Ce petit projet éducatif en Python permet de lire un fichier PDB d'ARN et d'en déduire sa structure secondaire au format **dot-bracket** (ex: `((((....))))`). 

## ✨ Fonctionnalités
- **Extraction de séquence** : Lit le fichier PDB et extrait la séquence primaire d'ARN.
- **Détection des paires de bases** : Identifie les paires canoniques (Watson-Crick A-U, G-C) et Wobble (G-U) via un seuil de distance 3D.
- **Nettoyage des conflits** : Résolution basique des bases appariées multiples en conservant la distance la plus courte.
- **Génération Dot-Bracket** : Exporte la structure sous la forme standard de points et de parenthèses.

## 🛠️ Prérequis
