# MetIDfyR V2 – BioTransformer

Nouvelle version de **MetIDfyR**, un package R pour prédire et détecter les métabolites des xénobiotiques en spectrométrie de masse haute résolution (HRMS).

Dans cette version, les métabolites sont prédits avec **BioTransformer** à partir du SMILES de la molécule parent.

La version 1 est décrite dans : Delcourt V. et al., *MetIDfyR, an open-source R package to decipher small-molecule drugs metabolism through high-resolution mass spectrometry*, Analytical Chemistry, 2020.

## Organisation du dossier

```
MetIDfyR_V2/
├── README.md
├── .gitignore
└── BioTransformer/
    ├── metID_V2.R                  script principal
    ├── diff.py                     fonctions Python (RDKit)
    ├── util_smiles-copie.R         fonctions (BioTransformer, chromatogrammes, MS2)
    ├── plot_functions.R            fonctions graphiques
    ├── biotransformer-1.1.5.jar    à télécharger (non inclus, voir Installation)
    ├── supportfiles/               à télécharger avec BioTransformer (non inclus)
    └── input/
        ├── config_LGD_smiles.R     exemple de fichier de configuration
        └── lgd_DIA_peak-picked_smiles.tsv   exemple de fichier d'entrée
```

## Installation

1. **R** avec les packages : pacman, optparse, BiocManager, MSnbase, Rdisop, reticulate, dplyr, plyr, doParallel, stringr, readr, tibble, arrangements, snow, svglite, ggplot2, ggpubr, ggrepel, htmlwidgets. Le script installe automatiquement ceux qui manquent.
2. **Java** pour lancer BioTransformer.
3. **BioTransformer 1.1.5** (non inclus dans ce dépôt) : téléchargez BioTransformer 1.1.5 depuis le site officiel (Wishart Lab). Placez le fichier `biotransformer-1.1.5.jar` **et** le dossier `supportfiles/` dans le dossier `BioTransformer/`.
4. **Environnement Python avec RDKit**, nommé `rdkit-env` :

```r
library(reticulate)
virtualenv_create("rdkit-env")
py_install("rdkit", envname = "rdkit-env")
```

## Utilisation

```bash
Rscript BioTransformer/metID_V2.R \
  -i BioTransformer/input/lgd_DIA_peak-picked_smiles.tsv \
  -o output \
  -c BioTransformer/input/config_LGD_smiles.R
```

Options :

- `-i` : fichier d'entrée (.tsv)
- `-o` : dossier de sortie
- `-c` : fichier de configuration (.R)

Le script peut être lancé depuis n'importe quel dossier.

## Fichier d'entrée (.tsv)

Une ligne par molécule, avec les colonnes :

| Colonne | Description |
|---|---|
| `name` | nom de la molécule |
| `smiles` | structure SMILES (obligatoire) |
| `ms_file` | chemin du fichier de données MS |
| `adduct_plus` | adduit en mode positif (ex. H), ou vide |
| `adduct_minus` | adduit en mode négatif (ex. H), ou vide |
| `ms2_reference_tsv` | fichier MS2 de référence (optionnel) |

## Fichier de configuration (.R)

Il définit les paramètres :

- `bool_phase_1`, `bool_phase_2` : phases de métabolisme à prédire. Phase I seule : mode `cyp450`. Phase II seule : mode `phaseII`. Les deux : mode `allHuman`.
- `do_smiles`
- `cores` : nombre de cœurs pour le calcul parallèle
- `lib_perso` : dossier de librairies R personnelles
- `min_peak_intensity`, `mz_ppm`, `nb_scan`, `minimum_mz`, `rt_windows`, `wdw_mz_ms2` : paramètres de détection des pics

## Résultats

Pour chaque molécule, un dossier `out_<molécule>_<fichier>_<date>/` contient :

- `POS/`, `NEG/`, `POSNEG/` : figures (.svg) et tables MS2 (.tsv) pour chaque métabolite détecté
- `images/` : structures des métabolites (.svg)
- `input/` : copie de la configuration et de la molécule analysée
- `out_<fichier>.tsv` : tableau final des métabolites, classés par score

## Auteur

Awa Bousso GUEYE – Laboratoire des Courses Hippiques (LCH), Verrières-le-Buisson.
