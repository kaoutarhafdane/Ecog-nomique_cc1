metabarcoding with dada2: environment installation
================

  - [update VM configuration](#update-vm-configuration)
  - [package install](#package-install)
  - [installation de package
    phangorn](#installation-de-package-phangorn)
  - [installation de DECIPHER](#installation-de-decipher)
  - [installation de gridExtra](#installation-de-gridextra)
  - [File description](#file-description)
  - [installation de package DESeq2](#installation-de-package-deseq2)
  - [installation de structSSI](#installation-de-structssi)
  - [installation de package
    phyloseq](#installation-de-package-phyloseq)
  - [installation de package ggplot2](#installation-de-package-ggplot2)
  - [installation de package
    Biostrings](#installation-de-package-biostrings)

# update VM configuration

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
```

# package install

Following instruction on
<https://benjjneb.github.io/dada2/dada-installation.html> cela permet
d’installer une première version de dada2 et puis obtenir la version
actuelle de dada2 en mettant à jour vers Bioconductor 3.11.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.11')
BiocManager::install("dada2", version = "3.11")
```

    ## Warning in install.packages(...): installation of package 'Biostrings' had non-
    ## zero exit status

    ## Warning in install.packages(...): installation of package 'dada2' had non-zero
    ## exit status

# installation de package phangorn

ce package contient des méthodes d’estimation des arbres et des réseaux
phylogénétiques

``` r
BiocManager::install("phangorn")
```

# installation de DECIPHER

Outils de conservation, d’analyse et de manipulation de séquences
biologiques

``` r
BiocManager::install("DECIPHER")
```

    ## Warning in install.packages(...): installation of package 'Biostrings' had non-
    ## zero exit status

    ## Warning in install.packages(...): installation of package 'DECIPHER' had non-
    ## zero exit status

# installation de gridExtra

Fournit un certain nombre de fonctions au niveau de l’utilisateur pour
travailler avec des graphiques «en grille», notamment pour organiser
plusieurs tracés basés sur une grille sur une page et dessiner des
tableaux.

``` r
install.packages("gridExtra")
```

# File description

Ici on descrit cran\_package en tant que Shiny…. on fait la même chose
pour github\_package on disant que c’est jfu…

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
```

on install cran\_package et github\_package et aussi bioc\_package

``` r
install.packages(.cran_packages)
devtools::install_github(.github_packages)
BiocManager::install(.bioc_packages)
library(vegan)
```

# installation de package DESeq2

Estimer la dépendance de la variance-moyenne dans les données de
dénombrement à partir de tests de séquençage à haut débit et tester
l’expression différentielle basée sur un modèle utilisant la
distribution binomiale négative

``` r
BiocManager::install("DESeq2")
```

# installation de structSSI

la fonction install\_local est vectorisée, elle permet d’installer
plusieurs packages en une seule commande.

``` r
library(devtools)
install_local("./structSSI_1.1.1.tar.gz")
```

``` bash
wget https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
```

# installation de package phyloseq

Le package phyloseq est un outil pour importer, stocker, analyser et
afficher graphiquement des données de séquençage phylogénétique
complexes qui ont déjà été regroupées en unités taxonomiques
opérationnelles (OTU), en particulier lorsqu’il existe des échantillons
de données associés, un arbre phylogénétique, et / ou attribution
taxinomique des OTU.

``` r
BiocManager::install('phyloseq')
```

    ## Warning in install.packages(...): installation of package 'Biostrings' had non-
    ## zero exit status

    ## Warning in install.packages(...): installation of package 'phyloseq' had non-
    ## zero exit status

# installation de package ggplot2

le package ggplot2 permet de créer une visualisation graphique élégante
des données.

``` r
BiocManager::install('ggplot2')
```

# installation de package Biostrings

le package biostrings permet une manipulation efficace des chaînes
biologiques.

``` r
BiocManager::install('Biostrings')
```

    ## Warning in install.packages(...): installation of package 'Biostrings' had non-
    ## zero exit status
