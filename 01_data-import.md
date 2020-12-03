01\_data-import
================

  - [Downloading the Standard Operating Procedure example
    data](#downloading-the-standard-operating-procedure-example-data)
  - [Decompress the data of Mothur MiSeq
    SOP](#decompress-the-data-of-mothur-miseq-sop)
  - [charger la base de donnée de
    silva](#charger-la-base-de-donnée-de-silva)

# Downloading the Standard Operating Procedure example data

Télécharger les données utilisées dans le Mothur MiSeq SOP. Ces fichiers
fastq ont été générés par séquençage d’amplicon 2x250 Illumina Miseq de
la région V4 du gène de l’ARNr 16S à partir d’échantillons d’intestin
prélevés longitudinalement à partir d’une souris après le sevrage.

``` bash
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
```

# Decompress the data of Mothur MiSeq SOP

``` bash
unzip miseqsopdata.zip
```

# charger la base de donnée de silva

SILVA fournit des ensembles de données complets,contrôlés par la qualité
et régulièrement mis à jour de séquences d’ARN ribosomal (ARNr) alignées
de petite (16S / 18S, SSU) et de grande sous-unité (23S / 28S, LSU) pour
les trois domaines de la vie (bactéries, archées et eucarya).

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

Lier l’annotation d’espèce à l’annotation taxonomique

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```
