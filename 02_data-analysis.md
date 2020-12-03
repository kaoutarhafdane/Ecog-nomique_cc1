Dada2 tutorial
================

  - [getting ready](#getting-ready)
  - [Inspect read quality profiles](#inspect-read-quality-profiles)
  - [Filter and trim](#filter-and-trim)
  - [appretissage des erreurs](#appretissage-des-erreurs)
  - [sample inference](#sample-inference)
  - [mairged paired reads](#mairged-paired-reads)
  - [construire la table
    d’observation](#construire-la-table-dobservation)
  - [Remove chimeras](#remove-chimeras)
  - [Track reads through the
    pipeline](#track-reads-through-the-pipeline)
  - [assignation de la taxonomie](#assignation-de-la-taxonomie)
  - [Evaluate accuracy](#evaluate-accuracy)

# getting ready

Nous chargeons d’abord le dada2package

``` r
library("dada2")
```

Définir la variable path suivante pour qu’elle pointe vers le répertoire
Miseq\_sop extrait sur la machine. Les données utilisées sont dans le
Mothur MiSeq SOP. Ces fichiers fastq ont été générés par séquençage
d’amplicon 2x250 Illumina Miseq de la région V4 du gène de l’ARNr 16S
à partir d’échantillons d’intestin prélevés longitudinalement à partir
d’une souris après le sevrage. Pour l’instant, considérez simplement les
fichiers fastq appariés à traiter.

``` r
path <- "~/github/ecog2_cc1_final/ecog2_cc1/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

cela permet de lire les noms des fichiers fastq et effectuons quelques
manipulations de chaînes pour obtenir des listes correspondantes des
fichiers fastq forward et reverse. On va créer une variable fnFs et on
lui assigne la valeur de résultat de la fonction sort qui va classer les
résultats de la fonction list.filesqui va lister le fichier
R1\_001.fastq et va faire la même chose pour R2\_001.fastq. En suite on
va extrairer les sample names avec la fonction strsplit(), en supposant
que les noms de fichiers ont le format: SAMPLENAME\_XXX.fastq

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

# Inspect read quality profiles

Nous commençons par visualiser les profils de qualité des Forward reads
en utilisant la fonction plotQualityProfile qui va permettre de tracer
un résumé visuel de la distribution des scores de qualité en fonction de
la position de la séquence pour le fichier: fnFs fastq d’entrée.

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_data-analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Les lectures avant sont de bonne qualité. mais on rogner quand même les
derniers nucléotides pour éviter des erreurs moins bien contrôlées qui
peuvent s’y produire. on va tronquer donc les lectures avant à la
position 240 (en coupant les 10 derniers nucléotides).

on visualise le profil de qualité des reverse reads en utilisant la
fonction plotQualityProfile qui va permettre de tracer un résumé visuel
de la distribution des scores de qualité en fonction de la position de
la séquence pour le fichier: fnRs fastq d’entrée.

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_data-analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Comentaire de graphe au dessus: sur ce graphe on remarque que la qualité
est moins bonne que celle du graphe d’avant, on tranque donc les
lectures inversées à la position 160 où la distribution de qualité se
bloque.

# Filter and trim

Attribuer les noms de fichiers aux fichiers fastq.gz filtrés. on crée
une variable filtFs et on met dedans \_F\_filt.fastq.gz et puis une
filtRs et on met dedans \_R\_filt.fastq.gz après on nome l’objet filtFs
en utilisant la fonction names on fait pareil pour filtRs et on leur
donne le nom sample.names qui est la valeur de la fonction.

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Nous utiliserons les paramètres de filtrage standard: maxN=0(DADA2 ne
nécessite aucun Ns) truncQ=2, rm.phix=TRUEet maxEE=2. Le aramètre
maxEEp définit le nombre maximum d ’«erreurs attendues» autorisées dans
une lecture. Ici, on crée une variale out et on lui assigne les valeurs
des résultats de la fonction filterAndTrim(), qui va filtrer et ajuster
les fichiers fnFs, filtFs, fnRs, filtRs de fastq d’entrée (peut être
compressé) en fonction de plusieurs critères: maxN=0, maxEE=c(2,2),
truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE , et génère des
fichiers fastq (compressés par défaut) contenant les lectures coupées
qui ont passé les filtres. Des fichiers fastq revers et forward
correspondants peuvent être fournis en entrée, auquel cas le filtrage
est effectué sur les lectures avant et arrière indépendamment, et les
deux lectures doivent passer pour que la paire de lecture soit sortie.
En suite on va utiliser la fonction head pour avoir un apperçu de
l’objet out

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

# appretissage des erreurs

dada2 calcul un model d’erreurs apartir des données de séquençage, cette
méthode sur les reads F Reverse. L’algorithme DADA2 utilise un modèle
d’erreur paramétrique ( err) et chaque jeu de données d’amplicon a un
ensemble différent de taux d’erreur. La learnErrors méthode apprend ce
modèle d’erreur à partir des données, en alternant l’estimation des taux
d’erreur et l’inférence de la composition de l’échantillon jusqu’à ce
qu’ils convergent vers une solution cohérente conjointement.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

dada2 calcul un model d’erreurs apartir des données de séquençage, cette
méthode sur les reads Reverse de la même manière que pour errF

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

le code suivant est pour vérifier si rien d’autre, de visualiser les
taux d’erreur estimés en utilisant la fonction plotErrors. Cette
fonction va tracer la fréquence observée de chaque transition (par
exemple A-\> C) en fonction du score de qualité associé. Il trace
également les taux d’erreur estimés finaux (s’ils existent). l’argument
nominalQ=TRUE va permettre de tracer les taux d’erreur attendus (ligne
rouge) si les scores de qualité correspondent exactement à leur
définition nominale: Q = -10 log10 (p\_err).

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_data-analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
Intérprétation de graphe précédent: graphe 1: la probabilité que A est A
doit être max ainsi de suite… les taux d’erreur pour chaque transition
sont indiquées sur les graphes (graphe 2: A–\>C etc) c2: probabilité
d’erreus de séq pour qu’un A soit un C. Ici, les taux d’erreur estimés
(ligne noire) correspondent bien aux taux observés (points), et les taux
d’erreur diminuent avec une qualité accrue comme prévu. Tout semble
raisonnable et nous procédons en toute confiance.

# sample inference

on crée une nouvelle variable dadaFs pour corriger les jeux de données
dada appliquée au donné Forward. La fonction dada supprime toutes les
erreurs de séquençage pour révéler les membres de la communauté
séquencée. l’argument multithread=TRUE: le multithreading est activé
et le nombre de threads disponibles est automatiquement déterminé. Si un
entier est fourni, le nombre de threads à utiliser est défini en
transmettant l’argument à setThreadOptions.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

dada appliqué au donné reverse on crée une nouvelle variable dadaRs et
on lui assigne la valeur de résultat de la fonction dada comme on a fait
pour les reverse à fin de corriger les jeux de données appliquées pour
les Reverse.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

dada a crée des objet de class dada: dadaFs et dada Rs on regarde ce qui
est dans le premier étagers du placard de dadaFS, on peut changer le 1
pour regarder à n’importe quel étage de dadaFS.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# mairged paired reads

on va merger Read 1 et read 2 pour obtenir les séquences entièrement
débruitées en utilisant la fonction mergePairs(). La fusion est
effectuée en alignant les lectures avant débruitées avec le complément
inverse des lectures inverses débruitées correspondantes, puis en
construisant les séquences «contig» fusionnées. verbose=TRUE:monter les
étape avec le texte head(mergers\[\[1\]\]):inspecter les résultats en
regardant les premières lignes

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6551 paired-reads (in 106 unique pairings) successfully merged out of 6907 (in 199 pairings) input.

    ## 5025 paired-reads (in 100 unique pairings) successfully merged out of 5188 (in 156 pairings) input.

    ## 4973 paired-reads (in 80 unique pairings) successfully merged out of 5268 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2756 (in 109 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3622 paired-reads (in 53 unique pairings) successfully merged out of 4103 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6515 (in 198 pairings) input.

    ## 3961 paired-reads (in 90 unique pairings) successfully merged out of 4384 (in 188 pairings) input.

    ## 14231 paired-reads (in 143 unique pairings) successfully merged out of 15358 (in 351 pairings) input.

    ## 10526 paired-reads (in 120 unique pairings) successfully merged out of 11166 (in 279 pairings) input.

    ## 11156 paired-reads (in 137 unique pairings) successfully merged out of 11799 (in 298 pairings) input.

    ## 4329 paired-reads (in 84 unique pairings) successfully merged out of 4788 (in 180 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7193 (in 187 pairings) input.

    ## 4430 paired-reads (in 67 unique pairings) successfully merged out of 4605 (in 127 pairings) input.

    ## 4574 paired-reads (in 100 unique pairings) successfully merged out of 4736 (in 172 pairings) input.

    ## 6094 paired-reads (in 109 unique pairings) successfully merged out of 6314 (in 172 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

Donc là on a les reads fusionnés avec succès.

# construire la table d’observation

à partir des merged, on crée un nouvelle objet seqtab la fonction va
permettre de construire une table de séquence (analogue à une table OTU)
à partir de la liste d’échantillon mergers. la fonction dim va permettre
de récupérer l’objet seqtab.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

regarder la distribution de la longeur des séquences. La table de
séquence est une matrix avec des lignes correspondant aux (et nommées
par) les échantillons, et des colonnes correspondant (et nommées par)
les variantes de séquence.

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

dans mon objet seqtab j’ai une séq qui fait 251 paires de bases, une
autre qui fait 252 etc

# Remove chimeras

une chimère: ça se passe pendant l’amplification par PCR donc par ex un
ADN 16s amplifier par un fragment reverse et forrward. Si on prend que
le forward,y aura élongation mais imaginant qu’elle s’arréte avant la
fin de la séq 16S. Donc on va avoir le fragment 16S qui n’a pas bouger
et un fragment non complet. donc àprès le 2ème cycle, le fragment non
complet va pouvoir s’hybrider avec un 16s d’une autre bactérie, et
l’élongation va continuer. on va avoir comme résultat au final un
fragment hybride qui provient du premier ARN 16 et du deuxième. cela
s’appelle chimère. On va éliminer ces chimères en utlisant la fonction
removeBimeraDenovo le système va regarder tout les séq rares dans le
début contig correspont au premier ARN et la fin au deuxième.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

    ## [1]  20 232

Donc le résultat montre qui il a détectrer 61 (=293-232) chimère sur
293(= 1+88+196+6+2).

calcul de ratio chimère qui est égale à la somme des sequences se
trouvant dans l’objet seqtab.mochim (c’est l’objet après remove des
chimères) / somme des séquences se traouvent dans l’objet seqtab(avant
d’enlever les chimères)

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.964263

y a (1-0.96)\*100 = 3.5% de séq chimériques dans notre jeu de données.
les chimères représentent environ 21% des variantes de séquence
fusionnées, mais lorsque nous tenons compte de l’abondance de ces
variantes, nous voyons qu’elles ne représentent qu’environ 4% des
lectures de séquence fusionnée.

# Track reads through the pipeline

résumé des fichiers qualité. construire une table on crée un nouvel
objet qui est getN c’est une variable qui prend le rôle d’une fonction
apliquer la fonnction get N qui est la somme des get uniq de dadaFS la
table track va être la concaténation de tout ce qui est entre parentèse.
head(track) : pour visualiser le tableau

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6996      6978   6551    6539
    ## F3D1    5869     5299      5227      5239   5025    5014
    ## F3D141  5958     5463      5339      5351   4973    4850
    ## F3D142  3183     2914      2799      2833   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4146      4224   3622    3483

On a conservé la majorité des lectures brutes et aucune baisse trop
importante n’est associée à une seule étape.

# assignation de la taxonomie

il va regarder dans le base de données et à partir des séq qui sont
proches, et va assigner une taxonomie.

c’est une façon d’attribuer une taxonomie aux séquences. La fonction
assignTaxonomy prend en entrée un ensemble de séquences à classer et un
ensemble d’apprentissage de séquences de référence avec une taxonomie
connue (silva en ce cas), et produit des affectations taxonomiques avec
au moins une minBootconfiance bootstrap.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

ajout d’espèces, téléchargement du fichier et le placer dans le
répertoire taxa contenant les fichiers fastq.

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

Inspecter les affectations taxonomiques:

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

les Bacteroidetes sont bien représentés parmi les taxons les plus
abondants dans ces échantillons fécaux. Peu d’attributions d’espèces ont
été faites, à la fois parce qu’il est souvent impossible de faire des
assignations d’espèces sans ambiguïté à partir de sous-segments du gène
16S, et parce qu’il y a étonnamment peu de couverture du microbiote
intestinal de souris indigène dans les bases de données de référence.

# Evaluate accuracy

L’un des échantillons inclus ici était une «communauté fictive», dans
laquelle un mélange de 20 souches connues a été séquencé (cette
communauté fictive est supposée être de 21 souches, mais P. acnes est
absent des données brutes). Les séquences de référence correspondant à
ces souches ont été fournies dans l’archive zip téléchargée. Nous
revenons à cet échantillon et comparons les variantes de séquence
inférées par DADA2 à la composition attendue de la communauté.

cat permet de faire sortir les objets en concaténant les représentations

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

On décrit match.ref

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

Cette fausse communauté contenait 20 souches bactériennes. DADA2 a
identifié 20 ASV qui correspondent tous exactement aux génomes de
référence des membres attendus de la communauté. Le taux d’erreur
résiduel après le pipeline DADA2 pour cet échantillon est donc de 0% .

``` r
save.image(file="02_data-analysis")
```
