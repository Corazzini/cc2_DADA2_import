analyse Dada2
================

  - [Preparation](#preparation)
  - [Score de qualité des reads](#score-de-qualité-des-reads)
  - [Filtration des données](#filtration-des-données)
  - [Modèle d’erreur](#modèle-derreur)
  - [Inférence d’échantillon](#inférence-déchantillon)
  - [Fusionner les reads appariées](#fusionner-les-reads-appariées)
  - [Construire une table de
    séquence](#construire-une-table-de-séquence)
  - [Supprimer les chimères](#supprimer-les-chimères)
  - [Suivre les reads dans la
    pipeline](#suivre-les-reads-dans-la-pipeline)
  - [Assigniation taxonomique](#assigniation-taxonomique)

# Preparation

``` r
# Appel la library
library(dada2)
```

    ## Loading required package: Rcpp

``` r
#Importer les jeux de données dans paths
path <- "~/cc2_DADA2_import/CC2"
list.files(path)
```

    ##  [1] "filtered"                            "Station5_Fond1_10sept14_R1.fastq"   
    ##  [3] "Station5_Fond1_10sept14_R2.fastq"    "Station5_Fond1_11mars15_R1.fastq"   
    ##  [5] "Station5_Fond1_11mars15_R2.fastq"    "Station5_Fond2_10sept14_R1.fastq"   
    ##  [7] "Station5_Fond2_10sept14_R2.fastq"    "Station5_Fond2_11mars15_R1.fastq"   
    ##  [9] "Station5_Fond2_11mars15_R2.fastq"    "Station5_Fond3_10sept14_R1.fastq"   
    ## [11] "Station5_Fond3_10sept14_R2.fastq"    "Station5_Median1_10sept14_R1.fastq" 
    ## [13] "Station5_Median1_10sept14_R2.fastq"  "Station5_Median2_10sept14_R1.fastq" 
    ## [15] "Station5_Median2_10sept14_R2.fastq"  "Station5_Surface1_10sept14_R1.fastq"
    ## [17] "Station5_Surface1_10sept14_R2.fastq" "Station5_Surface1_11mars15_R1.fastq"
    ## [19] "Station5_Surface1_11mars15_R2.fastq" "Station5_Surface2_10sept14_R1.fastq"
    ## [21] "Station5_Surface2_10sept14_R2.fastq" "Station5_Surface2_11mars15_R1.fastq"
    ## [23] "Station5_Surface2_11mars15_R2.fastq"

``` r
fnFs <- sort(list.files(path, pattern="1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="2.fastq", full.names = TRUE))
sample.namesfnFs <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
sample.namesfnRs <- sapply(strsplit(basename(fnRs), "\\."), `[`, 1)
sample.namesfnFs
```

    ##  [1] "Station5_Fond1_10sept14_R1"    "Station5_Fond1_11mars15_R1"   
    ##  [3] "Station5_Fond2_10sept14_R1"    "Station5_Fond2_11mars15_R1"   
    ##  [5] "Station5_Fond3_10sept14_R1"    "Station5_Median1_10sept14_R1" 
    ##  [7] "Station5_Median2_10sept14_R1"  "Station5_Surface1_10sept14_R1"
    ##  [9] "Station5_Surface1_11mars15_R1" "Station5_Surface2_10sept14_R1"
    ## [11] "Station5_Surface2_11mars15_R1"

``` r
sample.namesfnRs
```

    ##  [1] "Station5_Fond1_10sept14_R2"    "Station5_Fond1_11mars15_R2"   
    ##  [3] "Station5_Fond2_10sept14_R2"    "Station5_Fond2_11mars15_R2"   
    ##  [5] "Station5_Fond3_10sept14_R2"    "Station5_Median1_10sept14_R2" 
    ##  [7] "Station5_Median2_10sept14_R2"  "Station5_Surface1_10sept14_R2"
    ##  [9] "Station5_Surface1_11mars15_R2" "Station5_Surface2_10sept14_R2"
    ## [11] "Station5_Surface2_11mars15_R2"


On sépare maintenant les R1 et les R2, pour cela on créer une variable
FnFs qui contiendra les R1 et fnRs qui contiendra les R2.

# Score de qualité des reads

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_analyse_DADA2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_analyse_DADA2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> 


La fonction plotQualityProfile permet de créer une graphique permettant de
visualiser les scores de qualités.  
En abscisse nous avons la position des paires de bases allant de 0 à 250
pb. En ordonnée nous avons le score de qualité. La ligne en vert
correspond au score de qualité moyen pour chaque position. La ligne en
rouge correspond au seuil où le score de qualité est de 10. La ligne en
orange correspond a quartile de la distribution du score de qualité.
Concernant les reads Forward, on peut voir que en général le score de
qualité est plutot bon, on descend pas en dessous du Q30 avant la
position 240pb. Concernant les reads Revers, on descend en dessous du
Q30 a partir de la position 180pb.

# Filtration des données

``` r
# Placer les fichiers filtrés dans filtré
filtFs <- file.path(path, "filtered", paste0(sample.namesfnFs, "_R1.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.namesfnRs, "_R2.fastq"))
names(filtFs) <- sample.namesfnFs
names(filtRs) <- sample.namesfnRs
filtFs
```

    ##                                               Station5_Fond1_10sept14_R1 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond1_10sept14_R1_R1.fastq" 
    ##                                               Station5_Fond1_11mars15_R1 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond1_11mars15_R1_R1.fastq" 
    ##                                               Station5_Fond2_10sept14_R1 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond2_10sept14_R1_R1.fastq" 
    ##                                               Station5_Fond2_11mars15_R1 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond2_11mars15_R1_R1.fastq" 
    ##                                               Station5_Fond3_10sept14_R1 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond3_10sept14_R1_R1.fastq" 
    ##                                             Station5_Median1_10sept14_R1 
    ##  "~/cc2_DADA2_import/CC2/filtered/Station5_Median1_10sept14_R1_R1.fastq" 
    ##                                             Station5_Median2_10sept14_R1 
    ##  "~/cc2_DADA2_import/CC2/filtered/Station5_Median2_10sept14_R1_R1.fastq" 
    ##                                            Station5_Surface1_10sept14_R1 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface1_10sept14_R1_R1.fastq" 
    ##                                            Station5_Surface1_11mars15_R1 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface1_11mars15_R1_R1.fastq" 
    ##                                            Station5_Surface2_10sept14_R1 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface2_10sept14_R1_R1.fastq" 
    ##                                            Station5_Surface2_11mars15_R1 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface2_11mars15_R1_R1.fastq"

``` r
filtRs
```

    ##                                               Station5_Fond1_10sept14_R2 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond1_10sept14_R2_R2.fastq" 
    ##                                               Station5_Fond1_11mars15_R2 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond1_11mars15_R2_R2.fastq" 
    ##                                               Station5_Fond2_10sept14_R2 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond2_10sept14_R2_R2.fastq" 
    ##                                               Station5_Fond2_11mars15_R2 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond2_11mars15_R2_R2.fastq" 
    ##                                               Station5_Fond3_10sept14_R2 
    ##    "~/cc2_DADA2_import/CC2/filtered/Station5_Fond3_10sept14_R2_R2.fastq" 
    ##                                             Station5_Median1_10sept14_R2 
    ##  "~/cc2_DADA2_import/CC2/filtered/Station5_Median1_10sept14_R2_R2.fastq" 
    ##                                             Station5_Median2_10sept14_R2 
    ##  "~/cc2_DADA2_import/CC2/filtered/Station5_Median2_10sept14_R2_R2.fastq" 
    ##                                            Station5_Surface1_10sept14_R2 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface1_10sept14_R2_R2.fastq" 
    ##                                            Station5_Surface1_11mars15_R2 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface1_11mars15_R2_R2.fastq" 
    ##                                            Station5_Surface2_10sept14_R2 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface2_10sept14_R2_R2.fastq" 
    ##                                            Station5_Surface2_11mars15_R2 
    ## "~/cc2_DADA2_import/CC2/filtered/Station5_Surface2_11mars15_R2_R2.fastq"


Ici on va ranger les fichiers dans un dossier nommé filtered contenant
les objets filtFs et filtRs.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 21, truncLen=c(240,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                                    reads.in reads.out
    ## Station5_Fond1_10sept14_R1.fastq     159971    145448
    ## Station5_Fond1_11mars15_R1.fastq     175993    160423
    ## Station5_Fond2_10sept14_R1.fastq     197039    177018
    ## Station5_Fond2_11mars15_R1.fastq      87585     79989
    ## Station5_Fond3_10sept14_R1.fastq     117140    106150
    ## Station5_Median1_10sept14_R1.fastq   116519    106745


La fonction filterAndTrim permet de filtrer et couper les reads Forward
et les reverse.  
La fonction trimLeft = 21 permet d’éliminer les primers pour les reads
Forward et les reads Reverse. la fonction truncLen permet d’éliminer les
nucléotides en position 240pb et 200pb pour conserver le meilleure score
qualité pour les reads(au dessus du Q30). maxEE permet de recalculer le
Qscore moyen apres avoir coupé une partie du read incriminé. MaxN=0
permet d’enlever toutes les bases dans lesquelles il y aura un N (A,T,G
ou C) dans un read d’un jeu de données (le R1 et le R2). On peut voir
qu’on a pas perdu beaucoup de read après la filtration.

# Modèle d’erreur

DADA2 calcul un model d’erreur à partir des données de séquençage. On
applique cette méthode sur les reads forward puis reverse

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 105752691 total bases in 482889 reads from 3 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100755162 total bases in 562878 reads from 4 samples will be used for learning the error rates.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_analyse_DADA2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->


DADA2 analyse les variations de séquences et il va identifier et créer
un modèle d’erreur grâce a la fonction learnErrors. Ce modèle d’erreur
sera ensuite utiliser afin de corriger les reads du jeu de données.

Ce qu’on observe ici est un plot du modèle d’erreur généré par DADA2. En
abscisse nous avont le Qscore et en ordonner la probabilité. On obtient
donc la probabilité d’une mutation en fonction du Qscore. Pour A2A, la
pobabilité qu’un A devient un A est très forte. Pour A2C, lorsque le
Qscore est très élevé, la probabilité qu’un A devient un C est faible.
Si le Qscore est faible, la probablité qu’un A donne un C est élevé. A
l’inverse si le Qscore est élevé alors la probabilité qu’un A donne un
C est faible. La courbe en noir correspond au modèle d’erreur généré par
DADA2.

# Inférence d’échantillon

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 37907 unique sequences.
    ## Sample 2 - 160423 reads in 35863 unique sequences.
    ## Sample 3 - 177018 reads in 47212 unique sequences.
    ## Sample 4 - 79989 reads in 20356 unique sequences.
    ## Sample 5 - 106150 reads in 30255 unique sequences.
    ## Sample 6 - 106745 reads in 28836 unique sequences.
    ## Sample 7 - 98823 reads in 25824 unique sequences.
    ## Sample 8 - 107427 reads in 26733 unique sequences.
    ## Sample 9 - 71082 reads in 17976 unique sequences.
    ## Sample 10 - 78645 reads in 20422 unique sequences.
    ## Sample 11 - 91534 reads in 24487 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 45486 unique sequences.
    ## Sample 2 - 160423 reads in 41638 unique sequences.
    ## Sample 3 - 177018 reads in 55554 unique sequences.
    ## Sample 4 - 79989 reads in 23239 unique sequences.
    ## Sample 5 - 106150 reads in 34625 unique sequences.
    ## Sample 6 - 106745 reads in 31673 unique sequences.
    ## Sample 7 - 98823 reads in 29093 unique sequences.
    ## Sample 8 - 107427 reads in 28947 unique sequences.
    ## Sample 9 - 71082 reads in 21426 unique sequences.
    ## Sample 10 - 78645 reads in 22051 unique sequences.
    ## Sample 11 - 91534 reads in 28266 unique sequences.


L’objet dadaFs reçoit le modèle d’erreur pour les reads forward et
l’objet dadaRs reçoit le modèle d’erreur pour les reads revers Pour le
1er échantillon, on avait 145448 reads et 45486 séquence unique avant la
correction par DADA2.

# Fusionner les reads appariées

Aligner les R1 et les R2 en un contigs

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 117318 paired-reads (in 5196 unique pairings) successfully merged out of 141000 (in 21451 pairings) input.

    ## 138940 paired-reads (in 4296 unique pairings) successfully merged out of 156462 (in 15709 pairings) input.

    ## 142188 paired-reads (in 6989 unique pairings) successfully merged out of 171439 (in 27056 pairings) input.

    ## 67622 paired-reads (in 2721 unique pairings) successfully merged out of 77764 (in 9556 pairings) input.

    ## 83613 paired-reads (in 3458 unique pairings) successfully merged out of 102224 (in 16304 pairings) input.

    ## 86212 paired-reads (in 3348 unique pairings) successfully merged out of 103447 (in 14293 pairings) input.

    ## 80661 paired-reads (in 2727 unique pairings) successfully merged out of 95866 (in 12350 pairings) input.

    ## 89385 paired-reads (in 3073 unique pairings) successfully merged out of 104354 (in 12135 pairings) input.

    ## 59716 paired-reads (in 1939 unique pairings) successfully merged out of 68711 (in 7974 pairings) input.

    ## 66157 paired-reads (in 1763 unique pairings) successfully merged out of 76701 (in 8283 pairings) input.

    ## 75048 paired-reads (in 3149 unique pairings) successfully merged out of 88514 (in 12054 pairings) input.


La commande mergePairs permet la formation de contigs seulement quand
cela est possible. Le read 1 fait 240 pb et le read 2 fait 160 pb donc
nous avons un chevauchement entre les 2 séquences permettant ainsi la
formation des contigs.

# Construire une table de séquence

``` r
# Fait une table de sequence et l'affiche
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]    11 19426


On a créer un objet seqtable et de dans on y met une matrice
d’observation de l’objet mergers grâce a la fonction
makeSequenceTable. la fonction dim permet d’avoir la dimension du
tableau. Le nombre 11 correspond aux lignes et le nombre 19426
correspond aux colonnes.

``` r
# Inspecte la distribution des longueurs de séquence
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  352  353  362  363  364  365  366  367  368  369  370  371  372  373  374  375 
    ##    1    1    1    1    4  183   27  165  184 5608 3594 2312 2613 2738  126 1770 
    ##  376  377  378  382  386 
    ##   90    4    1    1    2


A partir de seqtab on va pouvoir savoir combien de fois on retrouve une
séquence a une certaine longueur en nucléotide. Les reads sont répartis
sur une plage assez resteinte. La majorité des reads ont une longueur de
369 pb.

# Supprimer les chimères

Une séquence chimère est une séquence d’ADN polymérisé par PCR mais qui
n’a pas fini de se polymériser.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 17869 bimeras out of 19426 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   11 1557


Les séquences chimériques doivent être éliminés du jeu de données sinon
cela peut entrainer des erreurs lors de nos analyses. l’objet
seqtab.nochim est créée dans lequel la fonction removeBineraDenovo
permet de supprimer les séquences chimériques. Ici nous avons identifié
17869 chimères sur les 19426 séquences.

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.7769154

``` r
1-sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.2230846

Il y avait 22% des séquences qui était des séquences chimérique, celle
ci ont donc était retiré du jeu de données.

# Suivre les reads dans la pipeline

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.namesfnFs
head(track)
```

    ##                               input filtered denoisedF denoisedR merged nonchim
    ## Station5_Fond1_10sept14_R1   159971   145448    142931    143292 117318   87962
    ## Station5_Fond1_11mars15_R1   175993   160423    158128    158473 138940  111552
    ## Station5_Fond2_10sept14_R1   197039   177018    173601    174591 142188  103668
    ## Station5_Fond2_11mars15_R1    87585    79989     78618     78926  67622   54711
    ## Station5_Fond3_10sept14_R1   117140   106150    103806    104338  83613   64259
    ## Station5_Median1_10sept14_R1 116519   106745    104811    105173  86212   65559


L’objet track correspond aux séquences après chaque étapes d’analyse
réalisé ici. Pour la station5\_fond1\_10\_sept14\_R1 on passe de 159971
à 87962 séquences. On a filtré un grand nombre de séquences.

# Assigniation taxonomique

Nous allons assigner une taxonomie à nos taxons grâce à silva. Cette
assignation taxonomique est déposée dans l’objet taxa.

``` r
# Assigniation Taxonomique
taxa <- assignTaxonomy(seqtab.nochim, "~/cc2_DADA2_import/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum             Class                 Order            
    ## [1,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [2,] "Bacteria" "Cyanobacteria"    "Cyanobacteriia"      "Synechococcales"
    ## [3,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [4,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [5,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [6,] "Bacteria" "Actinobacteriota" "Acidimicrobiia"      "Actinomarinales"
    ##      Family             Genus                    
    ## [1,] "Clade I"          "Clade Ia"               
    ## [2,] "Cyanobiaceae"     "Synechococcus CC9902"   
    ## [3,] "Clade I"          "Clade Ia"               
    ## [4,] "Clade I"          "Clade Ia"               
    ## [5,] "Clade II"         NA                       
    ## [6,] "Actinomarinaceae" "Candidatus Actinomarina"


Ici nous pouvons inspecter les affectations taxonomiques. Les
attributions vont rarement jusqu’à l’espèce car il est souvent
impossible de faire des assignations d’espèces sans ambiguité à partir
d’un fragment du gène 16S. L’assignation s’arrête donc souvent à la
famille et parfois au genre.

``` r
# Sauvegarde des données dans l'environnement afin de les réutiliser
save.image(file="02_analyse_DADA2_FinalEnv")
```
