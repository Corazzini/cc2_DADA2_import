Install packages
================

  - [mettre à jour la configuration de la
    VM](#mettre-à-jour-la-configuration-de-la-vm)
  - [installation du package](#installation-du-package)
  - [knitr](#knitr)

# mettre à jour la configuration de la VM

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install -y libglpk-dev
```

    ## sudo: unable to resolve host dc0402adb105: Name or service not known
    ## Hit:1 http://archive.ubuntu.com/ubuntu focal InRelease
    ## Get:2 http://archive.ubuntu.com/ubuntu focal-updates InRelease [114 kB]
    ## Get:3 http://archive.ubuntu.com/ubuntu focal-backports InRelease [101 kB]
    ## Get:4 http://security.ubuntu.com/ubuntu focal-security InRelease [109 kB]
    ## Get:5 http://archive.ubuntu.com/ubuntu focal-updates/universe amd64 Packages [881 kB]
    ## Get:6 http://archive.ubuntu.com/ubuntu focal-updates/main amd64 Packages [885 kB]
    ## Fetched 2,090 kB in 1s (2,575 kB/s)
    ## Reading package lists...
    ## sudo: unable to resolve host dc0402adb105: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libbz2-dev is already the newest version (1.0.8-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 47 not upgraded.
    ## sudo: unable to resolve host dc0402adb105: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## liblzma-dev is already the newest version (5.2.4-1ubuntu1).
    ## 0 upgraded, 0 newly installed, 0 to remove and 47 not upgraded.
    ## sudo: unable to resolve host dc0402adb105: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libglpk-dev is already the newest version (4.65-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 47 not upgraded.

# installation du package

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.11')
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("dada2", version = "3.11")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'dada2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("phangorn")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phangorn'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("DECIPHER")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DECIPHER'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("DESseq2")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DESseq2'

    ## Warning: package 'DESseq2' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("ggplot2")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'ggplot2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("phyloseq")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phyloseq'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

``` r
BiocManager::install("Biostrings", version = "3.11")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'Biostrings'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

    ## Old packages: 'broom', 'pkgbuild', 'RMariaDB', 'testthat', 'vctrs'

\#Phyloseq

``` r
install.packages ("gridExtra")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
install.packages('grid.arrange')
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

    ## Warning: package 'grid.arrange' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

# knitr

``` r
install.packages ("knitr")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)
