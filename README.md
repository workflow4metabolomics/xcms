NMR preprocessing for Galaxy
========================

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io) [![Build Status](https://travis-ci.org/workflow4metabolomics/nmr_preprocessing.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/nmr_preprocessing)

Our project
-----------
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data. It is based on the Galaxy platform.


NMR Preprocessing
-----------------

Preprocessing of 1D NMR spectra

Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)


Dependencies using Conda
------------------------
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)

The main recipe: [https://github.com/bioconda/bioconda-recipes/tree/master/recipes/r-ptw](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/r-ptw)

```
#To install miniconda2
#http://conda.pydata.org/miniconda.html
#To install the needed R library using conda:
conda install r-batch r-ptw
#To set an environment:
conda create -n nmr_preprocessing r-batch r-ptw r-Matrix r-ggplot2 r-gridExtra r-reshape2`

#To activate the environment:
. activate nmr_preprocessing
```

[Conda](http://conda.pydata.org/) is package manager that among many other things can be used to manage Python packages.

Travis
------
[![Build Status](https://travis-ci.org/workflow4metabolomics/nmr_preprocessing.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/nmr_preprocessing)

Test and Deploy with Confidence. Easily sync your GitHub projects with Travis CI and you'll be testing your code in minutes!

Historic contributors
---------------------
 - Manon Martin @mmartin - [Université catholique de Louvain (Belgique)]
 - Marie Tremblay-Franco @mtremblayfr - [French Metabolomics and Fluxomics Infrastructure (MetaboHUB)](http://www.metabohub.fr/en) - [MetaToul](http://www.metatoul.fr/)
 
