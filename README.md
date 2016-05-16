xcms for Galaxy
===============

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io) [![Build Status](https://travis-ci.org/workflow4metabolomics/xcms.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/xcms)

Our project
-----------
The [Workflow4Metabolomics](http://workflow4metabolomics.org), W4M in short, is a French infrastructure offering software tool processing, analyzing and annotating metabolomics data. It is based on the Galaxy platform.


xcms
----
LC/MS and GC/MS Data Analysis

Bioconductor version: Release (3.2)

Framework for processing and visualization of chromatographically separated and single-spectra mass spectral data. Imports from AIA/ANDI NetCDF, mzXML, mzData and mzML files. Preprocesses data for high-throughput, untargeted analyte profiling.

Author: Colin A. Smith <csmith at scripps.edu>, Ralf Tautenhahn <rtautenh at gmail.com>, Steffen Neumann <sneumann at ipb-halle.de>, Paul Benton <hpbenton at scripps.edu>, Christopher Conley <cjconley at ucdavis.edu>, Johannes Rainer <Johannes.Rainer at eurac.edu>

Maintainer: Steffen Neumann <sneumann at ipb-halle.de>

Citation (from within R, enter citation("xcms")):

Smith, C.A., Want, E.J., O'Maille, G., Abagyan,R., Siuzdak and G. (2006). “XCMS: Processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching and identification.” Analytical Chemistry, 78, pp. 779–787.

Tautenhahn R, Boettcher C and Neumann S (2008). “Highly sensitive feature detection for high resolution LC/MS.” BMC Bioinformatics, 9, pp. 504.

Benton HP, Want EJ and Ebbels TMD (2010). “Correction of mass calibration gaps in liquid chromatography-mass spectrometry metabolomics data.” BIOINFORMATICS, 26, pp. 2488.

Homepage: [https://bioconductor.org/packages/release/bioc/html/xcms.html](https://bioconductor.org/packages/release/bioc/html/xcms.html)


Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)


Dependencies using Conda
------------------------
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)

[Conda](http://conda.pydata.org/) is package manager that among many other things can be used to manage Python packages.

The main recipe: [https://github.com/bioconda/bioconda-recipes/tree/master/recipes/bioconductor-xcms](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/bioconductor-xcms)

```
#To install miniconda2
#http://conda.pydata.org/miniconda.html
#To install the xcms R library using conda:
conda install r-snow bioconductor-xcms r-batch
#To set an environment:
conda create -n bioconductor-xcms r-snow bioconductor-xcms r-batch`
#To activate the environment:
. activate bioconductor-xcms
```

Travis
------
[![Build Status](https://travis-ci.org/workflow4metabolomics/xcms.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/xcms)

Test and Deploy with Confidence. Easily sync your GitHub projects with Travis CI and you'll be testing your code in minutes!


Test Status
-----------

Planemo test using conda: passed

Planemo test using source env.sh: passed

Planemo shed_test: passed


Historic contributors
---------------------
 - Gildas Le Corguillé @lecorguille - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [UPMC](www.upmc.fr)/[CNRS](www.cnrs.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
 - Misharl Monsoor @mmonsoor - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[UPMC](www.upmc.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
 - Yann Guitton @yguitton - [LABERCA - Laboratory of Food Contaminants and Residue Analysis](http://www.laberca.org/) - Ecole Nationale Vétérinaire, Agroalimentaire et de l'Alimentation Nantes-Atlantique - France
 - Urszula Czerwinska [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[UPMC](www.upmc.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
 - Pierre Pericard - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [CNRS](www.cnrs.fr)/[UPMC](www.upmc.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
