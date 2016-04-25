
Changelog/News
--------------

**Version 2.1.5 - 21/04/2016**

- UPGRADE: upgrate the CAMERA version from 1.22.0 to 1.26.0


**Version 2.1.4 - 18/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies

**Version 2.1.3 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- BUGFIX: the conversion into minutes of the retention time was applied to the diffreport outputs (several conditions)

- IMPROVEMENT: when there are several conditions, the tool will generate individual datasets (tsv, pdf) instead of a zip file. The usual png (eic, boxplot) will from now be integrated in two pdf.

- UPDATE: refactoring of internal management of inputs/outputs


VERSION 2.1.0 - 09/10/2015**

- BUGFIX: There was a bug with the CAMERA.annotate (generating a bad dataMatrix (intensities which don't match with the metabolites))


VERSION 2.1.0 - 07/06/2015**

- IMPROVEMENT: new datatype/dataset formats (rdata.camera.positive, rdata.camera.negative, rdata.camera.quick ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.

- UPDATE: merged with annotateDiffreport. Some parameters are dedicated to experiences with several conditions



Test Status
-----------

Planemo test using conda: passed

Planemo test using source env.sh: passed

Planemo shed_test : passed
