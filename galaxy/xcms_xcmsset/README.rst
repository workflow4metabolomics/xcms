
Changelog/News
--------------

**Version 2.0.8 - 06/04/2016**

- TEST: refactoring to pass planemo test using conda dependencies


**Version 2.0.7 - 10/02/2016**

- BUGFIX: better management of errors. Datasets remained green although the process failed

- BUGFIX/IMPROVEMENT: New checking steps around the imported data in order to raise explicte error message before or after launch XCMS: checking of bad characters in the filenames, checking of the XML integrity and checking of duplicates which can appear in the sample names during the XCMS process because of bad characters

- BUGFIX/IMPROVEMENT: New step to check and delete bad characters in the XML: accented characters in the storage path of the mass spectrometer

- UPDATE: refactoring of internal management of inputs/outputs

- UPDATE: refactoring to feed the new report tool


**Version 2.0.2 - 18/01/2016**

- BUGFIX: Some zip files were tag as "corrupt" by R. We have changed the extraction mode to deal with thoses cases.


**Version 2.0.2 - 09/10/2015**

- BUGFIX: Some users reported a bug in xcms.xcmsSet. The preprocessing stops itself and doesn't import the whole dataset contained in the zip file without warning. But meanwhile, please check your samplemetadata dataset and the number of rows.


**Version 2.0.2 - 02/06/2015**

- NEW: The W4M workflows will now take as input a zip file to ease the transfer and to improve dataset exchange between tools and users. (See How_to_upload). The previous "Library directory name" is still available but we invite user to switch on the new zip system as soon as possible.

- IMPROVEMENT: new datatype/dataset formats (rdata.xcms.raw, rdata.xcms.group, rdata.xcms.retcor ...) will facilitate the sequence of tools and so avoid incompatibility errors.

- IMPROVEMENT: parameter labels have changed to facilitate their reading.


Test Status
-----------

Planemo test using conda: passed

