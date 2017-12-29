Instructions to integrate the ""NMR Annotation" tool into a local instance of Galaxy
Version December 2017 M Tremblay-Franco


## --- R bin and Packages : --- ##
R version 3.0.2 (2013-09-25) -- "Frisbee Sailing
Platform: x86_64-redhat-linux-gnu (64-bit)

Install the "batch" library, necessary for parseCommandArgs function and the "ASICS" library necessary for NMR spectra annotation:
 - Download package source (*.tar.gz file) from your favorite CRAN (http://www.r-project.org/)
For example: http://cran.univ-lyon1.fr/

 - Install package in your R session
install.packages("path/package_name.tar.gz",lib="path",repos=NULL)
For Example: install.packages("/usr/lib64/R/library/batch_1.1-4.tar",lib="/usr/lib64/R/library",repos=NULL)

 - Finally, load the packages into your R session
library(batch)
library(ASICS)


## --- Config : --- ##
 - Edit the file "/galaxy/dist/galaxy-dist/tool_conf.xml" and add 
<section id="id_name" name="Name">
  <tool file="path/asics_xml.xml" />
</section>
to create a new section containing the asics_xml tool
or add
  <tool file="path/asics_xml.xml" />
in an existing section

 - Put the two files asics_xml.xml, asics_wrapper.R, all the needeed R functions and the Library.RData (including compound reference spectra) in a same directory
For example, path=/galaxy/dist/galaxy-dist/tools/annotation

Finally, restart Galaxy