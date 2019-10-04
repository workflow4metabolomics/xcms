# massPix for Galaxy #

This is Galaxy tool for R package `massPix` which processes high resolution
mass spectrometry imaging data, performs multivariate statistics (PCA,
clustering) and lipid identification. For details, see
[R package massPix](https://github.com/hallz/massPix) and [Metabolomics
publication](https://link.springer.com/article/10.1007%2Fs11306-017-1252-5).

## Installation ##

- Install [Galaxy](https://github.com/galaxyproject/galaxy) under Linux.

- Install [conda](https://docs.conda.io/en/latest/miniconda.html) under
  Linux. `conda` is used to install `requirements` of this tool. 

- Use `git` to clone this tool

  ```bash
  git clone https://github.com/wanchanglin/massPix.git
  ```

- Add this tool's location into Galaxy' tool config file:
  `~/Galaxy/config/tool_conf.xml`. For example, one simplified
  `tool_conf.xml` looks like:

  ```xml
  <?xml version='1.0' encoding='utf-8'?>
  <toolbox monitor="true">
    
    <section id="getext" name="Get Data">
      <tool file="data_source/upload.xml" />
    </section>
    
    <section id="MyTools" name="My Tools">
      <tool file="/path/to/massPix/massPix.xml" />
    </section>

  </toolbox>
  ```

- The Galaxy tool uses [R](https://cran.r-project.org/) packages `optparse`,
  `WriteXLS`, `calibrate` and  `rJava`. This tool embeds a modified
  `massPix` package for Galaxy only.

- To test the tool, use data inside `test-data` or download test data (.ibd
  and .imzML files) from the
  [MetaboLights](https://www.ebi.ac.uk/metabolights/) repository accession
  number [MTBLS487](https://www.ebi.ac.uk/metabolights/MTBLS487). 

## Authors, contributors & contacts ##

- Wanchang Lin (wl361@cam.ac.uk), University of Cambridge 
- Zoe Hall (zlh22@cam.ac.uk), University of Cambridge 
- Julian L Griffin (jlg40@cam.ac.uk), University of Cambridge 

