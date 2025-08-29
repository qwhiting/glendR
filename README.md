# Welcome to glendR
<img src="https://github.com/qwhiting/glendR/blob/main/images/glendR.png" alt="glendR logo" style="width:30%; height:auto;">

glendR is an R package of useful functions intended to help with:
1) calculating recoveries (surrogate and native)
2) applying instrumental (CCVs, IBs, ISC) and laboratory (MBs, OPRs) QAQCs
3) determining the calibration range and marking which samples require dilution and re-analysis
4) determining the %RPD for duplicate samples
5) determining the matrix spike recoveries
6) appling qualifier codes
7) formating the dataset into the **G**reat **L**akes **EN**vironmental **DA**tabase (GLENDA)

General instructions can be found in the 'GLENDR HOW TO' html file.
## Stating from scratch
To download and use this R package, please follow the instructions:
-open Rstudio
- install the package 'remotes'
- ```r
  install.packages('remotes')
  ```
- use the install_github function to download the most current version of the package
- ```r
  install_github("qwhiting/glendR")
  ```

Once you have the R package downloaded and installed, check that ALL data is in the correct format. This includes the following:
- SCIEX output data (.txt)
- sample info data (.csv)
- MDL data (.csv)
- EIS and NIS data (.csv)
- Spiking data (.csv)

Further details on data format are in [GLENDR HOW TO](https://github.com/qwhiting/glendR/blob/main/GLENDR%20HOW%20TO.html)
Additionally, a video of how to use the glendR package is avalibe upon request or in the shared Box folder for the Ulrich Lab. 


### Statement on data avalibility and use
glendR is intended for high throughput data formatting in Rstudio. Datasets are required to be formatted and QAQC checked in a certian way for submission to GLENDA, where the data will eventually be made freely available to the public. glendR was made to be used as a tool for the Ulruch Lab at the Natural Resources Research Institute (University of Minnesota) and was made for a specific funtion (take SCIEX OS data > apply QAQC checks > convert it to a useable format). Due to the specific use of this package, functions use specific nomenclature that SCIEX OS outputs. Please freely use this package and make edits as your data requires, however please note that using this package ouside of the intended use without making changes to the functions may result errant data.
The R functions created in glendR all use the tidyverse package.

**Please cite as**: Quinn Whiting (2025). glendR. R package version 1.0. Gethub repository, https://www.github.com/qwhiting/glendR
