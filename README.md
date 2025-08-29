# Welcome to glendR
![glendR logo](https://github.com/qwhiting/glendR/blob/main/images/glendR.png)

glendR is an R package of useful functions to format PFAS analysis from SCIEX OS to GLENDA format. Used by the Ulrich Lab at the Natural Resources Research Institute (University of Minnesota). 
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


### Statement on data avalibility and use
glendR is intended for high throughput data formatting in Rstudio. Datasets are required to be formatted and QAQC checked in a certian way for submission to the Great Lakes Environmental Database (GLENDA), where this data will eventually be made freely available to the public. glendR is a tool to be used by the Ulruch Lab at NRRI and was made for a specific funtion (take SCIEX OS data, apply QAQC checks, and convert it to a useable format).
The R functions created in glendR all use the tidyverse package.

**Please cite as**: Quinn Whiting (2025). glendR. R package version 1.0. Gethub repository, https://www.github.com/qwhiting/glendR
