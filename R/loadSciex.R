# This function loads the exported .txt file from the processed SCIEX OS
# Required to have 'utils' and 'dplyr' packages installed, 'tidyverse' loads these automatically
#--#
#roxygen2 to create NAMESPACE- import required libraries for function, add export at end
#--#

#' @import utils
#' @import tidyr
#' @import dplyr
#' @export

loadSciex<-function(file.path){
  df<-read.table(file = file.path, header=T, fill = T, sep = "\t")
  df<-df%>%mutate_at(c("Actual.Concentration", "IS.Actual.Concentration", "Area", "IS.Area", "Area.Ratio", "Retention.Time","Signal...Noise","Calculated.Concentration", "Accuracy", "Ion.Ratio", "Concentration.Ratio"),as.numeric)
  return(df)
}
