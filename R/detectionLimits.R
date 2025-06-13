#this is a function to create a data frame of the detection limits calculated from the calibration curve from SCIEX OS


#' @import tidyr
#' @import dplyr

#' @export

detectionLimits<-function(df){
  stds<-df%>%filter(Sample.Type=="Standard")%>%filter(Used=="True")%>%filter(Component.Type=="Quantifiers")
  lod<-stds%>% group_by(Component.Name) %>% dplyr::summarise(LLOD=min(Actual.Concentration), ULOD=max(Actual.Concentration))
  return(lod)
}

