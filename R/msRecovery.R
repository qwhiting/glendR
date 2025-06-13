#this function calculates the matrix spike recovery

#' @import tidyr
#' @import dplyr
#' @export

msRecovery<-function(df){
  #remove "_MS" from sample names
  df$Sample.Name<-gsub("_MS","", df$Sample.Name)
  #get lm samples only
  ms<-df%>%filter(Sample.Type=="LMS")%>%filter(analyteType=="TRG")
  #get the associated sampels (unspiked)
  asms<-df%>%filter(Sample.Name %in% ms$Sample.Name)%>%filter(analyteType=="TRG")%>%filter(Sample.Type!="LMS")
  #create matching IDs to be used later
  ms$tempID<-paste(ms$Sample.Name, ms$Component.Name)
  asms$tempID<-paste(asms$Sample.Name, asms$Component.Name)
  #get the associated concentration (ng/g)
  ms$asconc<-asms$result[match(ms$tempID, asms$tempID)]
  #get the amount recovered from spike (subtract out ng/g of unspiked)
  ms$spikeresult<-ms$result-ms$asconc
  #calculate recovery
  ms$rec<-round(((ms$spikeresult*ms$weightVolumneAnalyzed)/ms$spikeLevel)*100, digits = 1)
  df$tempID<-paste(df$Sample.Name, df$Component.Name)
  #bind to original df
  df$MSrec<-ms$rec[match(df$tempID, ms$tempID)]
  #only have recovery for LMS sample
  df$MSrec<-if_else(df$Sample.Type=="LMS", df$MSrec, NA)
  df<-df[ , -which(names(df) %in% c("tempID"))]
  return(df)
}
