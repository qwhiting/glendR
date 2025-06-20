# this function is to deal with CCVs and associate them to their respective samples


#' @import tidyr
#' @import dplyr
#' @import lubridate
#' @export

CCVass<-function(df){
  #create df of only the CCVs target PFAS
  df.ccv<-filter(sampleType=="CLC")%>%filter(analyteType=="TRG")
  #get the number of CCVs in the batch
  cvv.num<-count(unique(df.ccv$field_id))

}


