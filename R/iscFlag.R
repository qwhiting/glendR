# this function is used to flag the ISC and any associated samples


#' @import tidyr
#' @import dplyr
#' @export

iscFlag<-function(df){
  df.isc<-df%>%filter(resultQcIdentifier=="LVM")%>%filter(analyteType=="TRG")
  #df.isc$iscFlag<-if_else(df.isc$recovery<70 | df.isc$recovery>130, "LVM", NA)
  df.isc$iscFlag<-if_else(df.isc$recovery<70, "LVM", NA)
  df.isc$iscFlag<-if_else(df.isc$recovery>130, "LVM", df.isc$iscFlag)

  #remove flag if isc is below ICAL LLOD
  df.isc$iscFlag<-if_else(str_detect(df.isc$limFlag, "LTL"), NA, df.isc$iscFlag)
  #match to all other relevant samples (exclude CAL, CLC) the ISC corresponds to entire batch
  df$iscFlag<-df.isc$iscFlag[match(df$Component.Name, df.isc$Component.Name)]
  df$iscFlag<-if_else(df$resultQcIdentifier %in% c("CAL","CLC"), NA, df$iscFlag)
  #remove flag if <MDL, <LTL, or UND
  df$iscFlag<-if_else(df$resultQcIdentifier=="LVM", df$iscFlag, if_else(str_detect(df$limFlag, "UND"), NA, df$iscFlag))
  df$iscFlag<-if_else(df$resultQcIdentifier=="LVM", df$iscFlag, if_else(str_detect(df$limFlag, "LTL"), NA, df$iscFlag))
  df$iscFlag<-if_else(df$resultQcIdentifier=="LVM", df$iscFlag, if_else(str_detect(df$limFlag, "MDL"), NA, df$iscFlag))

  return(df)
}
