#this function is to flag and associate CAL samples

#' @import tidyr
#' @import dplyr
#' @export

calFlag<-function(df){
  df.cal<-df%>%filter(resultQcIdentifier=="CAL")%>%filter(analyteType=="TRG")
  df.cal$calFlag<-if_else(df.cal$recovery <70 | df.cal$recovery>130, "FVM", NA)
  #remove flag if cal is below ICAL LLOD
  df.cal$calFlag<-if_else(str_detect(df.cal$limFlag, "LTL"), NA, df.cal$calFlag)
  #match to all other relevant samples (exclude LVM, CLC)
  df$calFlag<-df.cal$calFlag[match(df$Component.Name, df.cal$Component.Name)]
  df$calFlag<-if_else(df$resultQcIdentifier %in% c("LVM","CLC"), NA, df$calFlag)
  #remove flag if <MDL, <LTL, or UND
  df$calFlag<-if_else(df$resultQcIdentifier=="CAL", df$calFlag, if_else(str_detect(df$limFlag, "UND"), NA, df$calFlag))
  df$calFlag<-if_else(df$resultQcIdentifier=="CAL", df$calFlag, if_else(str_detect(df$limFlag, "LTL"), NA, df$calFlag))
  df$calFlag<-if_else(df$resultQcIdentifier=="CAL", df$calFlag, if_else(str_detect(df$limFlag, "MDL"), NA, df$calFlag))

  return(df)
}
