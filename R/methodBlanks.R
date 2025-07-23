#this function associates the method blanks and flags them accordingly

#' @import tidyr
#' @import dplyr
#' @export

methodBlanks<-function(df){
  #flag the MBs
  df$mbFlag<-if_else(df$resultQcIdentifier=="LMB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$mdl, "FPB", NA),NA), NA)
  #create id to match mb to samples with
  df$mbID<-paste(df$samplePrepDate, df$Component.Name)
  #get df of only the flagged MBs
  df.mb<-df%>%filter(resultQcIdentifier=="LMB")%>%filter(mbFlag=="FPB")
  #associate based on mbID and flag corresponding samples
  df$dibFlag<-df.mb$mbFlag[match(df$mbID, df.mb$mbID)]
  df$mbResult<-df.mb$result[match(df$mbID, df.mb$mbID)]
  df$dibmbFlag<-if_else(df$result>5*df$mbResult, NA, df$dibmbFlag)
  df$dibmbFlag<-gsub("FPB", "DIB", df$dibmbFlag)
  df$mbComment<-if_else(df$dibmbFlag=="DIB", "Detected in Method Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
df<-df[,-which(names(df) %in% c("mbResult", "mbID"))]
return(df)
  }
