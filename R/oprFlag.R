#this function is to flag and associate OPR samples

#' @import tidyr
#' @import dplyr
#' @export

oprFlag<-function(df){
  df.opr<-df%>%filter(resultQcIdentifier=="OPR")
  df.opr$batch<-if_else(str_detect(df.opr$Sample.Name, "1"), "b1", "b2")
  opr.time<-df.opr%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
  #get OPR flags
  df.opr$oprFlag<-if_else(df.opr$analyteType=="TRG", if_else(df.opr$recovery>40, if_else(df.opr$recovery<160, NA, "FLS"), "FLS"), NA)
  df.opr$hiFlag<-if_else(df.opr$analyteType=='TRG', if_else(df.opr$recovery>160, "HIB",NA),NA)
  df.opr$lowFlag<-if_else(df.opr$analyteType=='TRG', if_else(df.opr$recovery<40, "LOB",NA),NA)
  #get matching IDs
  df.opr$oprID<-paste(df.opr$batch, df.opr$Component.Name)
  #associate samples to OPRs
  ## this gets the times for the OPR injections, everything before b2 is batch1 the rest are in batch 2
  b2cutoff<-df.opr%>%group_by(batch)%>%summarise(batchTime=min(analysisDateTime, na.rm=T))
  b2cutoff<-b2cutoff%>%filter(batch=="b2")
  #assocate the batch
  df$oprBatch<-if_else(df$analysisDateTime>=b2cutoff$batchTime, "b2", "b1")
  #assocate the analyte
  assOPR<-df.opr%>%group_by(oprID)%>%reframe(flag=oprFlag, hi=hiFlag, low=lowFlag)%>%distinct(.)

  #get batchID
  df$batchID<-paste(df$oprBatch, df$Component.Name)
  #associate to all samples
  df$oprFlag<-assOPR$flag[match(df$batchID, assOPR$oprID)]
  df$oprHi<-assOPR$hi[match(df$batchID, assOPR$oprID)]
  df$oprLo<-assOPR$low[match(df$batchID, assOPR$oprID)]
 #remove flags from CAL/CLC/CLB/FVM samples
  df$oprFlag<-if_else(df$resultQcIdentifier%in%c("CLC","CAL","CLB","LVM"), NA, df$oprFlag)

#remove flag from OPR sample if it is from the associated OPR (as they should not be associated)
  df$oprFlagKeep<-if_else(df$resultQcIdentifier=="OPR", if_else(df$analyteType=="TRG", if_else(df$recovery>40, if_else(df$recovery<160, NA, "FLS"), "FLS"), NA),NA)
  df$oprFlag<-if_else(df$resultQcIdentifier=="OPR", df$oprFlagKeep, df$oprFlag)
  df$oprHi<-if_else(is.na(df$oprFlag), NA, df$oprHi)
  df$oprLo<-if_else(is.na(df$oprFlag), NA, df$oprLo)
  #remove Hi/Lo flags from OPR samples
  df$oprHi<-if_else(df$resultQcIdentifier=="OPR", NA, df$oprHi)
  df$oprLo<-if_else(df$resultQcIdentifier=="OPR", NA, df$oprLo)
  #remove HIB if <MDL
  df$oprHi<-if_else(df$result<=df$mdl, NA, df$oprHi)
  #rename and remove columns
  colnames(df)[colnames(df)=="oprFlag"]<-"fls"
  colnames(df)[colnames(df)=="oprHi"]<-"flsHi"
  colnames(df)[colnames(df)=="oprLo"]<-"flsLo"

  col.remove<-c("oprBatch", "batchID", "oprFlagKeep")
  df<-df[,!names(df) %in% col.remove]

  return(df)
}
