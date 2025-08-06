#this function is to flag the IBs and associated samples
#' @import tidyr
#' @import dplyr
#' @export

instrumentBlanks<-function(df){

  IB2<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    #get into one column
    df$ib<-paste0(df$ib1,df$ib2) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }
  IB3<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  IB4<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime

    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }
    IB5<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
      df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
      #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime
    df$ib5time<-ib.time$IB5$time-df$analysisDateTime

    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    df$ib5<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time>=0, "IB5", NA), NA), NA), NA), NA)
    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4,df$ib5) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  IB6<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime
    df$ib5time<-ib.time$IB5$time-df$analysisDateTime
    df$ib6time<-ib.time$IB6$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    df$ib5<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time>=0, "IB5", NA), NA), NA), NA), NA)
    df$ib6<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time>=0, "IB6", NA), NA), NA), NA), NA), NA)
    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4,df$ib5, df$ib6) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  IB7<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime
    df$ib5time<-ib.time$IB5$time-df$analysisDateTime
    df$ib6time<-ib.time$IB6$time-df$analysisDateTime
    df$ib7time<-ib.time$IB7$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    df$ib5<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time>=0, "IB5", NA), NA), NA), NA), NA)
    df$ib6<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time>=0, "IB6", NA), NA), NA), NA), NA), NA)
    df$ib7<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time >=0, "IB7", NA), NA), NA), NA), NA), NA), NA)
    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4,df$ib5, df$ib6, df$ib7) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  IB8<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime
    df$ib5time<-ib.time$IB5$time-df$analysisDateTime
    df$ib6time<-ib.time$IB6$time-df$analysisDateTime
    df$ib7time<-ib.time$IB7$time-df$analysisDateTime
    df$ib8time<-ib.time$IB8$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    df$ib5<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time>=0, "IB5", NA), NA), NA), NA), NA)
    df$ib6<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time>=0, "IB6", NA), NA), NA), NA), NA), NA)
    df$ib7<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time >=0, "IB7", NA), NA), NA), NA), NA), NA), NA)
    df$ib8<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time <0, if_else(df$ib8time>=0, "IB8", NA), NA), NA), NA), NA), NA), NA),NA)

    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4,df$ib5, df$ib6, df$ib7, df$ib8) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  IB9<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime
    df$ib5time<-ib.time$IB5$time-df$analysisDateTime
    df$ib6time<-ib.time$IB6$time-df$analysisDateTime
    df$ib7time<-ib.time$IB7$time-df$analysisDateTime
    df$ib8time<-ib.time$IB8$time-df$analysisDateTime
    df$ib9time<-ib.time$IB9$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    df$ib5<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time>=0, "IB5", NA), NA), NA), NA), NA)
    df$ib6<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time>=0, "IB6", NA), NA), NA), NA), NA), NA)
    df$ib7<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time >=0, "IB7", NA), NA), NA), NA), NA), NA), NA)
    df$ib8<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time <0, if_else(df$ib8time>=0, "IB8", NA), NA), NA), NA), NA), NA), NA),NA)
    df$ib9<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time <0, if_else(df$ib8time<0, if_else(df$ib9time>=0, "IB9", NA), NA), NA), NA), NA), NA), NA),NA), NA)

    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4,df$ib5, df$ib6, df$ib7, df$ib8, df$ib9) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  IB10<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the ibs in the main df
    df$ibFlag<-if_else(df$resultQcIdentifier=="CLB", if_else(df$analyteType=="TRG", if_else(df$result>0.5*df$quantificationLimit, "FNB",NA),NA),NA)
    #create df of only the ibs target PFAS
    df.ib<-df%>%filter(resultQcIdentifier=="CLB")
    #get the number of ibs in the batch as a character
    ib.num<-as.character(n_distinct(df.ib$Sample.Name))
    #correlate the flagged ibs to the associated sampels
    #get the inj. time for each ib and make a new columns of them
    ib.time<-df.ib%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to ib#
    ib.time$Sample.Name<-paste0("IB", row.names(ib.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLB", ib.time$Sample.Name[match(df$analysisDateTime, ib.time$time)], df$Sample.Name)
    #make into a list
    ib.time<-split(ib.time, ib.time$Sample.Name)
    #get the time differece of injections
    df$ib1time<-ib.time$IB1$time-df$analysisDateTime
    df$ib2time<-ib.time$IB2$time-df$analysisDateTime
    df$ib3time<-ib.time$IB3$time-df$analysisDateTime
    df$ib4time<-ib.time$IB4$time-df$analysisDateTime
    df$ib5time<-ib.time$IB5$time-df$analysisDateTime
    df$ib6time<-ib.time$IB6$time-df$analysisDateTime
    df$ib7time<-ib.time$IB7$time-df$analysisDateTime
    df$ib8time<-ib.time$IB8$time-df$analysisDateTime
    df$ib9time<-ib.time$IB9$time-df$analysisDateTime
    df$ib10time<-ib.time$IB10$time-df$analysisDateTime
    #GET ASSOCIATED IBS
    df$ib1<- if_else(df$ib1time>=0, "IB1", NA)
    df$ib2<- if_else(df$ib1time<0, if_else(df$ib2time>=0, "IB2", NA), NA)
    df$ib3<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time>=0, "IB3", NA), NA), NA)
    df$ib4<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time>=0, "IB4", NA), NA), NA), NA)
    df$ib5<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time>=0, "IB5", NA), NA), NA), NA), NA)
    df$ib6<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time>=0, "IB6", NA), NA), NA), NA), NA), NA)
    df$ib7<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time >=0, "IB7", NA), NA), NA), NA), NA), NA), NA)
    df$ib8<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time <0, if_else(df$ib8time>=0, "IB8", NA), NA), NA), NA), NA), NA), NA),NA)
    df$ib9<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time <0, if_else(df$ib8time<0, if_else(df$ib9time>=0, "IB9", NA), NA), NA), NA), NA), NA), NA),NA), NA)
    df$ib10<- if_else(df$ib1time<0, if_else(df$ib2time<0, if_else(df$ib3time<0, if_else(df$ib4time<0, if_else(df$ib5time<0, if_else(df$ib6time<0, if_else(df$ib7time <0, if_else(df$ib8time<0, if_else(df$ib9time<0, if_else(df$ib10time>=0, "IB10", NA), NA), NA), NA), NA), NA), NA),NA), NA),NA)

    #get into one column
    df$ib<-paste0(df$ib1,df$ib2,df$ib3,df$ib4,df$ib5, df$ib6, df$ib7, df$ib8, df$ib9, df$ib10) #get into single column
    df$ib<-gsub("NA","", df$ib) #remove NAs
    df$ib<-if_else(df$ib=="", NA, df$ib) #create NA for blanks
    df$ib<-if_else(is.na(df$ib), df$Sample.Name, df$ib)
    #create matching IDs
    df$ibID<-paste0(df$ib, df$Component.Name)
    df.ib$ibID<-paste0(df.ib$Sample.Name, df.ib$Component.Name)
    #associate
    df$ibFlag<-df.ib$ibFlag[match(df$ibID, df.ib$ibID)]
    df$ibResult<-df.ib$result[match(df$ibID, df.ib$ibID)]
    #remove falg if >5* the concentration in the blank
    df$ibFlag<-if_else(df$result*df$weightVolumeAnalyzed>5*df$ibResult, NA, df$ibFlag)
    #associate the non-IB samples with DIB flags
    df$dibibFlag<-df$ibFlag
    df$dibibFlag<-if_else(df$resultQcIdentifier=="CLB", NA, gsub("FNB", "DIB", df$dibibFlag))
    #add comment
    df$ibComment<-if_else(df$dibibFlag=='DIB', "Detected in Instrument Blank >0.5MDL. Sample result < 5x blank result; therefore result considered non-detect.", NA)
    #remove columns
    remove.cols<-c("ibID", "ib","ib1","ib2","ib3","ib4","ib5","ib6", "ib7","ib8","ib9","ib10","ib1time","ib2time","ib3time","ib4time","ib5time","ib6time", "ib7time","ib8time","ib9time","ib10time", "ibResult")
    df<-df[,!names(df) %in% remove.cols]
    return(df)
  }

  #this picks the correct function above to run based on the unique IBs
  numberOfIb<-df%>%filter(resultQcIdentifier=="CLB")
  numberOfIb<-as.character(n_distinct(numberOfIb$Sample.Name))
  # run each function according to how many unique Ibs there are (max 10 Ibs)
  if(numberOfIb=="1"){
    df<-IB1(df)
  } else if (numberOfIb=="2"){
    df<-IB2(df)
  } else if (numberOfIb=="3"){
    df<-IB3(df)
  } else if (numberOfIb=="4"){
    df<-IB4(df)
  } else if (numberOfIb=="5"){
    df<-IB5(df)
  } else if (numberOfIb=="6"){
    df<-IB6(df)
  } else if (numberOfIb=="7"){
    df<-IB7(df)
  } else if (numberOfIb=="8"){
    df<-IB8(df)
  } else if (numberOfIb=="9"){
    df<-IB9(df)
  } else{
    df<-IB10(df)
  }


 return(df)

}
