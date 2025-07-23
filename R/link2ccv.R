# this function is to deal with CCVs and associate them to their respective samples


#' @import tidyr
#' @import dplyr
#' @import lubridate
#' @export

link2ccv<-function(df){


  CCV1<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1)

    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"
  return(df)
  }

  CCV2<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1,df$ccv2)
    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

    return(df)
  }

  CCV3<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
    df$ccv3time<-ccv.time$CCV3$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
    df$ccv3<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time>=0, "CCV3", NA), NA), NA)
    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1,df$ccv2,df$ccv3)
    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

    return(df)
  }

  CCV4<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
    df$ccv3time<-ccv.time$CCV3$time-df$analysisDateTime
    df$ccv4time<-ccv.time$CCV4$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
    df$ccv3<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time>=0, "CCV3", NA), NA), NA)
    df$ccv4<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time>=0, "CCV4", NA), NA), NA), NA)
    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1,df$ccv2,df$ccv3,df$ccv4)
    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

    return(df)
  }

  CCV5<-function(df){
  #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
  df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
  #create df of only the CCVs target PFAS
  df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
  #get the number of CCVs in the batch as a character
  cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
  #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
  ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
  #change the sample name to CCV#
  ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
  #add the changed names back into original df
  df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
  #make into a list
  ccv.time<-split(ccv.time, ccv.time$Sample.Name)
  ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
  # get the time difference between CCV and sample injections

  df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
  df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
  df$ccv3time<-ccv.time$CCV3$time-df$analysisDateTime
  df$ccv4time<-ccv.time$CCV4$time-df$analysisDateTime
  df$ccv5time<-ccv.time$CCV5$time-df$analysisDateTime
  #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
  df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
  df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
  df$ccv3<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time>=0, "CCV3", NA), NA), NA)
  df$ccv4<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time>=0, "CCV4", NA), NA), NA), NA)
  df$ccv5<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time>=0, "CCV5", NA), NA), NA), NA), NA)
  #get ccvs into a single column
  df$ccv<-paste0(df$ccv1,df$ccv2,df$ccv3,df$ccv4, df$ccv5)

  df$ccv<-gsub("NA","", df$ccv)
  df$ccv<-if_else(df$ccv=="", NA, df$ccv)
  df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
  #remove CALs from ccv
  df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
  #create a matching ID
  df$ccvID<-paste0(df$ccv, df$Component.Name)
  df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
  #associate
  df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
  # remove flags for samples <MDL
  df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
  #columns to remove
  col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
  #remove columns
  df<-df[,!names(df) %in% col.remove]
  #rename flag column
  colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

  return(df)
  }

  CCV6<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
    df$ccv3time<-ccv.time$CCV3$time-df$analysisDateTime
    df$ccv4time<-ccv.time$CCV4$time-df$analysisDateTime
    df$ccv5time<-ccv.time$CCV5$time-df$analysisDateTime
    df$ccv6time<-ccv.time$CCV6$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
    df$ccv3<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time>=0, "CCV3", NA), NA), NA)
    df$ccv4<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time>=0, "CCV4", NA), NA), NA), NA)
    df$ccv5<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time>=0, "CCV5", NA), NA), NA), NA), NA)
    df$ccv6<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time<0, if_else(df$ccv6time>=0, "CCV6", NA), NA), NA), NA), NA),NA)

    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1,df$ccv2,df$ccv3,df$ccv4, df$ccv5, df$ccv6)
    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

    return(df)
  }

  CCV7<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
    df$ccv3time<-ccv.time$CCV3$time-df$analysisDateTime
    df$ccv4time<-ccv.time$CCV4$time-df$analysisDateTime
    df$ccv5time<-ccv.time$CCV5$time-df$analysisDateTime
    df$ccv6time<-ccv.time$CCV6$time-df$analysisDateTime
    df$ccv7time<-ccv.time$CCV7$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
    df$ccv3<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time>=0, "CCV3", NA), NA), NA)
    df$ccv4<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time>=0, "CCV4", NA), NA), NA), NA)
    df$ccv5<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time>=0, "CCV5", NA), NA), NA), NA), NA)
    df$ccv6<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time<0, if_else(df$ccv6time>=0, "CCV6", NA), NA), NA), NA), NA),NA)
    df$ccv7<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time<0, if_else(df$ccv6time<0, if_else(df$ccv7time>=0, "CCV7", NA), NA), NA), NA), NA),NA),NA)
    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1,df$ccv2,df$ccv3,df$ccv4, df$ccv5, df$ccv6, df$ccv7)
    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

    return(df)
  }

  CCV8<-function(df){
    #flag the TRG compounds outside the 70-130% recovery range in the CCVs in the main df
    df$ccvFlag<-if_else(df$resultQcIdentifier=="CLC", if_else(df$analyteType=="TRG", if_else(df$recovery>70, if_else(df$recovery<130, NA, "FCC"),"FCC"),NA),NA)
    #create df of only the CCVs target PFAS
    df.ccv<-df%>%filter(resultQcIdentifier=="CLC")
    #get the number of CCVs in the batch as a character
    cvv.num<-as.character(n_distinct(df.ccv$Sample.Name))
    #correlate the flagged ccvs to the associated sampels
    #get the inj. time for each CCV and make a new columns of them
    ccv.time<-df.ccv%>%group_by(Sample.Name)%>%summarise(time=unique(analysisDateTime))
    #change the sample name to CCV#
    ccv.time$Sample.Name<-paste0("CCV", row.names(ccv.time))
    #add the changed names back into original df
    df$Sample.Name<-if_else(df$resultQcIdentifier=="CLC", ccv.time$Sample.Name[match(df$analysisDateTime, ccv.time$time)], df$Sample.Name)
    #make into a list
    ccv.time<-split(ccv.time, ccv.time$Sample.Name)
    ## this depends on the total CCVs in the df, an if_else statement allows to have different numbers of CCVs in the function
    # get the time difference between CCV and sample injections

    df$ccv1time<-ccv.time$CCV1$time-df$analysisDateTime
    df$ccv2time<-ccv.time$CCV2$time-df$analysisDateTime
    df$ccv3time<-ccv.time$CCV3$time-df$analysisDateTime
    df$ccv4time<-ccv.time$CCV4$time-df$analysisDateTime
    df$ccv5time<-ccv.time$CCV5$time-df$analysisDateTime
    df$ccv6time<-ccv.time$CCV6$time-df$analysisDateTime
    df$ccv7time<-ccv.time$CCV7$time-df$analysisDateTime
    df$ccv8time<-ccv.time$CCV8$time-df$analysisDateTime
    #associate the CCVs (including the CCVs themselves [i.e., CCV1 is associated with CCV1])
    df$ccv1<- if_else(df$ccv1time>=0, "CCV1", NA)
    df$ccv2<- if_else(df$ccv1time<0, if_else(df$ccv2time>=0, "CCV2", NA), NA)
    df$ccv3<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time>=0, "CCV3", NA), NA), NA)
    df$ccv4<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time>=0, "CCV4", NA), NA), NA), NA)
    df$ccv5<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time>=0, "CCV5", NA), NA), NA), NA), NA)
    df$ccv6<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time<0, if_else(df$ccv6time>=0, "CCV6", NA), NA), NA), NA), NA),NA)
    df$ccv7<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time<0, if_else(df$ccv6time<0, if_else(df$ccv7time>=0, "CCV7", NA), NA), NA), NA), NA),NA),NA)
    df$ccv8<- if_else(df$ccv1time<0, if_else(df$ccv2time<0, if_else(df$ccv3time<0, if_else(df$ccv4time<0, if_else(df$ccv5time<0, if_else(df$ccv6time<0, if_else(df$ccv7time<0, if_else(df$ccv8time>=0, "CCV8", NA), NA), NA), NA), NA),NA),NA),NA)
    #get ccvs into a single column
    df$ccv<-paste0(df$ccv1,df$ccv2,df$ccv3,df$ccv4, df$ccv5, df$ccv6, df$ccv7, df$ccv8)
    df$ccv<-gsub("NA","", df$ccv)
    df$ccv<-if_else(df$ccv=="", NA, df$ccv)
    df$ccv<-if_else(is.na(df$ccv), df$Sample.Name, df$ccv)
    #remove CALs from ccv
    df$ccv<-if_else(df$resultQcIdentifier %in% c("CAL", "LVM"), NA, df$ccv)
    #create a matching ID
    df$ccvID<-paste0(df$ccv, df$Component.Name)
    df.ccv$ccvID<-paste0(df.ccv$Sample.Name, df.ccv$Component.Name)
    #associate
    df$assCCVflag<-df.ccv$ccvFlag[match(df$ccvID, df.ccv$ccvID)]
    # remove flags for samples <MDL
    df$assCCVflag2<-if_else(df$result<df$mdl, NA, df$assCCVflag)
    #columns to remove
    col.remove<-c("ccv1time","ccv2time","ccv3time","ccv4time","ccv5time","ccv6time","ccv7time","ccv8time", "ccv1","ccv2","ccv3","ccv4","ccv5","ccv6","ccv7","ccv8", "ccv","ccvID", "assCCVflag", "ccvFlag")
    #remove columns
    df<-df[,!names(df) %in% col.remove]
    #rename flag column
    colnames(df)[which(colnames(df)=="assCCVflag2")]<-"ccvFlag"

    return(df)
  }

  numberOfCcv<-df%>%filter(resultQcIdentifier=="CLC")
  numberOfCcv<-as.character(n_distinct(numberOfCcv$Sample.Name))
 # run each function according to how many unique CCVs there are (max 8 CCVs)
  if(numberOfCcv=="1"){
    df<-CCV1(df)
  } else if (numberOfCcv=="2"){
    df<-CCV2(df)
  } else if (numberOfCcv=="3"){
    df<-CCV3(df)
  } else if (numberOfCcv=="4"){
    df<-CCV4(df)
  } else if (numberOfCcv=="5"){
    df<-CCV5(df)
  } else if (numberOfCcv=="6"){
    df<-CCV6(df)
  } else if (numberOfCcv=="7"){
    df<-CCV7(Df)
  } else{
    df<-CCV8(df)
    }

 return(df)

  #df<-if_else(numberOfCcv=="1", CCV1, if_else(numberOfCcv=="2", CCV2, if_else(numberOfCcv=="3", CCV3, if_else(numberOfCcv=="4", CCV4, if_else(numberOfCcv=="5", CCV5, if_else(numberOfCcv=="6", CCV6, if_else(numberOfCcv=="7", CCV7, CCV8)))))))
  #return(df)
}


