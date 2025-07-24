# cleans up the df into the GLENDA format

#' @import tidyr
#' @import dplyr
#' @export

glendaFormat<-function(df, glendaInfo){
  columns.to.remove<-c("IS.Name","IS.Actual.Concentration","Actual.Concentration","Calculated.Concentration","Acquisition.Date...Time","Accuracy","Retention.Time","matchID","eisConc","eisRec","nisRec","trgConc","trgLCSrec","matchID2","eisrec2","nis","matchID3","nisrec2", "rpd","ionratio","rtdelta", "MSrec")
  df<-df[,!names(df) %in% columns.to.remove]
  #resultComment
  df$resultComments<-paste(df$rpdCom, df$ibComment, df$mbComment)
  df$resultComments<-gsub("NA", "", df$resultComments)
  df$resultComments<-trimws(df$resultComments)
  df$resultComments<-if_else(df$resultComments=="", NA, df$resultComments)
  #qualifier
  df$qualifier<-paste(df$limFlag, df$eisFlag1, df$eisFlag2, df$eisLo, df$eisHi, df$nisFlag1, df$nisFlag2, df$dupFlag, df$msFlag, df$msLo, df$msHi, df$IRflag, df$rtFlag, df$mbFlag, df$dibmbFlag, df$ibFlag, df$dibibFlag, df$ccvFlag, df$fls, df$flsLo, df$flsHi, df$calFlag, df$iscFlag, sep = ", ")
  df$qualifier<-gsub("NA, ", "", df$qualifier)
  df$qualifier<-gsub(", NA", "", df$qualifier)
  df$qualifier<-gsub("NA", "", df$qualifier)
  df$qualifier<-trimws(df$qualifier)
  df$qualifier<-if_else(df$qualifier=="", NA, df$qualifier)
  #remove columns
  more.cols.to.remove<-c("rpdCom","ibComment","mbComment","limFlag","eisFlag1","eisFlag2","eisLo","eisHi","nisFlag1","nisFlag2","dupFlag","msFlag","msLo","msHi","IRflag","rtFlag","mbFlag","dibmbFlag","ibFlag","dibibFlag", "ccvFlag", "fls", "flsLo", "flsHi","calFlag","iscFlag")
  df<-df[,!names(df) %in% more.cols.to.remove]
  #change remaining SCIEX names to GLENDA names
  colnames(df)[colnames(df)=="Sample.Name"]<-"field_id"
  colnames(df)[colnames(df)=="Sample.Type"]<-"sampleType"
  colnames(df)[colnames(df)=="Component.Name"]<-"analyte"
  colnames(df)[colnames(df)=="sampleTime"]<-"time_utc"
  #create field_id

  df$field_id<-if_else(df$resultQcIdentifier %in% c("RFS", "LD","LMS"), if_else(df$sedType=="Core", paste0(sampleYr, "-", substr(df$field_id,1,4),":1C_", as.character(df$topOfInterval), "-", as.character(df$bottomOfInterval)),paste0(sampleYr, "-", substr(df$field_id,1,4),":SS")),df$field_id)
  df$field_id<-gsub("_0-1","_00.0-01.0", df$field_id)
  df$field_id<-gsub("_1-2","_01.0-02.0", df$field_id)
  df$field_id<-gsub("_2-3","_02.0-03.0", df$field_id)
  df$field_id<-gsub("_3-4","_03.0-04.0", df$field_id)
  df$field_id<-gsub("_4-5","_04.0-05.0", df$field_id)
  df$field_id<-gsub("_5-6","_05.0-06.0", df$field_id)
  df$field_id<-gsub("_6-7","_06.0-07.0", df$field_id)
  df$field_id<-gsub("_7-8","_07.0-08.0", df$field_id)
  df$field_id<-gsub("_8-9","_08.0-09.0", df$field_id)
  df$field_id<-gsub("_9-10","_09.0-10.0", df$field_id)
  df$field_id<-gsub("_10-11","_10.0-11.0", df$field_id)
  df$field_id<-gsub("_11-12","_11.0-12.0", df$field_id)
  df$field_id<-gsub("_12-13","_12.0-13.0", df$field_id)
  df$field_id<-gsub("_13-14","_13.0-14.0", df$field_id)
  df$field_id<-gsub("_14-15","_14.0-15.0", df$field_id)
  df$field_id<-gsub("_15-16","_15.0-16.0", df$field_id)
  df$field_id<-gsub("_16-17","_16.0-17.0", df$field_id)
  df$field_id<-gsub("_17-18","_17.0-18.0", df$field_id)
  df$field_id<-gsub("_18-19","_18.0-19.0", df$field_id)
  df$field_id<-gsub("_19-20","_19.0-20.0", df$field_id)
  df$field_id<-gsub("_20-21","_20.0-21.0", df$field_id)
  df$field_id<-gsub("_21-22","_21.0-22.0", df$field_id)
  df$field_id<-gsub("_22-23","_22.0-23.0", df$field_id)
  df$field_id<-gsub("_23-24","_23.0-24.0", df$field_id)
  df$field_id<-gsub("_24-25","_24.0-25.0", df$field_id)
  df$field_id<-gsub("_25-26","_25.0-26.0", df$field_id)
  df$field_id<-gsub("_26-27","_26.0-27.0", df$field_id)
  df$field_id<-gsub("_27-28","_27.0-28.0", df$field_id)
  df$field_id<-gsub("_28-29","_28.0-29.0", df$field_id)
  df$field_id<-gsub("_28-29","_28.0-29.0", df$field_id)
  df$field_id<-gsub("_29-30","_29.0-30.0", df$field_id)
  df$field_id<-gsub("_30-31","_30.0-31.0", df$field_id)
  df$field_id<-gsub("_31-32","_31.0-32.0", df$field_id)
  df$field_id<-gsub("_32-33","_32.0-33.0", df$field_id)
  df$field_id<-gsub("_33-34","_33.0-34.0", df$field_id)
  df$field_id<-gsub("_34-35","_34.0-35.0", df$field_id)
  df$field_id<-gsub("_35-36","_35.0-36.0", df$field_id)
  df$field_id<-gsub("_36-37","_36.0-37.0", df$field_id)
  df$field_id<-gsub("_37-38","_37.0-38.0", df$field_id)
  df$field_id<-gsub("_38-39","_38.0-39.0", df$field_id)
  df$field_id<-gsub("_39-40","_39.0-40.0", df$field_id)
  #half cm intervals
  df$field_id<-gsub("_0-0.5","_00.0-00.5", df$field_id)
  df$field_id<-gsub("_0.5-1","_00.5-01.0", df$field_id)
  df$field_id<-gsub("_1-1.5","_01.0-01.5", df$field_id)
  df$field_id<-gsub("_1.5-2","_01.5-02.0", df$field_id)
  df$field_id<-gsub("_2-2.5","_02.0-02.5", df$field_id)
  df$field_id<-gsub("_2.5-3","_02.5-03.0", df$field_id)
  df$field_id<-gsub("_3-3.5","_03.0-03.5", df$field_id)
  df$field_id<-gsub("_3.5-4","_03.5-04.0", df$field_id)
  df$field_id<-gsub("_4-4.5","_04.0-04.5", df$field_id)
  df$field_id<-gsub("_4.5-5","_04.5-05.0", df$field_id)
  df$field_id<-gsub("_5-5.5","_05.0-05.5", df$field_id)
  df$field_id<-gsub("_5.5-6","_05.5-06.0", df$field_id)
  df$field_id<-gsub("_6-6.5","_06.0-06.5", df$field_id)
  df$field_id<-gsub("_6.5-7","_06.5-07.0", df$field_id)
  df$field_id<-gsub("_7-7.5","_07.0-07.5", df$field_id)
  df$field_id<-gsub("_7.5-8","_07.5-08.0", df$field_id)
  df$field_id<-gsub("_8-8.5","_08.0-08.5", df$field_id)
  df$field_id<-gsub("_8.5-9","_08.5-09.0", df$field_id)
  df$field_id<-gsub("_9-9.5","_09.0-09.5", df$field_id)
  df$field_id<-gsub("_9.5-10","_09.5-10.0", df$field_id)
#add LD to duplicates
  df$field_id<-if_else(df$sampleType=="LD1", paste0(df$field_id, "LD1"), if_else(df$sampleType=="LD2", paste0(df$field_id,"LD2"), df$field_id), df$field_id)
#add site ID column
  df$site_id<-if_else(df$resultQcIdentifier %in% c("RFS", "LD", "LMS"), substr(df$field_id,1,9), NA)
#add repnumber
  df$repNumber<-if_else(df$sampleType=="LD2", 1, 0)
#add analytic method
  df$analyticMethod<-sop
  #add submissionId
  df$submissionId<-paste0(sampleYr, LkAbrr, "-T:PFAS_NRRI", as.character(format(Sys.Date(), "%Y%m%d")))
#add sordOrder
  df$analyteSortOrder<-sortOrder
#change field_id with the order
  df$field_id<-if_else(df$resultQcIdentifier %in% c("CLB","CLC","LVM","CAL","LMB","OPR"), paste0(df$field_id,"-",df$analyteSortOrder), df$field_id)
#add labsampleID
  df$labSampleId<-NA
  #add labBatchId
  df$labBatchId<-paste0(sampleYr, "-",df$fieldBlank, "_", df$samplePrepDate,"_",df$analyteSortOrder)
  df$labBatchId<-if_else(df$resultQcIdentifier %in% c("CLB","CLC","LVM","CAL","LMB","OPR"), NA, df$labBatchId)
#add unique_id
  df$unique_id<-df$field_id
  #add recovery limits
  df$recoveryUl<-if_else(df$analyteType=="SUR", 150, if_else(df$analyteType=="IS", NA, if_else(df$sampleType %in% c("OPR","LMS"), 150, if_else(df$sampleType %in% c("LVM","CLC","CAL"), 130, NA))))
  df$recoveryLl<-if_else(df$analyteType=="SUR", 20, if_else(df$analyteType=="IS", 30, if_else(df$sampleType %in% c("OPR","LMS"), 40, if_else(df$sampleType %in% c("LVM","CLC","CAL"), 70, NA))))
 #add result stat type
  df$resultStatType<-"ACT"
   #add sig figs
  df$resultSignFigures<-if_else(df$result>df$quantificationLimit, 3, 2)
  df$resultSignFigures<-if_else(df$analyteType=="TRG", df$resultSignFigures, NA)
  #column order
  col.order<-c("site_id","field_id","submissionId","topOfInterval","bottomOfInterval","depthUnit","sampleType","sampleDate","time_utc","resultQcIdentifier","repNumber","analyte","analyteCasNum","analyteType","analyteSortOrder","analyticMethod","result","resultUnits","resultSignFigures","resultStatType", "mdl","quantificationLimit", "recovery","recoveryLl","recoveryUl","spikeLevel","spikeUnits","dilutionFactor","samplePrepDate","analysisDateTime","qualifier","weightVolumeAnalyzed","weightVolumeUnits","weightBasis","resultComments","lat","long","depth_m","labBatchId","labSampleId", "unique_id")
  #rearrange
  df<-select(df, col.order)

  return(df)
  }
