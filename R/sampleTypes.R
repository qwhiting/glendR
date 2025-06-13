#this function creates the OCIDs for each sample in the data frame

#' @import tidyr
#' @import dplyr
#' @export


sampleTypes<-function(df, MS.Samples, LD1.Samples, LD2.Samples, Travel.Blanks, Field.Blanks){
  df$Sample.Type<-if_else(df$Sample.Name %in% MS.Samples, "LMS", NA)
  df$Sample.Type<-if_else(df$Sample.Name %in% LD1.Samples, "LD1", df$Sample.Type)
  df$Sample.Type<-if_else(df$Sample.Name %in% LD2.Samples, "LD2", df$Sample.Type)
  #df$Sample.Type<-if_else(df$Sample.Type=="LD", if_else(str_ends(df$Sample.Name, "_2"), "LD2", "LD1"), df$Sample.Type)
  df$Sample.Type<-if_else(str_detect(df$Sample.Name, "CCV"), "CLC", df$Sample.Type)
  df$Sample.Type<-if_else(str_detect(df$Sample.Name, "ISC"), "LVM", df$Sample.Type)
  df$Sample.Type<-if_else(str_detect(df$Sample.Name, "(?i)cal"), "CAL", df$Sample.Type)
  df$Sample.Type<-if_else(str_detect(df$Sample.Name, "IB"), "CLB", df$Sample.Type)
  df$Sample.Type<-if_else(str_detect(df$Sample.Name, "MB"), "LMB", df$Sample.Type)
  df$Sample.Type<-if_else(str_detect(df$Sample.Name, "OPR"), "OPR", df$Sample.Type)
  df$Sample.Type<-if_else(df$Sample.Name %in% Travel.Blanks, "FTB", df$Sample.Type)
  df$Sample.Type<-if_else(df$Sample.Name %in% Field.Blanks, "FRB", df$Sample.Type)
  df$Sample.Type<-if_else(df$Sample.Type %in% c("FRB","FTB","CAL","OPR","LMB","CLB","LVM","CLC","LD1","LD2","LMS"), df$Sample.Type, "RFS")

  df$resultQcIdentifier<-df$Sample.Type
  df$resultQcIdentifier<-if_else(df$Sample.Type=="LD1", "RFS", df$resultQcIdentifier)
  df$resultQcIdentifier<-if_else(df$Sample.Type=="LD2", "LD", df$resultQcIdentifier)

  return(df)
}
