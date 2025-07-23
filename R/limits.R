#this function is to flag any limits (UND/LTL/MDL/BQL/GTL)

#' @import tidyr
#' @import dplyr
#' @export

limits<-function(df, cal){
  df$undFlag<-if_else(df$analyteType=="TRG", if_else(df$result=0.0, "UND", NA), NA)
  cal$Component.Name<-gsub(" ","", cal$Component.Name)
  df$ltlFlag<-if_else(df$analyteType=="TRG", if_else(df$result*df$weightVolumeAnalyzed<cal$LLOD[match(df$Component.Name, cal$Component.Name)], "LTL", NA),NA)
  df$gtlFlag<-if_else(df$analyteType=="TRG", if_else(df$result*df$weightVolumeAnalyzed>cal$ULOD[match(df$Component.Name, cal$Component.Name)], "GTL", NA),NA)
  df$bqlFlag<-if_else(df$analyteType=="TRG", if_else(df$result>mdl, if_else(df$result<=df$quantificationLimit, "BQL", NA), NA),NA)
  df$mdlFlag<-if_else(df$analyteType=="TRG", if_else(df$result<= df$mdl, "MDL", NA), NA)
  df$mdlFlag<-if_else(df$undFlag=="UND", NA, df$mdlFlag)
  df$mdlFlag<-if_else(df$ltlFlag=="LTL", NA, df$mdlFlag)
  #combine limit flags into one column
  df$limFlag<-unite(df$undFlag, df$ltlFlag, df$mdlFlag, df$bqlFlag, df$gtlFlag, sep = ", ", na.rm = T)
  return(df)
  }
