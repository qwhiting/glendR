#changes to analyte type from Component.Type

#' @import tidyr
#' @import dplyr
#' @export

analyteType<-function(df, NIS){
  df$Component.Type<-if_else(df$Component.Type=="Quantifiers", "TRG", if_else(df$Component.Name %in% NIS, "IS", "SUR"))
  colnames(df)[which(names(df) == "Component.Type")] <- "analyteType"
  return(df)
}
