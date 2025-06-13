# calculated the relative percent difference between the two duplicate samples

#' @import tidyr
#' @import dplyr
#' @export

rpd<-function(df){
  #remove the "_2" from the end of the SampleName
  df$Sample.Name<-if_else(df$Sample.Type=="LD2", substr(df$Sample.Name, 1, nchar(df$Sample.Name)-2), df$Sample.Name)
  df$Sample.Name<-if_else(df$Sample.Type=="LD2", paste(df$Sample.Name, "LD2", sep = ""), df$Sample.Name)
  #subset into LD1 and LD2 and create a 'mtching ID' to use later
  ld1<-df%>%filter(Sample.Type=="LD1")
  ld1$tempID<-paste(ld1$Sample.Name, ld1$Component.Name)
  ld1<-ld1%>%filter(analyteType=="TRG")
  ld2<-df%>%filter(Sample.Type=='LD2')
  ld2$tempID<-paste(ld2$Sample.Name, ld2$Component.Name)
  ld2$tempID<-gsub("LD2","", ld2$tempID)
  ld2<-ld2%>%filter(analyteType=="TRG")
  #get conc <MDL to be 0.5*MDL
  ld1$concMDL<-if_else(ld1$result>ld1$mdl, ld1$result, 0.5*ld1$mdl)
  ld2$concMDL<-if_else(ld2$result>ld2$mdl, ld2$result, 0.5*ld2$mdl)
  #get a df with ld1 and ld2 results
  rpd<-ld1[,c("tempID","concMDL")]
  rpd$concMDL2<-ld2$concMDL[match(rpd$tempID, ld2$tempID)]
  rpd$RPD<-abs(rpd$concMDL-rpd$concMDL2)/((rpd$concMDL+rpd$concMDL2)/2)*100
  #combine with original df
  df$tempID<-paste(df$Sample.Name, df$Component.Name)
  df$tempID<-gsub("LD2", "", df$tempID)
  df$rpd<-rpd$RPD[match(df$tempID, rpd$tempID)]
  df$rpd<-if_else(df$Sample.Type %in% c("LD1","LD2"), df$rpd, NA)
  df$rpd<-round(df$rpd, digits = 1)
  #add comments
  df$rpdCom<-paste0("RPD= ",as.character(df$rpd)," %")
  df$rpdCom<-gsub("RPD= NA %", "", df$rpdCom)
  #remove columns
  df<-df[ , -which(names(df) %in% c("tempID"))]
  return(df)
}
