#this calculates the deviation in retention time for each analyte as it compares to the retention times in the calset
#might have to do something about NA values (make 0 then if-else to get back to NA?)
#add flag if >0.25 min difference (15 sec)

#' @import dplyr
#' @import tidyr
#' @export

retention<-function(df){
  #get average RT for each PFAS in the calset
  cal<-df%>%filter(Sample.Type=="Standard")%>%filter(Used=="True")
  calRT<-cal%>%group_by(Component.Name)%>%summarise(RT=mean(Retention.Time, na.rm=T))

  rt.list<-split(df, df$Sample.Name, drop=F)
  for(i in 1:length(rt.list)){
    rt.list[[i]]$deltaRT<- round(abs(rt.list[[i]]$Retention.Time - calRT$RT[match(rt.list[[i]]$Component.Name, calRT$Component.Name)]), 2)
  }
  rt<-do.call("rbind", rt.list)

  rt$deltaRT<-if_else(rt$Sample.Type=="Standard", NA, rt$deltaRT)
  rt$matchId<-paste(rt$Sample.Name, rt$Component.Name)
  rt<-rt%>%filter(Sample.Type %in% c("Unknown","Blank","Quality Control"))
  rt2<-rt[,c("Sample.Name","Component.Name","Retention.Time","deltaRT","matchId")]
  return(rt2)

}
