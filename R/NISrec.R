#this function calculates the NIS recovery for each sample


#' @import dplyr
#' @import tidyr

#' @export


NISrec<-function(df, NIS){
  #filter for NIS and get average are in calset
  c.nis<-df%>%filter(Component.Name %in% NIS)%>%filter(Sample.Type=="Standard")%>%filter(Used=="True")
  c.nis<-c.nis%>%group_by(Component.Name)%>%summarise(area=mean(Area, na.rm=T))
  #add the recovery to the samples in the dataframe
  rec<-df%>%filter(Sample.Type!="Standard")%>%filter(Sample.Type!="Double Blank")%>%filter(Sample.Type!="Solvent")
  rec$nisRec<-if_else(rec$Sample.Type %in% c("Quality Control", "Unknown", "Blank"),
                      round((rec$Area/c.nis$area[match(rec$Component.Name, c.nis$Component.Name)])*100, digits=1),
                      NA)
  rec<-rec%>%filter(Component.Name %in% NIS)
  rec$matchId<-paste(rec$Sample.Name, rec$Component.Name)
  rec2<-rec[,c("Sample.Name","Component.Name","nisRec","matchId")]
  return(rec2)
}




