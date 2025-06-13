#this function calculates the NIS recovery for each sample


#' @import dplyr
#' @import tidyr

#' @export


NISrec<-function(df, NIS){
  #filter for NIS and get average are in calset
  c.nis<-df%>%filter(Component.Name %in% NIS)%>%filter(Sample.Type=="Standard")%>%filter(Used=="True")
  c.nis<-c.nis%>%group_by(Component.Name)%>%summarise(area=mean(Area))
  #add the recovery to the samples in the dataframe
  rec<-df
  rec$nisRec<-if_else(rec$Sample.Type %in% c("Quality Control", "Unknown", "Blank"),
                      round((rec$Area/c.nis$area[match(rec$Component.Name, c.nis$Component.Name)])*100, digits=1),
                      NA)
  rec<-rec%>%filter(Component.Name %in% NIS)
  return(rec)
}




