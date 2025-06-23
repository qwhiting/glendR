# this function calculated the EIS recovery for each sample (excluding double blanks, solvents, and standard samplt types)


#' @import tidyr
#' @import dplyr
#' @export


EISrec<-function(df, EIS, NIS, matchDF){
  #filter for all samples that have EIS and NIS such that we can calculate the EIS recovery later
  rf.all<-subset(df, Sample.Type!="Solvent" & Sample.Type != "Double Blank")
  rf.all$Component.Name<-gsub(" ","", rf.all$Component.Name)
  rf.eis<-rf.all%>%filter(Component.Name %in% EIS)
  rf.nis<-rf.all%>%filter(Component.Name %in% NIS)
  #match NIS to corresponding EIS
  rf.eis$IS.Name<-matchDF$NIS[match(rf.eis$Component.Name, matchDF$EIS)]
  #create list of cals to add in the NIS areas to the corresponding EIS
  nis.list<- split(rf.nis, rf.nis$Sample.Name)
  eis.list<-split(rf.eis, rf.eis$Sample.Name)
  #loop to add in the areas for each sample
  for(i in 1:length(eis.list)){
    eis.list[[i]]$IS.Area <- nis.list[[i]]$Area[match(eis.list[[i]]$IS.Name, nis.list[[i]]$Component.Name)]
    eis.list[[i]]$IS.Actual.Concentration<-nis.list[[i]]$Actual.Concentration[match(eis.list[[i]]$IS.Name, nis.list[[i]]$Component.Name)]
  }

  rf<-do.call("rbind", eis.list)
  #calculate RF for each sample and EIS/NIS combo
  c.rf<-rf%>%filter(Sample.Type=="Standard")%>%filter(Used=="True")
  c.rf$rf<- (c.rf$Area * c.rf$IS.Actual.Concentration)/(c.rf$IS.Area * c.rf$Actual.Concentration)
  #average RF in calset
  RF <- c.rf%>%group_by(Component.Name)%>%summarise(RF=mean(rf, na.rm=T), stdev=sd(rf, na.rm=T))
  #EIS recovery for non standards (unknown, quality control, blanks)
  EISrec<- subset(rf, Sample.Type!="Standard")
  #get the average RF for each EIS/NIS combo
  EISrec$rf<-RF$RF[match(EISrec$Component.Name, RF$Component.Name)]
  #calculate EIS concentration based off of area and RF
  EISrec$Calculated.Concentration <- (EISrec$Area * EISrec$IS.Actual.Concentration)/
    (EISrec$IS.Area * EISrec$rf)
  #EIS recovery (%)
  EISrec$sr<- round((EISrec$Calculated.Concentration/EISrec$Actual.Concentration) *100, digits = 1)
  #create unique temp ID to match back to main df with
  EISrec$matchId<-paste(EISrec$Sample.Name, EISrec$Component.Name)
  #filter to only have sample name, component name, surrogatee recovery, and matching id
  EISrec2<-EISrec[,c("Sample.Name", "Component.Name", "Calculated.Concentration","sr","matchId")]

  return(EISrec2)
}
