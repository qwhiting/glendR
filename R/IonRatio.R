#this function calculates the ion ratio between the quantifier and qualifier for each analyte that had them
#only requred on samples that are 2x the MDL


#' @import dplyr
#' @import tidyr
#' @export


IonRatio<-function(df){
  quant<-df%>%filter(Component.Type=="Quantifiers")
  qual<-df%>%filter(Component.Type=="Qualifiers")
  #make lists by sample
  quant.list <- split(quant, quant$Sample.Name, drop = F)
  qual.list <- split(qual, qual$Sample.Name, drop =F)
  #loop through the lists
  for(i in 1:length(quant.list)){
    quant.list[[i]]$Ion.Ratio <- quant.list[[i]]$Area /
      qual.list[[i]]$Area[match(quant.list[[i]]$Component.Group.Name, qual.list[[i]]$Component.Group.Name)]
  }
  #bind list back together that now has the ion ratios for all analytes (drop the qalifiers bc no longer need them)
  ir<-do.call("rbind", quant.list)
  #get the IS so that we dont loose any data
  standards<-df%>%filter(Component.Type=="Internal Standards")
  standards$Ion.Ratio<-"NA"
  ir2<-rbind(ir,standards)
  #get ion ratios for the cal standards
  c.ir<-subset(ir2, Sample.Type=="Standard" & Ion.Ratio!= "NA")
  #c.ir$zscore<- abs(((c.ir$Ion.Ratio - mean(c.ir$Ion.Ratio)) / sd(c.ir$Ion.Ratio)))
  c.ir<-c.ir %>% group_by(Component.Name) %>% summarise(IR=mean(as.numeric(Ion.Ratio), na.rm=T))
  #check the ion ratios (>30% of what is in the cal set)
  ir2$Standard.Ion.Ratio<-if_else(ir2$Sample.Type %in% c("Quality Control", "Unknown", "Blank"),c.ir$IR[match(ir2$Component.Name, c.ir$Component.Name)],NA)

  ir2$Ion.Ratio<-if_else(ir2$Ion.Ratio=="NA", NA, ir2$Ion.Ratio)
  ir2$Ion.Ratio<-as.numeric(ir2$Ion.Ratio)
  ir2$Ion.Ratio.Check<-round((ir2$Ion.Ratio/ir2$Standard.Ion.Ratio)*100, digits = 1)

  return(ir2)
}








