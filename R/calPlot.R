#this is a function to plot the calibration curve and to see the R2 for each PFAS

#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import ggpmisc
#' @export

#the function
calPlot<-function(df, Title, Subtitle){
  cal<-df%>%filter(Sample.Type=="Standard")%>%filter(Used=="True")%>%filter(Component.Type=="Quantifiers")

  p<-ggplot(cal, aes(Concentration.Ratio, Area.Ratio, weight=1/Concentration.Ratio))+
    stat_poly_line(method = "lm", col="coral2")+
    geom_point(col="dodgerblue")+
    stat_poly_eq(use_label(c("eq")), size=3, label.y = 0.95)+
    stat_poly_eq(use_label(c("R2")), size=4,  rr.digits = 4, label.y =0.75 )+
    facet_wrap(~Component.Name, scales = "free")+
    labs(x="Concentration Ratio", y="Area Ratio", title=Title,
         subtitle = Subtitle)+
    theme(plot.title = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  p
}
