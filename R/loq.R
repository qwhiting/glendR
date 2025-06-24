#this function is to calculate the LOQ from the MDL which should have already been imported to the df

#' @import tidyr
#' @import dplyr
#' @export



loq<-function(df){
  df$LOQ<-signif(df$mdl*3.18, digits = 4)
  df$quantificationLimit<-if_else(df$LOQ>7.5, 10, if_else(df$LOQ>3.5, 5, if_else(df$LOQ>1.5, 2, if_else(df$LOQ>0.75, 1, if_else(df$LOQ>0.35,0.5,if_else(df$LOQ>0.15, 0.2, if_else(df$LOQ>.075,0.1, if_else(df$LOQ>0.035,0.05, if_else(df$LOQ>0.015, 0.02, if_else(df$LOQ>0.0075, 0.01, if_else(df$LOQ>0.0035, 0.005, if_else(df$LOQ>0.0015, 0.002,0))))))))))))
  df<-df[,-which(names(df) %in% c("LOQ"))]
  return(df)
}
