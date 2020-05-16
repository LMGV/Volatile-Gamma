# basic functions

## plotting Functiions ----

value_at_risk_empirical = function(data, significance_level, time) {
  # data = returns in correct time scale
  # signif level= bound  for var, signifince
  # time = scaling factor for VaR with different timeframe. Requires normality
  if(missing(time)) { time=1}
  value_at_risk = -quantile(data,significance_level, na.rm=TRUE)*sqrt(time)
  return(value_at_risk)
}

expected_shortfall_empirical = function(data, significance_level, time) {
  # data = returns in correct time scale
  # signif level= bound  for var, signifince
  # time = scaling factor for VaR with different timeframe. Requires normality
  value_at_risk = -quantile(data,significance_level, na.rm=TRUE)
  expected_shortfall = -mean(subset(data, data < -value_at_risk), na.rm=TRUE)*sqrt(time)
  return(expected_shortfall)
}

## plotting function ----
line_plot_multiple = function(title, x,xlab, ylab, names_y, y_percent, y_discrete, legend, y1, y2,y3,y4,y5,y6,y7,y8,y9,y10) {
  # plots up to 10 lines with same scale 
  # insert vectors for x, y1, y2... y2-10 are optional
  # give vector for legend entry names_y, ...
  # if names_y na or missing colors are still given
  # legend: boolean if legend should be created, Default is TRUE
  if(missing(xlab)) {xlab=NULL}
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  if(missing(names_y)==T | all(is.na(names_y))==T) {names_y=c(1:10)}
  if(missing(legend)){legend = T}
  if(missing(y_percent)){y_percent = FALSE}
  
  plot = ggplot(data=NULL, aes(x=x))+
    geom_line(aes(y=y1, col = names_y[1])) 
  
  # only add layers if value provided
  if (missing(y2)==F) {
    plot = plot + geom_line(aes(y =y2 , col =names_y[2]))
  }  
  if (missing(y3)==F) {
    plot = plot + geom_line(aes(y =y3 , col =names_y[3]))
  }  
  if (missing(y4)==F) {
    plot = plot + geom_line(aes(y =y4 , col =names_y[4]))
  }  
  if (missing(y5)==F) {
    plot = plot + geom_line(aes(y =y5 , col =names_y[5]))
  }  
  if (missing(y6)==F) {
    plot = plot + geom_line(aes(y =y6 , col =names_y[6]))
  }  
  if (missing(y7)==F) {
    plot = plot + geom_line(aes(y =y7 , col =names_y[7]))
  }  
  if (missing(y8)==F) {
    plot = plot + geom_line(aes(y =y8 , col =names_y[8]))
  }  
  if (missing(y9)==F) {
    plot = plot + geom_line(aes(y =y9 , col =names_y[9]))
  }  
  if (missing(y10)==F) {
    plot = plot + geom_line(aes(y =y10 , col =names_y[10]))
  }  
  
  
  plot = plot + theme_bw()+
    labs(x = xlab)+
    labs(y = ylab)
  
  # add legend if Legend is true are given
  if(legend==T) {
    plot = plot + theme(legend.title = element_blank(), legend.position = "bottom", legend.box.background = element_rect(colour = "black"))
  } else {
    plot = plot + theme(legend.position = "none")
  }
  
  # add units to yaxis (e.g. percent)
  if(y_percent==T) {
    plot = plot + scale_y_continuous(labels = scales::percent_format(accuracy = 2))
  }
  
  if(y_discrete==T) {
    plot = plot + scale_y_discrete()
  }
  
  # other plot options
  plot = plot +  ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  
  print(plot)
  return(plot)
}

line_point_plot_multiple = function(title, x,xlab, ylab, names_y, y_percent,  y_discrete, legend, y1, y2,y3,y4,y5,y6,y7,y8,y9,y10) {
  # plots up to 10 lines with same scale 
  # insert vectors for x, y1, y2... y2-10 are optional
  # give vector for legend entry names_y, ...
  # if names_y na or missing colors are still given
  # legend: boolean if legend should be created, Default is TRUE
  # y6-y10 are geom point
  if(missing(xlab)) {xlab=NULL}
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  if(missing(names_y)==T | all(is.na(names_y))==T) {names_y=c(1:10)}
  if(missing(legend)){legend = T}
  if(missing(y_percent)){y_percent = FALSE}
  
  #first geom point
  plot = ggplot(data=NULL, aes(x=x))+
    geom_point(aes(y=y1, col = names_y[1])) 
  
  # geom lines
  if (missing(y2)==F) {
    plot = plot + geom_line(aes(y =y2 , col =names_y[2]))
  }  
  if (missing(y3)==F) {
    plot = plot + geom_line(aes(y =y3 , col =names_y[3]))
  }  
  if (missing(y4)==F) {
    plot = plot + geom_line(aes(y =y4 , col =names_y[4]))
  }  
  if (missing(y5)==F) {
    plot = plot + geom_line(aes(y =y5 , col =names_y[5]))
  }  
  
  # geom points
  if (missing(y6)==F) {
    plot = plot + geom_point(aes(y =y6 , col =names_y[6]))
  }  
  if (missing(y7)==F) {
    plot = plot + geom_point(aes(y =y7 , col =names_y[7]))
  }  
  if (missing(y8)==F) {
    plot = plot + geom_point(aes(y =y8 , col =names_y[8]))
  }  
  if (missing(y9)==F) {
    plot = plot + geom_point(aes(y =y9 , col =names_y[9]))
  }  
  if (missing(y10)==F) {
    plot = plot + geom_point(aes(y =y10 , col =names_y[10]))
  }  
  
  
  plot = plot + theme_bw()+
    labs(x = xlab)+
    labs(y = ylab)
  
  # add legend if Legend is true are given
  if(legend==T) {
    plot = plot + theme(legend.title = element_blank(), legend.position = "bottom", legend.box.background = element_rect(colour = "black"))
  } else {
    plot = plot + theme(legend.position = "none")
  }
  
  # add units to yaxis (e.g. percent)
  if(y_percent==T) {
    plot = plot + scale_y_continuous(labels = scales::percent_format(accuracy = 2))
  }
  
  # y discrete=
  if(y_discrete==T) {
    plot = plot + scale_y_discrete()
  }
  
  # other plot options
  plot = plot +  ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+ 
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  
  print(plot)
  return(plot)
}


# times series plot function
time_series_plot_basic = function(date,title, y, ylab, min_date, max_date) {
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  
  plot = ggplot(data=NULL, aes(x=as.Date(date,format = "%m/%d/%Y"), y = y))+
    geom_line(colour = "cornflowerblue", size = 1)+
    scale_x_date(limits = c(min=min_date, max=max_date)) + 
    theme_bw()+
    labs(x = "Date")+
    labs(y = ylab)+
    theme(legend.position = "bottom", legend.box.background = element_rect(colour = "black"))+
    ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  plot
  return(plot)
}



time_series_plot_multiple = function(date,title, y, ylab, seperating_var,min_date, max_date) {
  if(missing(ylab)) {ylab=NULL}
  if(missing(title)) {title=NULL}
  
  plot = ggplot(data=NULL, aes(x=as.Date(date,format = "%m/%d/%Y"), y = y))+
    geom_line(aes(color = seperating_var), size = 1)+
    scale_x_date(limits = c(min=min_date, max=max_date)) + 
    theme_bw()+
    labs(x = "Date")+
    labs(y = ylab)+
    theme(legend.position = "bottom", legend.box.background = element_rect(colour = "black"))+
    ggtitle(paste(title,sep=" ")) +
    theme(plot.title = element_text(size=10, face="bold"))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    ggsave(file=paste0("output/",title,".png"), width=6, height=4, dpi=600)
  
  print(plot)
  return(plot)
}
