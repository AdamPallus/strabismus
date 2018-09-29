#plan: 
#1 - choose appropriate kernerls for convolving to find velocity and spike density functions
#2 - mark slow vergence movements
#3 - Train the model using just those movements 
#4 - Use the model to predict for the whole data set

#Functions----
modelSlowVerg<- function(x,model.form='verg.velocity~sdf20+verg.angle',
                    lagsdf=31,parabolic_n=15,sdf_sd=25,
                    returnmodel=FALSE){
  source('markSlowVerg.R')
  
  if (!('time' %in% names(x))){
    x<- mutate(x,time=row_number())
  }
  
  x %>%
    mutate(verg.velocity=parabolicdiff(lep-rep,parabolic_n),
           rev=parabolicdiff(rep,parabolic_n),
           lev=parabolicdiff(lep,parabolic_n),
           revV=parabolicdiff(repV,parabolic_n),
           levV=parabolicdiff(lepV,parabolic_n),
           sdf=spikedensity(rasters,sd=sdf_sd),
           sdf20=dplyr::lag(sdf,lagsdf),
           conj.velocity=sqrt(((rev+lev)/2)^2+((revV+levV)/2)^2),
           slowverg=markSlowVerg(verg.velocity,
                                 highthresh = 12,
                                 lowthresh=4,
                                 rejectthresh=40,
                                 min.dur=150)) ->
    x
  
  x<-  mutate(x,saccadic=!is.na(markSaccadesDouble(conj.velocity,markFixations = FALSE)))
  # x<-  mutate(x,saccadic=!is.na(markSaccadesDoubleTEST(conj.velocity,
  #                                                      threshold1=50,
  #                                                      threshold2=20,
  #                                                      driftcorrect=FALSE,
  #                                                      markFixations=FALSE)))
  
  mod<- lm(model.form,data=filter(x,!is.na(slowverg),!saccadic))
  
  # mod<- lm(model.form,data=filter(x,!is.na(slowverg)))
  
  
  # mod<- lm(model.form,data=filter(x,saccadic | !is.na(slowverg)))
  
  if (returnmodel) return(mod)
  print(tidy(mod))
 
  message(paste('R-squared =', summary(mod)$r.squared))
  # x<- dplyr::select(x,-sacnum,-counter)
  x<- mutate(ungroup(x),predV=predict(mod,newdata=x),
             showrasters=replace(rasters,rasters<1,NA))
  
}


MakeFigure<-function(d,labelloc=288000,zoom=TRUE){
  require(gridExtra)
  # d<-mutate(d,showsaccadic=replace(saccadic,!saccadic,NA))
  
  Hpositions<- ggplot(d)+
    geom_line(aes(time,(rep)),color='red',size=1)+
    geom_line(aes(time,(lep)),color='blue',size=1)+
    geom_line(aes(time,thp),color='black',linetype=2)+
    theme_minimal()+xlab('')+theme(axis.text.x = element_blank())+
    ylab('Horizontal\nEye Positon (deg)')
  
  Vpositions<- ggplot(d)+
    geom_line(aes(time,(((repV+lepV)/2))),color='violet',size=1)+
    geom_line(aes(time,tvp),color='black',linetype=2)+
    theme_minimal()+xlab('')+theme(axis.text.x = element_blank())+
    ylab('Vertical\nEye Position (deg)')
  
  VergPositions<- ggplot(d)+
    geom_line(aes(time,lep-rep),color='darkgreen',size=1)+
    theme_minimal()+xlab('')+theme(axis.text.x = element_blank())+
    ylab('Vergence Position\n(deg)')
  
  VergVelocities<-ggplot(d)+theme_minimal()+
    geom_line(aes(time,verg.velocity),color='blue',alpha=1)+
    geom_line(aes(time,predV),color='orange')+
    # geom_line(aes(time,predV2),color='magenta')+
    # geom_point(aes(time,as.numeric(showsaccadic)),alpha=0.2)+
    # ylim(c(-10,30))+
    xlab('')+theme(axis.text.x = element_blank())+
    ylab('Vergence Velocity\n(deg/s)')
  # annotate('text',x=labelloc,y=-5,label='Predicted',color='orange')
  
  if (zoom){
    # VergVelocities<-VergVelocities+ylim(c(-10,30))
    VergVelocities<-VergVelocities+ylim(c(-30,30))
  }
  
  
  sdf<- ggplot(d)+theme_minimal()+
    geom_area(aes(time,sdf),alpha=1)+
    geom_point(aes(time,showrasters+5),shape='|',color='lightgrey',size=3)+
    ylab('Firing Rate\n(spk/s)')
  xlab('Time (ms)')
  
  # multiplot(Hpositions,Vpositions,VergPositions,VergVelocities,sdf)
  # p<- list(Hpositions,Vpositions,VergPositions,VergVelocities,sdf)
  grid.arrange(arrangeGrob(Hpositions,Vpositions,VergPositions,VergVelocities,sdf,ncol=1,
                           heights=c(.5,.3,.5,1,.5)),padding=0.2)
  
  # grid.arrange(arrangeGrob(Hpositions,Vpositions,VergPositions,VergVelocities,VergVelocities+ylim(c(-10,30)),sdf,ncol=1,
  #                          heights=c(.5,.3,.5,0.5,1,.5)),padding=0.2)
}
  
  
MakeFigureSDF<-function(d,labelloc=288000,zoom=TRUE){
    require(gridExtra)
    # d<-mutate(d,showsaccadic=replace(saccadic,!saccadic,NA))
    
    Hpositions<- ggplot(d)+
      geom_line(aes(time,(rep)),color='red',size=1)+
      geom_line(aes(time,(lep)),color='blue',size=1)+
      geom_line(aes(time,thp),color='black',linetype=2)+
      theme_minimal()+xlab('')+theme(axis.text.x = element_blank())+
      ylab('Horizontal\nEye Positon (deg)')
    
    Vpositions<- ggplot(d)+
      geom_line(aes(time,(((repV+lepV)/2))),color='violet',size=1)+
      geom_line(aes(time,tvp),color='black',linetype=2)+
      theme_minimal()+xlab('')+theme(axis.text.x = element_blank())+
      ylab('Vertical\nEye Position (deg)')
    
    VergPositions<- ggplot(d)+
      geom_line(aes(time,lep-rep),color='darkgreen',size=1)+
      theme_minimal()+xlab('')+theme(axis.text.x = element_blank())+
      ylab('Vergence Position\n(deg)')
    
    VergVelocities<-ggplot(d)+theme_minimal()+
      geom_line(aes(time,verg.velocity),color='blue',alpha=1)+
      geom_hline(yintercept = 0)+
      # geom_line(aes(time,predV2),color='magenta')+
      # geom_point(aes(time,as.numeric(showsaccadic)),alpha=0.2)+
      # ylim(c(-10,30))+
      xlab('')+theme(axis.text.x = element_blank())+
      ylab('Vergence Velocity\n(deg/s)')
    # annotate('text',x=labelloc,y=-5,label='Predicted',color='orange')
    
    if (zoom){
      # VergVelocities<-VergVelocities+ylim(c(-10,30))
      VergVelocities<-VergVelocities+ylim(c(-30,30))
    }
  
  
  sdf<- ggplot(d)+theme_minimal()+
    geom_area(aes(time,sdf),alpha=1)+
    geom_line(aes(time,predV),color='orange')+
    geom_point(aes(time,showrasters+5),shape='|',color='lightgrey',size=3)+
    ylab('Firing Rate\n(spk/s)')
  xlab('Time (ms)')
  
  # multiplot(Hpositions,Vpositions,VergPositions,VergVelocities,sdf)
  # p<- list(Hpositions,Vpositions,VergPositions,VergVelocities,sdf)
  grid.arrange(arrangeGrob(Hpositions,Vpositions,VergPositions,VergVelocities,sdf,ncol=1,
                           heights=c(.5,.3,.5,.5,1)),padding=0.2)
  
  # grid.arrange(arrangeGrob(Hpositions,Vpositions,VergPositions,VergVelocities,VergVelocities+ylim(c(-10,30)),sdf,ncol=1,
  #                          heights=c(.5,.3,.5,0.5,1,.5)),padding=0.2)
}

#Running----
t<- loadnewcsv2(path="C:/Users/setup/Desktop/NRTP Vergence/SOASTRAB/patos/")

t %>% group_by(neuron) %>%
  mutate(conj.vertical=(repV+lepV)/2)->
  t

z<- modelSlowVerg(filter(t,neuron=='Patos-103'),
                  # model.form='verg.velocity~sdf20+verg.angle+conj.vertical',
                  model.form='sdf20 ~ verg.velocity+verg.angle+conj.vertical',
                  # model.form='sdf20 ~ verg.velocity+verg.angle',
                  # model.form='verg.velocity~sdf20+verg.angle',
                         lagsdf=20,parabolic_n=10,sdf_sd=10)

window_size=5000
manipulate({
  d=filter(z,time>=window,time<window+window_size)
  d<- group_by(d,time) %>% summarize_all(funs(first))
  MakeFigureSDF(d,zoom=FALSE)
  # MakeFigure(d,zoom=TRUE)
},
window=slider(window_size,max(z$time-window_size),step=window_size))


movements=unique(z$slowverg)

manipulate({
  d=filter(z,slowverg==movements[chosenMove])
  MakeFigureSDF(d,zoom=TRUE)
},
chosenMove=slider(2,length(movements))
)

#One of the free parameters is the maximum velocity to accept as a slow vergence movement
MakeFigure(filter(z,time>25000,time<35000),zoom=FALSE)

MakeFigureSDF(filter(z,time>25000,time<35000),zoom=FALSE)


MakeFigureSDF(filter(z,time>23500,time<26000),zoom=FALSE)

MakeFigureSDF(filter(z,time>37000,time<47000),zoom=FALSE)

MakeFigureSDF(filter(z,time>39000,time<45000),zoom=FALSE)



MakeFigureSDF(filter(z,time>660500,time<678000),zoom = FALSE)


z<- modelSlowVerg(filter(tt,neuron=='Bee-33'),
                  model.form='verg.velocity~sdf20+verg.angle',
                  
                  lagsdf=30,parabolic_n=10,sdf_sd=15)


MakeFigure(filter(z,time>110000,time<112000),zoom = FALSE)

tt %>% 
  # filter(neuron == 'Bee-15') %>%
  group_by(neuron) %>%
  summarize(fit=as.numeric(glance(modelSlowVerg(.,model.form='verg.velocity~sdf20+verg.angle',
                              lagsdf=30,parabolic_n=20,sdf_sd=25,returnmodel = TRUE))$r.squared))->
  tfit

arrange(tfit,desc(fit))
