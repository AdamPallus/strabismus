
#'I was working on this at the end of March 2018 because I noticed that my drift correct
#'algorithm was cutting of both ends of the saccades rather than just the ending. It also seemed
#'to me like it was cutting off too much of the ends of saccades and I had a whole discussion
#'with Mark about what is appropriate and I decided to stick with his algorithm since it's published
#'and only fails on tiny saccades that shouldn't really be analyzed anyway.
#'
#'I ended up includeing the function MarkSaccadesDoubleDRIFT in Adamhelperfunctions.R as the new
#'default MarkSaccadesDouble function. The previosu version is marked OLD

xx<- read.csv('SaccadeBootstrapCIs.csv')
neurons<- unique(xx$neuron)


ptb %>%
  filter(saccade.type %in% c('left.only','right.only'))%>%
  group_by(neuron,saccade.type) %>%
  do(r2=getR2(.))->
  ptb.r2

ptb %>%
  filter(saccade.type %in% c('left.only','right.only')) %>%
  group_by(neuron,saccade.type) %>%
  summarize(n=length(unique(sacnum))) ->
  ptb.n

fits<-left_join(ptb.r2,ptb.n)

plotr2<- spread(ptb.r2,saccade.type,r2)

qplot(right.only,left.only,data=plotr2)


plotr2<- separate(plotr2,neuron,c('monkey','cellnum'),remove = FALSE)


ptb.r2 %>%
  spread(saccade.type,r2) %>%
  separate(neuron,c('monkey','cellnum'),remove=FALSE) %>%
  mutate(eso=monkey %in% c('Kopachuck','QT'))->
  plotr2

qplot(right.only,left.only,data=filter(plotr2,neuron %in% goodcells$neuron),color=eso)+
  geom_abline()+
  theme_minimal()+
  coord_fixed()+
  labs(color='Esotrope')+
  xlab('Right Eye Viewing')+
  ylab('Left Eye Viewing')+
  ggtitle('Goodness of fit during Saccades')


arrange(plotr2,desc(left.only))


kk<-filter(ptb,neuron=='Kopachuck-117',saccade.type=='left.only')

ggplot(kk)+
  geom_line(aes(time,verg.angle))+
  geom_line(aes(time,sdf20),color='hotpink')

lm(kk)


#Show model fit
qt<- filter(ptb,neuron=='QT-107')
mod<- lm(sdf20~verg.angle+verg.velocity+conj.vertical,data=qt)
summary(mod)
qt$pred<-predict(mod,newdata=qt)


qt %>%
  filter(sacnum==22) %>%
  ggplot()+
  geom_line(aes(time,conj.velocity))

qt %>%
  group_by(sacnum) %>%
  summarize(startvel=conj.velocity[1])->
  qts

qt<- mutate(qt,sacnum2=markSaccadesDouble(conj.velocity,threshold1=50,threshold2=10,
                                          min.dur=20,driftcorrect = TRUE,markFixations = FALSE))

qt %>%
  group_by(sacnum2) %>%
  summarize(startvel=conj.velocity[1])->
  qts

qplot(startvel,data=qts)

qq<- filter(t,neuron=='QT-107')
qq<- filter(t,neuron=='Patos-103')

sn<-markSaccadesDoubleDRIFT(qq$conj.velocity,threshold1=40,threshold2=10,
                           min.dur=20,driftcorrect = TRUE,markFixations = FALSE)
qq$issaccade<-sn

manipulate(ggplot(filter(qq,time>startTime,time<startTime+5000))+
             geom_line(aes(time,conj.velocity))+
             geom_line(aes(time,conj.velocity*(issaccade*0+1)),color='hotpink',size=2,alpha=0.5)+
             ylim(0,1000),
           startTime=slider(1,nrow(qq),step=4000))




qq<- mutate(qq,sacnum2=markSaccadesDoubleTEST(conj.velocity,threshold1=40,threshold2=10,
                                          min.dur=20,driftcorrect = TRUE,markFixations = FALSE))
qq %>%
  group_by(sacnum2) %>%
  summarize(startvel=conj.velocity[1])->
  qts

qplot(startvel,data=qts)

markSaccadesDoubleDRIFT<- function(v, threshold1=60,threshold2=20,min.dur=5,maxreject=1000,
                              driftcorrect=FALSE,markFixations=TRUE){
  #This function is an R implementation of a two-threshold event marker
  #The algorithm works like this: Find all the times when velocity is above the high threshold
  #Extend this out until velocity is below the lower threshold
  
  #in practice, I'm identifying all the times that the saccades cross the low threshold and then rejecting 
  #any that don't meet the higher threshold
  #I'm also rejecting events that are below a certain duration
  
  #this algorithm also assigns positive ID numbers to the saccades and 
  #negative ID numbers to the non-saccades (fixations?)
  #after running this function, you can group_by(event) and measure the fixations or saccades as you wish
  require(dplyr)
  require(data.table) #for rbindlist - a fast version of do.call('rbind') that uses data.table
  datalength<-length(v)
  i<-which(abs(v)>threshold2) #find all the times when speed is above the lower threshold
  
  #For continuous numbers, diff=1. If there is a larger jump, it means there's a gap. 
  #That indicates another saccade
  sacoff<-which(diff(i)>1) #sacoff now contains the indices of the ends of all the saccades
  #sacon has all the indices of when the saccades start
  #After an offset, the next index will be the onset of the next saccade.
  #find the onsets by looking at the next index after the offset
  sacon<-c(1,sacoff+1) #first saccade always starts at first index
  sacoff<-c(sacoff,length(i)) #end of last saccade is always at the end
  event.onset<-i[sacon] #Convert from the indices to actual times
  event.offset<-i[sacoff] 
  
  #event.onset now has the time (in sample time) of all saccade onsets
  #set up stimes as a data.frame with two columns. Onset and offset. 
  stimes<- data.frame(event.onset,event.offset) 
  
  #this is a little function that works with the weirdness of R's "apply" family of functions
  #it just takes the onset and offset and returns the whole range. 
  #if you give it [10, 20], it will return [10 11 12 13 14 15 16 17 18 19 20]
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 4: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$event<- stimes[[4]]
    df$event.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes %>%
    mutate(dur=event.offset-event.onset, #calculate duration of each saccade
           s=row_number())-> #assign an ID number to each saccade
    stimes
  
  #Use "apply" to run the function "jsac" (above) on each line of the "stimes" data.frame
  #the result is the times of when all the saccades happen
  x<-rbindlist(apply(stimes,1,jsac)) 
  
  v<- data.frame(v=v) #Make the original velocity trace into a data.frame
  v<- mutate(v, time=row_number()) #add time to keep track
  
  #join the marked saccades and the velocity
  #the result is the velocity trace plus a row that just indicates whether there's a saccade
  #each saccade is identified by it's unique marker "event" that comes from df$event<- stimes[[4]] above 
  xx<- left_join(v,x,by='time')
  
  if (driftcorrect){
    
    xx%>%
      mutate(event2=event, #you're not allowed to change grouping variables, so make a copy first
             acc=parabolicdiff(v,7)) %>%
      group_by(event2) %>%
      mutate(realstart=first(time[v>50]),
             total.end=last(time),
             driftcut=first(time[time>(realstart+10)&abs(v)<100&abs(acc)<10000]),
             rejectdrift=time>driftcut,
             event=replace(event,time>driftcut,NA)) %>%
      ungroup()->
      xx
    
    # xx %>%
    #   ungroup() %>%
    #   mutate(acc=parabolicdiff(v,7)) %>%
    #   mutate(event=replace(event,abs(v)<100 & abs(acc)<10000,NA)) ->
    #   xx
  }
  
  xx %>%
    group_by(event) %>% #This means we analyze each saccade individually
    summarize(max.vel=max(abs(v)), #calculate max velocity
              dur=n()) %>% #calculate duration
    filter(max.vel>threshold1, #reject all saccades that fail to exceed the large threshold
           max.vel<maxreject, #reject all saccades over the max threshold
           dur>min.dur)-> #reject all saccades that fail to exceed the minimum duration
    xm #xm is a summary which means it just lists the saccades and their measured values
  
  xx %>% #go back to the full data set and now reject all the saccades that were rejected above
    filter(event %in% unique(xm$event)) %>% 
    dplyr::select(time,event) -> #All we need is the time and the eventID
    g
  
  if (markFixations){
    #this next part goes through and assigns an ID to all the non-saccade portions of the data
    stimes2<- filter(stimes,s %in% unique(xm$event))
    
    ftimes<-data.frame(fix.onset=c(1,stimes2$event.offset+1),
                       fix.offset=c(stimes2$event.onset-1,datalength))
    
    ftimes %>%
      filter(fix.onset>0,fix.offset>0)%>%
      mutate(dur=fix.offset-fix.onset,
             s=row_number()) %>%
      filter(fix.onset<datalength)->
      ftimes
    
    f<-rbindlist(apply(ftimes,1,jsac))
    
    f<- select(f,-event.dur)
    
    #the code below isn't very elegant, but it just combines the fixations and saccades
    #and assigns negative IDs to the fixations and positive IDs to the saccades
    f$issaccade=FALSE 
    g$issaccade=TRUE
    fg<-rbind(f,g)
    fg<- arrange(fg,time)
    
    fg$event[!fg$issaccade]=fg$event[!fg$issaccade]*-1
    
    #This is a debugging message in case the result isn't the correct length
    #we have to return a vector of the same length as the input
    if (length(fg$event)!=datalength){
      message('FAIL')
      message(length(fg$event))
      message(datalength)
      # message(paste('FAILED: ',length(fg$event,datalength,sep='-')))
    }
    else{
      # message('SUCCESS')
    }
    
    fg$event #return just an array of the IDs of saccades and fixations
  }else{
    g<- left_join(select(xx,time),g,by='time')
    return(g$event)
  }
  
}

markSaccadesDoubleTEST<- function(v, threshold1=60,threshold2=20,min.dur=5,maxreject=1000,
                              driftcorrect=FALSE,markFixations=TRUE){
  #This function is an R implementation of a two-threshold event marker
  #The algorithm works like this: Find all the times when velocity is above the high threshold
  #Extend this out until velocity is below the lower threshold
  
  #in practice, I'm identifying all the times that the saccades cross the low threshold and then rejecting 
  #any that don't meet the higher threshold
  #I'm also rejecting events that are below a certain duration
  
  #this algorithm also assigns positive ID numbers to the saccades and 
  #negative ID numbers to the non-saccades (fixations?)
  #after running this function, you can group_by(event) and measure the fixations or saccades as you wish
  require(dplyr)
  require(data.table) #for rbindlist - a fast version of do.call('rbind') that uses data.table
  datalength<-length(v)
  i<-which(abs(v)>threshold2) #find all the times when speed is above the lower threshold
  
  #For continuous numbers, diff=1. If there is a larger jump, it means there's a gap. 
  #That indicates another saccade
  sacoff<-which(diff(i)>1) #sacoff now contains the indices of the ends of all the saccades
  #sacon has all the indices of when the saccades start
  #After an offset, the next index will be the onset of the next saccade.
  #find the onsets by looking at the next index after the offset
  sacon<-c(1,sacoff+1) #first saccade always starts at first index
  sacoff<-c(sacoff,length(i)) #end of last saccade is always at the end
  event.onset<-i[sacon] #Convert from the indices to actual times
  event.offset<-i[sacoff] 
  
  #event.onset now has the time (in sample time) of all saccade onsets
  #set up stimes as a data.frame with two columns. Onset and offset. 
  stimes<- data.frame(event.onset,event.offset) 
  
  #this is a little function that works with the weirdness of R's "apply" family of functions
  #it just takes the onset and offset and returns the whole range. 
  #if you give it [10, 20], it will return [10 11 12 13 14 15 16 17 18 19 20]
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 4: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$event<- stimes[[4]]
    df$event.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes %>%
    mutate(dur=event.offset-event.onset, #calculate duration of each saccade
           s=row_number())-> #assign an ID number to each saccade
    stimes
  
  #Use "apply" to run the function "jsac" (above) on each line of the "stimes" data.frame
  #the result is the times of when all the saccades happen
  x<-rbindlist(apply(stimes,1,jsac)) 
  
  v<- data.frame(v=v) #Make the original velocity trace into a data.frame
  v<- mutate(v, time=row_number(),
             acc=parabolicdiff(v,4)) #add time to keep track
  
  #join the marked saccades and the velocity
  #the result is the velocity trace plus a row that just indicates whether there's a saccade
  #each saccade is identified by it's unique marker "event" that comes from df$event<- stimes[[4]] above 
  xx<- left_join(v,x,by='time')
  

  
  xx %>%
    group_by(event) %>% #This means we analyze each saccade individually
    summarize(max.vel=max(abs(v)), #calculate max velocity
              dur=n()) %>% #calculate duration
    filter(max.vel>threshold1, #reject all saccades that fail to exceed the large threshold
           max.vel<maxreject, #reject all saccades over the max threshold
           dur>min.dur)-> #reject all saccades that fail to exceed the minimum duration
    xm #xm is a summary which means it just lists the saccades and their measured values
  
  xx %>% #go back to the full data set and now reject all the saccades that were rejected above
    filter(event %in% unique(xm$event)) ->
    xx
  
  
  
  if (driftcorrect){
    xx%>%
      group_by(event) %>%
      summarize(realstart=first(time[v>50]),
                total.end=last(time),
                driftcut=first(time[time>(realstart+10)&abs(v)<100&abs(acc)<10000])) %>%
      mutate(total.cut=total.end-driftcut)->
      xs

    events= unique(xs$event)
    # manipulate(
    #   ggplot(filter(xx,event==events[chosenEvent]))+
    #     geom_line(aes(time,v))+
    #     geom_line(aes(time,acc/20),linetype=2)+
    #     geom_vline(xintercept = xs$driftcut[xs$event==events[chosenEvent]]),
    #   chosenEvent=slider(1,length(events),step=1))
    
    xx%>%
      group_by(event) %>%
      mutate(realstart=first(time[v>50]),
                total.end=last(time),
                driftcut=first(time[time>(realstart+10)&abs(v)<100&abs(acc)<5000]),
             rejectdrift=time>driftcut) %>%
      ungroup() %>%
      mutate(event=replace(event,rejectdrift,NA)) %>%
      select(time,event)->
      xx
  

  }
  
  
  
  xx %>% #go back to the full data set and now reject all the saccades that were rejected above
    # filter(event %in% unique(xm$event)) %>% 
    dplyr::select(time,event) -> #All we need is the time and the eventID
    g
  
  if (markFixations){
    #this next part goes through and assigns an ID to all the non-saccade portions of the data
    stimes2<- filter(stimes,s %in% unique(xm$event))
    
    ftimes<-data.frame(fix.onset=c(1,stimes2$event.offset+1),
                       fix.offset=c(stimes2$event.onset-1,datalength))
    
    ftimes %>%
      filter(fix.onset>0,fix.offset>0)%>%
      mutate(dur=fix.offset-fix.onset,
             s=row_number()) %>%
      filter(fix.onset<datalength)->
      ftimes
    
    f<-rbindlist(apply(ftimes,1,jsac))
    
    f<- select(f,-event.dur)
    
    #the code below isn't very elegant, but it just combines the fixations and saccades
    #and assigns negative IDs to the fixations and positive IDs to the saccades
    f$issaccade=FALSE 
    g$issaccade=TRUE
    fg<-rbind(f,g)
    fg<- arrange(fg,time)
    
    fg$event[!fg$issaccade]=fg$event[!fg$issaccade]*-1
    
    #This is a debugging message in case the result isn't the correct length
    #we have to return a vector of the same length as the input
    if (length(fg$event)!=datalength){
      message('FAIL')
      message(length(fg$event))
      message(datalength)
      # message(paste('FAILED: ',length(fg$event,datalength,sep='-')))
    }
    else{
      # message('SUCCESS')
    }
    
    fg$event #return just an array of the IDs of saccades and fixations
  }else{
    g<- left_join(select(xx,time),g,by='time')
    return(g$event)
  }
  
}

tt<- 