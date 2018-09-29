library(cladoRcpp)
library(dplyr)
library(data.table)

spikedensity<-function (rasters,sd=100) {
  gsize<- sd*10
  g<-dnorm(-gsize:gsize,mean=0,sd=sd)
  sdf<-rcpp_convolve(rasters,g)
  sdf<-sdf[gsize:(length(sdf)-(gsize+1))]*1000
  sdf
}

dynamiclead<-function(p,lags=seq(10,300,by=10),formula='rev+lev') {
  
  formula=paste('sdflag~',formula)
  
  rsq<-NULL
  for (i in 1:length(lags)) {
    if (lags[i] > 0){
      p$sdflag<-dplyr::lag(p$sdf,lags[i])
    }
    else{
      p$sdflag<-dplyr::lead(p$sdf,lags[i]*-1)
    }
    
    rsq[i]<- summary(lm(formula=formula,data=p))$r.squared
  }
  #return(rsq)
  # return(lags[rsq==max(rsq)])
  bestlag=lags[rsq==max(rsq)]
  p$dynamiclead<- bestlag
  if (bestlag > 0){
    p$sdflag<-dplyr::lag(p$sdf,bestlag)
  }
  else{
    p$sdflag<-dplyr::lead(p$sdf,bestlag*-1)
  }
  return(p)
  
}

findSaccades<-function(ev,threshold=40){
  mindur<-1
  i<-which(abs(ev)>threshold) #find all the times when speed > threshold
  sacoff<-which(diff(i)>mindur) #minimum duration of an accepted saccade
  sacon<-c(1,sacoff+1) #first saccade
  sacoff<-c(sacoff,length(i)) #end of last saccade
  saccade.onset<-i[sacon] #get actual times
  saccade.offset<-i[sacoff] 
  return(data.frame(saccade.onset,saccade.offset))
}

markSaccades<-function(ev,buffer=15,threshold=40){
  #this function finds and marks saccades given a velocity input
  stimes<-findSaccades(ev,threshold)
  
  #remove saccades without enough data at the end of the file, based on buffer size
  toolong<- stimes$saccade.offset> length(ev)-buffer
  tooshort<- stimes$saccade.onset<buffer+1
  stimes<- filter(stimes, !tooshort, !toolong)
  
  nsaccades=nrow(stimes)
  
  stimes$saccade.onset=stimes$saccade.onset-buffer
  stimes$saccade.offset=stimes$saccade.offset+buffer
  
  s<-1:length(ev)*0
  
  for (k in 1:nsaccades){
    s[stimes$saccade.onset[k]:stimes$saccade.offset[k]]<- k
    if(k>1){
      s[stimes$saccade.offset[k-1]:stimes$saccade.onset[k]]<-(k*-1)
    }
  }
  s[1:stimes$saccade.onset[1]]<- -1
  s[stimes$saccade.offset[nrow(stimes)]:length(s)]<- (nrow(stimes)*-1)-1
  return(s)
}

parabolicdiff <- function(pos,n=7){
  q <- sum(2*((1:n)^2))
  convoutput<- rcpp_convolve(pos,c(-n:-1, 1:n))
  convoutput<- convoutput[(n*2):(length(pos)-((n*2)+1))]
  vels<- c(array(convoutput[1],dim=n*2),convoutput,array(convoutput[length(convoutput)],dim=n*2))
  vels <- vels/q*-1000
}

maxabs<- function(x){
  m1<-max(x,na.rm=T)
  m2<-min(x,na.rm=T)
  if (abs(m1)>abs(m2)) {
    return(m1)
  } else{
    return(m2)
  }
}


minabs<- function(x){
  m1<-max(x,na.rm=T)
  m2<-min(x,na.rm=T)
  if (abs(m1)>abs(m2)) {
    return(m2)
  } else{
    return(m1)
  }
}


loadnewcsv<- function(referencefile=NULL,path="C:/Users/setup/Desktop/NRTP Vergence/",
                      keeptarget=FALSE){
  require(stringr)
  require(dplyr)
  #This function loads .csv files in a particular folder. They must have the same columns for rbind
  #Saves time by only reading the csv when necessary
  
  #get names of all files in path
  files <- list.files(path=path,pattern='*.csv')
  #extract neuron name eg. Bee-01
  names<-sapply(files, str_match,"^[a-zA-Z]+-[0-9]+",USE.NAMES=FALSE)
  # check for new cells
  if (!is.null(referencefile)){
    files<-files[!names %in% referencefile$neuron] #comparison
  }
  
  nfiles<-length(files)
  
  if (nfiles>0){
    message(c('New Files: ',files))
    loadedfiles <- lapply(paste(path,files,sep=''),read.csv)
    t<-data.frame()
    # t.old<-NULL
    for (i in 1:nfiles) {
      f<- files[i]
      message(paste('Loading:',f))
      temp=loadedfiles[[i]]
      names<-str_match(f,"(^[a-zA-Z]+)-([0-9]+)")
      temp$neuron<-names[1]
      temp$monkey<-names[2]
      temp$cellnum<-as.numeric(names[3])
      temp$sdf<-spikedensity(temp$rasters,sd=10)
      # leadtime<-dynamiclead(temp)
      message(paste('ms of data:',nrow(temp)))
      mutate(temp,
             # sdflag=lag(sdf,leadtime),
             conj.velocity=sqrt((rev^2+lev^2)/2)+sqrt((revV^2+levV^2)/2),
             # s=markSaccades(conj.velocity,buffer=10,threshold=10),
             # slong=markSaccades(conj.velocity,buffer=longbuffer,threshold=10),
             time=row_number(),
             verg.angle=lep-rep,
             # verg.velocity=parabolicdiff(verg.angle,7),
             verg.velocity=lev-rev)->
        temp
      
      t <-rbind(t,temp)
    }
    if (!keeptarget){
      t<- dplyr::select(t, -thp,-tvp,-time)
      if('tvp2' %in% names(t)){
        t<- dplyr::select(t,-thp2,tvp2)
      }
    }
  }else{
    message('********NO NEW CELLS********')
    t<-NULL
  }
  return(t)
}

joinsaccadesuniformOLD<-function(t,buffer=20,threshold=30,saccade.length=150){
  
  findSaccades<-function(ev,threshold=40){
    mindur<- 1
    i<-which(abs(ev)>threshold) #find all the times when speed > threshold
    sacoff<-which(diff(i)>mindur) #minimum duration of an accepted saccade
    sacon<-c(1,sacoff+1) #first saccade
    sacoff<-c(sacoff,length(i)) #end of last saccade
    saccade.onset<-i[sacon] #get actual times
    saccade.offset<-i[sacoff] 
    
    return(data.frame(saccade.onset,saccade.offset))
  }
  
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 2: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$sacnum<- stimes[[4]]
    df$saccade.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes<-findSaccades(t$conj.velocity,threshold)
  stimes %>%
    mutate(dur=saccade.offset-saccade.onset,
           s=row_number(),
           saccade.onset=saccade.onset-buffer,
           ###########HERE IS WHERE i MAKE SACCADES UNIFORM #######
           saccade.offset=saccade.onset+saccade.length+2*buffer)->
    stimes
  
  x<- do.call('rbind',apply(stimes,1,jsac))
  x %>%
    group_by(sacnum) %>%
    mutate(counter=time-first(time)) ->
    x
  left_join(t,x,by='time')
}

joinsaccades<-function(t,buffer=20,threshold=30){
  findSaccades<-function(ev,threshold=40){
    mindur<- 1
    i<-which(abs(ev)>threshold) #find all the times when speed > threshold
    sacoff<-which(diff(i)>mindur) #minimum duration of an accepted saccade
    sacon<-c(1,sacoff+1) #first saccade
    sacoff<-c(sacoff,length(i)) #end of last saccade
    saccade.onset<-i[sacon] #get actual times
    saccade.offset<-i[sacoff] 
    
    return(data.frame(saccade.onset,saccade.offset))
  }
  
  jsac<- function(stimes){
    #input should be an array of length 2: c(onsettime,offsettime)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$sacnum<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes<-findSaccades(t$conj.velocity,threshold)
  stimes %>%
    mutate(s=row_number(),
           saccade.onset=saccade.onset-buffer,
           saccade.offset=saccade.offset+buffer)->
    stimes
  x<- as.data.frame(rbindlist(apply(stimes,1,jsac)))
  # x<- do.call('rbind',apply(stimes,1,jsac))
  x %>%
    group_by(sacnum) %>%
    mutate(counter=time-first(time)) ->
    x
  left_join(t,x,by='time')
}

markEnhancement<- function(v, threshold1=15,threshold2=8){
  require(dplyr)
  
  mindur<- 1
  i<-which(abs(v)>threshold2) #find all the times when speed is above the lower threshold
  sacoff<-which(diff(i)>mindur) #minimum duration of an accepted saccade
  sacon<-c(1,sacoff+1) #first saccade
  sacoff<-c(sacoff,length(i)) #end of last saccade
  event.onset<-i[sacon] #get actual times
  event.offset<-i[sacoff] 
  
  stimes<- data.frame(event.onset,event.offset)
  nsaccades=nrow(stimes)
  
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 2: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$enhancenum<- stimes[[4]]
    df$enhance.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes %>%
    mutate(dur=event.offset-event.onset,
           s=row_number())->
    stimes
  
  x<- do.call('rbind',apply(stimes,1,jsac))
  
  
  v<- data.frame(v=v)
  v<- mutate(v, time=row_number())
  
  xx<- left_join(v,x,by='time')
  
  xx %>%
    group_by(enhancenum) %>%
    summarize(max.vel=max(abs(v))) %>%
    filter(max.vel>threshold1)->
    xm
  
  xx %>%
    filter(enhancenum %in% unique(xm$enhancenum)) %>%
    dplyr::select(time,v,enhancenum) ->
    g
  
  
}

makeRelImp<- function(summaryforplot,formula='mean.Spikerate~mean.Verg.Angle+mean.C.Ver',normalize=TRUE){
  require(tidyr)
  summaryforplot %>% 
    group_by(neuron) %>%
    do(m=lm(formula=formula,data=.))-> 
    mm
  
  r<- data.frame()
  r2 <- NULL
  if (length(mm$m[[1]]$coefficients)<3){
    message('Not enough parameters for importance calculation...')
    message('Returning r.squared only.')
    for (i in 1:nrow(mm)){
      r2[i]<- summary(mm$m[[i]])$r.squared
    }
    r<- data.frame(neuron=mm$neuron,R2=r2)
    r<- separate(r,neuron,c('monkey','cellnum'),remove=FALSE)
    return(r)
  }  
  for (i in 1:nrow(mm)){
    # message('Trying...')
    bb<- relaimpo::calc.relimp(mm$m[[i]],rela=normalize)
    b<- bb$lmg
    r<-rbind(r,b)
    r2<- c(r2,bb$R2)
    
  }
  r$neuron<-mm$neuron
  r$R2<- r2
  names(r)<- c(names(b),'neuron','R2')
  r<- separate(r,neuron,c('monkey','cellnum'),remove=FALSE)
  r
}

joinsaccadesuniform<-function(t,buffer=20,threshold=30,saccade.length=150){
  
  findSaccades<-function(ev,threshold=40){
    mindur<- 1
    i<-which(abs(ev)>threshold) #find all the times when speed > threshold
    sacoff<-which(diff(i)>mindur) #minimum duration of an accepted saccade
    sacon<-c(1,sacoff+1) #first saccade
    sacoff<-c(sacoff,length(i)) #end of last saccade
    saccade.onset<-i[sacon] #get actual times
    saccade.offset<-i[sacoff] 
    
    return(data.frame(saccade.onset,saccade.offset))
  }
  
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 2: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$sacnum<- stimes[[4]]
    df$saccade.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes<-findSaccades(t$conj.velocity,threshold)
  stimes %>%
    mutate(dur=saccade.offset-saccade.onset,
           s=row_number(),
           saccade.onset=saccade.onset-buffer,
           ###########HERE IS WHERE i MAKE SACCADES UNIFORM #######
           saccade.offset=saccade.onset+saccade.length+2*buffer)->
    stimes
  
  x<- as.data.frame(rbindlist(apply(stimes,1,jsac)))
  x %>%
    group_by(sacnum) %>%
    mutate(counter=time-first(time)) ->
    x
  left_join(t,x,by='time')
}

dynamiclead2<-function(p,lags=seq(20,80,by=2),formula='verg.velocity~sdflag+verg.angle') {
  #rewrote this function to work with variables other than sdflag as the prediction
  rsq<-NULL
  for (i in 1:length(lags)) {
    if (lags[i] > 0){
      p$sdflag<-dplyr::lag(p$sdf,lags[i])
    }
    else{
      p$sdflag<-dplyr::lead(p$sdf,lags[i]*-1)
    }
    
    rsq[i]<- summary(lm(formula=formula,data=p))$r.squared
  }
  # return(rsq)
  return(lags[rsq==max(rsq)])
  # bestlag=lags[rsq==max(rsq)]
  
}

loadnewcsv2<- function(referencefile=NULL,path="C:/Users/setup/Desktop/NRTP Vergence/"){
  require(stringr)
  require(dplyr)
  #This function loads .csv files in a particular folder. They must have the same columns for rbind
  #Saves time by only reading the csv when necessary
  
  #get names of all files in path
  files <- list.files(path=path,pattern='*.csv')
  #extract neuron name eg. Bee-01
  names<-sapply(files, str_match,"^[a-zA-Z]+-[0-9]+",USE.NAMES=FALSE)
  # check for new cells
  if (!is.null(referencefile)){
    files<-files[!names %in% referencefile$neuron] #comparison
  }
  
  nfiles<-length(files)
  
  if (nfiles>0){
    message(c('New Files: ',files))
    loadedfiles <- lapply(paste(path,files,sep=''),read.csv)
    t<-data.frame()
    temp<- NULL
    # t.old<-NULL
    for (i in 1:nfiles) {
      f<- files[i]
      message(paste('Loading:',f))
      # temp[[i]]=loadedfiles[[i]]
      names<-str_match(f,"(^[a-zA-Z]+)-([0-9]+)")
      loadedfiles[[i]]$neuron<-names[1]
      loadedfiles[[i]]$monkey<-names[2]
      loadedfiles[[i]]$cellnum<-as.numeric(names[3])
      # loadedfiles[[i]]$sdf<-spikedensity(temp$rasters,sd=10)
      # leadtime<-dynamiclead(temp)
      message(paste('ms of data:',nrow(temp)))
      mutate(loadedfiles[[i]],
             # sdflag=lag(sdf,leadtime),
             conj.velocity=sqrt(((rev+lev)/2)^2+((revV+levV)/2)^2),
             # s=markSaccades(conj.velocity,buffer=10,threshold=10),
             # slong=markSaccades(conj.velocity,buffer=longbuffer,threshold=10),
             time=row_number(),
             verg.angle=lep-rep,
             # verg.velocity=parabolicdiff(verg.angle,7),
             verg.velocity=lev-rev)->
        loadedfiles[[i]]
      
      
    }
    t <-rbindlist(loadedfiles,fill=TRUE)
    # t<- dplyr::select(t, -thp,-tvp,-time)
  }else{
    message('********NO NEW CELLS********')
    t<-NULL
  }
  return(t)
}

markSaccadesDoubleOLD<- function(v, threshold1=60,threshold2=20,min.dur=5,maxreject=1000,
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
    xx %>%
      ungroup() %>%
      mutate(acc=parabolicdiff(v,7)) %>%
      mutate(event=replace(event,abs(v)<100 & abs(acc)<10000,NA)) ->
      xx
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
markSaccadesDouble<- function(v, threshold1=60,threshold2=20,min.dur=5,maxreject=1000,
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
    #'The way that drift correct works is that it finds when the velocity 
    #'exceeds the higher threshold and then look for the first time after that
    #'when the velocity is below 100 deg/s AND the acceleration is below 10,000
    #'This comes from Mark's published work. He determined these values by comparing
    #'the firing of MLBs and matching the duration of the burst to the duration of the saccades
    #'The reason for this correction is that monkeys with strabismus often show a post-saccadic
    #'drift where the eyes seem to re-accelerate after the end of the saccade which can cause
    #'velocity-only algorithms to overestimate the duration of the saccade. 
    xx%>%
      mutate(event2=event, #you're not allowed to change grouping variables, so make a copy first
             acc=parabolicdiff(v,7)) %>%
      group_by(event2) %>%
      mutate(realstart=first(time[v>threshold1]),
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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}