removeBlink<-function(x,buffer=5){
  #This function finds periods of missing data
  #and removes additional points nearby, defined by buffer
  
  #This is necessary beacuse there is a period where data are not missing, 
  #but they are not accurate. 
  
  
  i<-which(is.na(x)) #find all the times when the data are missing
  sacoff<-which(diff(i)>1)
  sacon<-c(1,sacoff+1) #first saccade
  sacoff<-c(sacoff,length(i)) #end of last saccade
  missing.onset<-i[sacon] #get actual times
  missing.offset<-i[sacoff] 
  
  blinks<- data.frame(onset=missing.onset,offset=missing.offset)
  
  xnew<- x
  
  blinks %>%
    mutate(expanded.onset=onset-buffer,
           expanded.offset=offset+buffer,
           expanded.onset=replace(expanded.onset,expanded.onset<0,0),
           expanded.offset=replace(expanded.offset,expanded.offset>length(x),length(x)))->
    blinks
  
  for (k in 1:nrow(blinks)){
    xnew[blinks$expanded.onset[k]:blinks$expanded.offset[k]]<- NA
  }  
  
  xnew
}