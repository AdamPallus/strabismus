bootstrapSaccadesOK<- function(n){
  #n is just a label that is added. This function is run multiple times
  #using different random samples each time 
  get('z') #get from global environment
  #randomly pick from the list of saccades of each type. Replace=true means that
  #we can get the same number multiple times
  samp <- sample(unique(z$sacnum[z$verg.bins=='Converging']), 
                 length(unique(z$sacnum[z$verg.bins=='Converging'])), replace = TRUE)
  
  samp <- c(samp,sample(unique(z$sacnum[z$verg.bins=='Diverging']), 
                        length(unique(z$sacnum[z$verg.bins=='Diverging'])), replace = TRUE))
  #convert to data.table for faster processing
  z <- as.data.table(z)
  setkey(z, "sacnum") #like group_by
  # create the new data set
  z <- z[J(samp), allow.cartesian = TRUE] #replicate data set based on above sample 
  
  z %>% #make the model
    filter(abs(verg.velocity)<200) %>%
    do(tidy(lm(sdf20~verg.angle+verg.velocity,data=.))) %>%
    # do(tidy(lm(sdf20~lep+rep+lev+rev,data=.))) %>%
    mutate(repN=n)-> #add the number of the bootstrap iteration. 
    z
}
source('preparetoBOOTmarksaccs.R')
#step 1 is to run the first two chunks from SOArevision.Rmd
#that just loads the necessary libraries, helper functions and the main data file: ('SOA-NRTP.RDS')
nreps=1999 #number of bootstrap iterations
n<- matrix(1:nreps) #set this up in the proper form to work with apply below
# t<- readRDS('SOA.RDS')

neurons=unique(t$neuron)

xx<- NULL
for (i in 1:length(neurons)){
  # for(i in 1:2){
  message(paste('Processing: ',neurons[i]))
  z <- filter(t,neuron==neurons[i])
  z<-preparetoBOOTmarksaccs(z) #marks and measures saccades
  x<- as.data.frame(rbindlist(apply(n,1,bootstrapSaccadesOK))) #calls the function n times
  x<- mutate(x,neuron=neurons[i]) 
  xx[[i]]<- x #add to list for combining in next phase
}

xx<- rbindlist(xx) #combine efficiently (from data.table)

# saveRDS(xx,'SOA simple Bootstrap analysis 1999reps.RDS')
xx %>%
  group_by(neuron,term) %>%
  summarize(ciLOW=quantile(estimate,probs=0.025),
            ciHIGH=quantile(estimate,probs=0.975),
            m=mean(estimate)) %>%
  mutate(zerocross=(ciLOW*ciHIGH)<0)->
  confints

confints %>%
  select(neuron,term,m) %>%
  spread(term,m)->
  ciplot

confints %>%
  select(neuron,term,ciLOW) %>%
  spread(term,ciLOW)->
  ci.low
names(ci.low)<-paste(names(ci.low),'low',sep='.')

confints %>%
  select(neuron,term,ciHIGH) %>%
  spread(term,ciHIGH)->
  ci.high
names(ci.high)<-paste(names(ci.high),'high',sep='.')

ciplot<-cbind(ciplot,ci.low,ci.high)


ggplot(ciplot,aes(verg.velocity,verg.angle))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5)+
  geom_point(size=1,color='hotpink')+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')+
  geom_text(aes(label=neuron))

ggplot(ciplot,aes(verg.velocity,verg.angle))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5)+
  geom_point(size=1,color='hotpink')+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')


#remove coil noise from ozette-123
coilnoise<-c(65000:71000,
            80500:85000,
            95500:98500,
            101000:105500,
            151500:156000,
            187000:191500,
            226000:230500,
            242000:nrow(o123))

o123 %>%
  mutate(coil.noise=time %in% coilnoise)->
  o123
  