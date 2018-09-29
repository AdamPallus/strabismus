
#http://influentialpoints.com/Training/bootstrap_confidence_intervals.htm#bias

bootstrap.stat<- function(samp,z,chosenParam='sacnum'){
  
  #convert to data.table for faster processing
  z <- as.data.table(z)
  setkey(z, 'sacnum') #like group_by but in data.tables jargon
  # create the new data set
  z <- z[J(samp), allow.cartesian = TRUE] #replicate data set based on above sample 
  
  mod<- lm(sdf20~verg.angle+verg.velocity+conj.vertical,data=z)
  coef(mod)
}

bca.ci<-function(xx,ptb,chosenParam='sacnum'){
  
  # y<- unique(ptb$chosenParam)
  y<- unique(ptb[[chosenParam]])
  obs<-bootstrap.stat(y,ptb)
  
  #convert to tibble format (messy)
  
  obst<- as_tibble(obs)
  obst$term<-names(obs)
  obs<-select(obst,term,observation=value)
  
  xx %>% 
    select(term,estimate,repN,neuron) %>%
    left_join(obs,by='term')->
    boot
  
  
  r<-max(boot$repN)
  
  # estimate bias in std. norm deviates
  # b=qnorm((sum(boot > obs)+sum(boot==obs)/2)/r)
  
  boot %>% 
    group_by(term) %>%
    summarize(b=qnorm((sum(estimate > observation)+sum(estimate==observation)/2)/r)) %>%
    left_join(obs,.,by='term')->
    obs
  
  
  # estimate acceleration constant
  n=length(y)
  n1=n-1
  obsn=obs$observation*n
  # pv=i=0 ; while(i < n){i=i+1 ; pv[i]=obsn-n1*bootstrap.stat(y[-i],ptb,term)}
  
  pv<-tibble()
  for (i in 1:n){
    pv=bind_rows(pv,obsn-n1*bootstrap.stat(y[-i],ptb))
  }
  
  pv %>%
    gather() %>%
    rename(term=key) %>%
    group_by(term) %>%
    mutate(je=mean(value)-value) %>%
    summarize(a=sum(je^3)/(6*sum(je^2))^(3/2)) %>%
    left_join(obs,.,by='term')->
    obs
  # 
  # je=mean(pv)-pv
  # a=sum(je^3)/(6*sum(je^2))^(3/2)
  
  alpha=0.05 # 95% limits
  z=qnorm(c(alpha/2,1-alpha/2)) # Std. norm. limits
  
  terms<- unique(boot$term)
  q<-tibble()
  for (i in seq_along(terms)){
    obs_loop<- filter(obs,term==terms[i])
    b<- obs_loop$b
    a<- obs_loop$a
    p<-pnorm((z-b)/(1-a*(z-b))-b)
    bootterms<- filter(boot,term==terms[i])$estimate
    
    q.temp<-quantile(bootterms,p=p)
    names(q.temp)<-c('ciLOW','ciHIGH')
    q.temp$term=terms[i]
    q.temp$neuron=xx$neuron[1]
    q<-bind_rows(q,q.temp)
  }
  
  q<- left_join(q,select(obs,term,observation),by='term')
  
}



ptb<- readRDS('SOAnormalANDstrabismus.RDS')


xxBOTH<- readRDS('BootstrapSaccadesBoth3-12-2018.RDS')
neurons<- unique(xxBOTH$neuron)


ptb %>%
  filter(neuron %in% neurons) %>%
  group_by(neuron) %>%
  mutate(time=row_number()) %>%
  do(preparetoBOOTmarksaccs(.)) ->
  ptb

ptb %>%
  filter(saccade.type %in% c('left.only','right.only'))->
  ptb

ptb %>%
  group_by(neuron) %>%
  summarize(n=length(unique(sacnum)))->
  ptb.n

neurons<- filter(ptb.n,n>10)$neuron

cis<-list()

for (i in seq_along(neurons)){
  cis[[i]]<-bca.ci(filter(xxBOTH,neuron==neurons[i]),filter(ptb,neuron==neurons[i]))
}

cis<-rbindlist(cis)

cis %>%
  select(neuron,term,observation) %>%
  spread(term,observation)->
  ciplot

cis %>%
  select(neuron,term,ciLOW) %>%
  spread(term,ciLOW)->
  ci.low
names(ci.low)<-paste(names(ci.low),'low',sep='.')

cis %>%
  select(neuron,term,ciHIGH) %>%
  spread(term,ciHIGH)->
  ci.high
names(ci.high)<-paste(names(ci.high),'high',sep='.')

ciplot<-cbind(ciplot,ci.low,ci.high)

ciplot %>%
  separate(neuron,c('monkey','cellnum'),remove=FALSE) %>%
  mutate(strabismic=monkey %in% c('Kopachuck','Pilchuck','Patos','QT'),
         exo=monkey %in% c('Pilchuck','Patos'),
         eso= strabismic & !exo)->
  ciplot


getR2<-function(t,chosenParam='slowverg'){
  summary(lm(sdf20~verg.angle+verg.velocity+conj.vertical,data=t))$r.squared
}

ptb %>%
  group_by(neuron) %>%
  do(r2=getR2(.))->
  ptb.r2

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


ptb.r2$r2<-as.numeric(ptb.r2$r2)

citable<-left_join(ptb.n,ptb.r2,by='neuron')

citable<- left_join(citable,ciplot,by='neuron')

citable$neuron.low<-NULL
citable$neuron.high<-NULL

issignif<-function(x){
  s='none'
  velocity<-FALSE
  position<-FALSE
  if ((x$verg.velocity.low*x$verg.velocity.high)>0){
    velocity<-TRUE
    s='velocity'
  }
  if ((x$verg.angle.low*x$verg.angle.high)>0){
    position<-TRUE
    s='position'
  }
  if (velocity & position){
    s='both'
  }
  x$significant<-s
  return(x)
}


citable %>%
  filter(n>10) %>%
  group_by(neuron) %>%
  do(issignif(x=.))->
  citable

citable %>%
  group_by(significant) %>%
  summarize(n=n())

library(stringr)

nn<-names(citable)
names(citable)<- str_replace_all(nn,'[.]','_')
names(citable)<- nn

write.csv(citable,'SaccadeBootstrapCIs.csv')

ggplotly(
  ggplot(citable,aes(verg.velocity,verg.angle,label=neuron))+
    geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0,alpha=0.3)+
    geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5,alpha=0.3)+
    geom_point(aes(size=r2),color='hotpink',alpha=0.5)+
    geom_vline(xintercept = 0)+
    # coord_fixed()+
    # geom_text(aes(label=neuron))+
    theme_minimal()+
    xlab('Sensitivity to Vergence Velocity')+
    ylab('Sensitivity to Vergence Angle')
)

ggplot(citable,aes(verg.velocity,verg.angle,label=neuron))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0,alpha=0.3)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5,alpha=0.3)+
  geom_point(aes(size=r2),color='hotpink',alpha=0.5)+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')

ggplot(citable,aes(verg.velocity,verg.angle,label=neuron))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0,alpha=0.3)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5,alpha=0.3)+
  geom_point(aes(size=r2,color=significant),alpha=0.5)+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')
  # xlim(c(-2.88,4.61))+
  # ylim(c(-4.48,16.9))

ggplot(citableSLOW,aes(verg.velocity,verg.angle,label=neuron))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0,alpha=0.3)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5,alpha=0.3)+
  geom_point(aes(size=r2,color=significant),alpha=0.5)+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')

ggplot(citableSLOW,aes(verg.velocity,verg.angle,label=neuron))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0,alpha=0.3)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5,alpha=0.3)+
  geom_point(aes(color=significant),size=2,alpha=0.5)+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')

ggplot(citable,aes(verg.velocity,verg.angle,label=neuron))+
  geom_errorbar(aes(ymin=verg.angle.low,ymax=verg.angle.high),size=0.5,width=0,alpha=0.3)+
  geom_errorbarh(aes(xmin=verg.velocity.low,xmax=verg.velocity.high),height=0,size=0.5,alpha=0.3)+
  geom_point(aes(color=significant),size=2,alpha=0.5)+
  geom_vline(xintercept = 0)+
  # coord_fixed()+
  # geom_text(aes(label=neuron))+
  theme_minimal()+
  xlab('Sensitivity to Vergence Velocity')+
  ylab('Sensitivity to Vergence Angle')

