
xx<-readRDS('BootstrapSaccadesBoth3-12-2018.RDS')
t %>%
  filter(neuron=='Bee-101') %>%
  group_by(neuron) %>%
  do(preparetoBOOTslow(.)) ->
  ptb



bootstrap.stat<- function(samp,z){
  
  #convert to data.table for faster processing
  z <- as.data.table(z)
  setkey(z, "slowverg") #like group_by but in data.tables jargon
  # create the new data set
  z <- z[J(samp), allow.cartesian = TRUE] #replicate data set based on above sample 
  
  mod<- lm(sdf20~verg.angle+verg.velocity+conj.vertical,data=z)
  coef(mod)[['verg.velocity']]
}



bca.ci<- function(ptb,R=199){
  y<- unique(ptb$slowverg)
  obs<-bootstrap.stat(y,ptb)
  # R=1999 # use 5000 bootstraps
  boot=r=0 ; while(r < R){r=r+1 # for each of R samples
  boot[r]=bootstrap.stat(sample(y,replace=TRUE),ptb)} # get bootstrap stat.
  
  
  
  # estimate bias in std. norm deviates
  b=qnorm((sum(boot > obs)+sum(boot==obs)/2)/r)
  
  # estimate acceleration constant
  n=length(y) ; n1=n-1 ; obsn=obs*n
  pv=i=0 ; while(i < n){i=i+1 ; pv[i]=obsn-n1*bootstrap.stat(y[-i],ptb)}
  je=mean(pv)-pv
  a=sum(je^3)/(6*sum(je^2))^(3/2)
  
  alpha=0.05 # 95% limits
  z=qnorm(c(alpha/2,1-alpha/2)) # Std. norm. limits
  p=pnorm((z-b)/(1-a*(z-b))-b) # correct & convert to proportions
  
  quantile(boot,p=p) # ABC percentile lims.
}

bca.ci(ptb,19)

