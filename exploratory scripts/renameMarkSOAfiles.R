
library(stringr)

path="C:/Users/setup/Desktop/NRTP Vergence/SOA Strabismus/"

#get name of all files
files <- list.files(path=path,pattern='*.csv')
#extract monkey name
monkeys<- lapply(files,function(x) str_split(x,'_')[[1]][1])
#unique monkey names
monkeys<- unique(unlist(monkeys))
#set up counters for each monkey 
cellnums<- (1:length(monkeys))*0+101
newFiles=list()

for (i in 1:length(files)){
  f<- files[[i]]
  f<- str_split(f,'_')
  f<- f[[1]]
  monkey<-f[1]
  f[1]=paste0(f[1],'-',cellnums[monkeys==monkey])
  cellnums[monkeys==monkey]=cellnums[monkeys==monkey]+1
  newFiles[[i]]=str_c(f,collapse='_')

}

for (i in seq_along(files)){
  
  file.rename(from = paste0(path,files[[i]]),
              to=paste0(path,newFiles[[i]]))
  
}

#THIS ACTUALLY WORKED 3-1-2018