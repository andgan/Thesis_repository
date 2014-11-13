## analysis

library(VIM)
num <- as.numeric(save[[1]][,1])
nam <- save[[1]][,2]
aggregate(num, list(nam), countNA)
GOF <- aggregate(num, list(nam), mean, na.rm=T)


expr <- aggregate(save[[1]],list(save[[1]][,6]),mean, na.rm=T)
obsr <- aggregate(save[[2]],list(save[[2]][,6]),mean, na.rm=T)
mean(data.frame(save[[3]]), na.rm=T)


expr <- aggregate(save[[4]],list(save[[4]][,6]),mean, na.rm=T)
obsr <- aggregate(save[[5]],list(save[[5]][,6]),mean, na.rm=T)
mean(data.frame(save[[6]]), na.rm=T)

expr <- aggregate(save[[7]],list(save[[7]][,6]),mean, na.rm=T)
obsr <- aggregate(save[[8]],list(save[[8]][,6]),mean, na.rm=T)
mean(data.frame(save[[9]]), na.rm=T)


#### CLINICAL CLASSESS
### QUINTILES
expr <- NULL
obsr <- NULL
for (i in seq(1,nrow(save[[4]]),by=3)){ 
  exp <- save[[4]][i:(i+2),]
  obs <- save[[5]][i:(i+2),]	
  obsr <- rbind(obsr,colSums(obs))
  expr <- rbind(expr,colSums(exp))
}




### QUINTILES
expr <- NULL
obsr <- NULL
for (i in 1:seq(1,nrow(save[[1]]),by=5))
{ exp <- save[[1]][i:(i+4),]
  obs <- save[[2]][i:(i+4),]	
  obsr <- c(obsr,obs)
  expr <- c(expr,exp)
}


### QUINTILES
expr <- NULL
obsr <- NULL
for (i in seq(1,nrow(save[[7]]),by=5))
{ exp <- save[[7]][i:(i+4),]
  obs <- save[[8]][i:(i+4),]	
  obsr <- rbind(obsr,colSums(obs))
  expr <- rbind(expr,colSums(exp))
}
