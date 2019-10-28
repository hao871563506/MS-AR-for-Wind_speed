
##load data
finaldf<-read.csv(file = "D:/columbia/research/californiawind/finalHMM_biggerwide.csv")
# 0 can not do boxcox, so we make it to 0.01
finaldf[finaldf == 0] <- 0.01

##take the month 3,4,5,6
library(lubridate)
finaldf$date<-as.Date(finaldf$date, format = "%m/%d/%y")
mjdftotal<- subset(finaldf, month(date) == '3' |  month(date) == '4' |  month(date) == '5' |  month(date) == '6')

#mjdftotal<-finaldf
numberofrecord<-nrow(mjdftotal)
#predict n observation
n<-0
##Now let's take out the last 11 to do our test
t1 = numberofrecord-n
t1plus1 = t1+1
mjdftest<-mjdftotal[t1plus1:numberofrecord,]
mjdftotal<-mjdftotal[1:t1,]

##BOXCOX
##Now let's change the whole dataframe
library(MASS)
mjdftotal2<-mjdftotal[,1:31]
lambdalist<-c()
for (i in 2:31){
  dummy = data.matrix(mjdftotal[i])
  Box = boxcox(dummy ~ 1,lambda = seq(-6,6,0.1))
  Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
  Cox2[1,]                                  # Display the lambda with the greatest log likelihood
  lambda = Cox2[1, "Box.x"]                 # Extract that lambda
  lambdalist<-c(lambdalist,lambda)
  T_box = (dummy ^ lambda - 1)/lambda   # Transform the original data
  mjdftotal2[i]<-T_box
}


#Do the NHMSAR 
library(NHMSAR)
mydata<- mjdftotal2[,2:31]
T = dim(mydata)[1]
d = dim(mydata)[2]
N.samples = 1 #dim(mydata)[3]
new.arr <- array( unlist(mydata), dim=c(T, N.samples,d) ) 

#number of regimes
M = 4
#order of AR processes
order = 2
theta.init = init.theta.MSAR(new.arr,M=M,order=order,label="HH")
mod.hh = fit.MSAR(new.arr,theta.init,verbose=TRUE,MaxIter=100)


##Simulate it
Bsim = 1000
Ksim = Bsim*N.samples
#Y0 = array(new.arr[t1minus1:t1,1,],c(2,Ksim,d))
Y0 = array(new.arr[1:2,1,],c(2,Ksim,d))
T = 8662
Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = Ksim)
###Box-cox back transfrom
for (i in 1:30){
  T_box = (dummy ^ lambda - 1)/lambda 
  Y.sim$Y[,,i]<-(Y.sim$Y[,,i]*lambdalist[i] + 1)^(1/lambdalist[i])
}
Y.sim$Y[is.nan(Y.sim$Y)]<-0.1
medianYsim<-Y.sim$Y[,1,]
for(location in 1:30){
  medianYsim[,location]<-apply(Y.sim$Y[,,location],1, median, na.rm = TRUE)
}


#calculate real 10% wind speed
mydata2<-mjdftotal[,2:31]
lowreal<-c()
highreal<-c()
for( location in 1:30){
  highlow<-quantile(mydata2[,location], probs = c(0.1, 0.9))
  lowreal<-c(lowreal,highlow[1])
  highreal<-c(highreal,highlow[2])
}

consecutiveObsmatrix<-matrix(0,4,30)
#consecutive in a row
for( location in 1:30){
thresholdObs <- which(mydata2[,location] < lowreal [location] )
result <- rle(diff(thresholdObs))
for(number in 2:5){
consecutiveObs<-sum(result$lengths>=number & result$values==1)
consecutiveObsmatrix[(number-1),location]<-consecutiveObs
}}




simulation<-Y.sim$Y
#consecutiveObsmatrix2<-matrix(0,4,30)
#for(i in 1:100){
#for( location in 1:30){
#thresholdObs <-which(simulation[,i,location]< lowreal [location] )
#result <- rle(diff(thresholdObs))
#for(number in 2:5){
#  consecutiveObs<-sum(result$lengths>=number & result$values==1)
#  consecutiveObsmatrix2[(number-1),location]<-consecutiveObsmatrix2[(number-1),location] + consecutiveObs
#}}}

#consecutiveObsmatrix
#round(consecutiveObsmatrix2/100,0)

number = 2
test=array(0,c(1000,30))
for(i in 1:1000){
  for( location in 1:30){
  thresholdObs <-which(simulation[,i,location]< lowreal [location] )
  result <- rle(diff(thresholdObs))
  consecutiveObs<-sum(result$lengths>=number & result$values==1)
  test[i,location]<-consecutiveObs
}}

lowtest<-c()
mediantest<-c()
hightest<-c()
for( location in 1:30){
  highlow<-quantile(test[,location], probs = c(0.05,0.5, 0.95))
  lowtest<-c(lowtest,highlow[1])
  mediantest<-c(mediantest,highlow[2])
  hightest<-c(hightest,highlow[3])
}

real<-consecutiveObsmatrix[(number-1),]
df<-data.frame(lowtest,mediantest,hightest,real)
df$dummy<-c(1:30)
library(reshape2)
dfm <- melt(df,id.vars = "dummy")
library(ggplot2)
p1<- ggplot(dfm,aes(x = dummy,y = value, width=0.75)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ xlab('location')+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))



testdf<-data.frame(test)
realdf<-data.frame(real)
colnames(testdf)<-c(1:30)
names<-colnames(testdf)
library(reshape)
meltData <- melt(testdf)
p1<- ggplot() + geom_boxplot(data = meltData, aes(x=variable,y=value),outlier.shape = NA)+theme_classic()+xlab('location')+ylab('occurance')+
  geom_point(data = realdf, aes(x=names,y=real),colour="red")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+scale_x_discrete(breaks=seq(0, 30, 5))














number = 3
test=array(0,c(1000,30))
for(i in 1:1000){
  for( location in 1:30){
    thresholdObs <-which(simulation[,i,location]< lowreal [location] )
    result <- rle(diff(thresholdObs))
    consecutiveObs<-sum(result$lengths>=number & result$values==1)
    test[i,location]<-consecutiveObs
  }}


lowtest<-c()
mediantest<-c()
hightest<-c()
for( location in 1:30){
  highlow<-quantile(test[,location], probs = c(0.05,0.5, 0.95))
  lowtest<-c(lowtest,highlow[1])
  mediantest<-c(mediantest,highlow[2])
  hightest<-c(hightest,highlow[3])
}

real<-consecutiveObsmatrix[(number-1),]
df<-data.frame(lowtest,mediantest,hightest,real)
df$dummy<-c(1:30)
library(reshape2)
dfm <- melt(df,id.vars = "dummy")
library(ggplot2)
p2<- ggplot(dfm,aes(x = dummy,y = value, width=0.75)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ xlab('location')+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))




testdf<-data.frame(test)
realdf<-data.frame(real)
colnames(testdf)<-c(1:30)
names<-colnames(testdf)
library(reshape)
meltData <- melt(testdf)
p2<- ggplot() + geom_boxplot(data = meltData, aes(x=variable,y=value),outlier.shape = NA)+theme_classic()+xlab('location')+ylab('occurance')+
  geom_point(data = realdf, aes(x=names,y=real),colour="red")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+scale_x_discrete(breaks=seq(0, 30, 5))










number = 4
test=array(0,c(1000,30))
for(i in 1:1000){
  for( location in 1:30){
    thresholdObs <-which(simulation[,i,location]< lowreal [location] )
    result <- rle(diff(thresholdObs))
    consecutiveObs<-sum(result$lengths>=number & result$values==1)
    test[i,location]<-consecutiveObs
  }}


lowtest<-c()
mediantest<-c()
hightest<-c()
for( location in 1:30){
  highlow<-quantile(test[,location], probs = c(0.05,0.5, 0.95))
  lowtest<-c(lowtest,highlow[1])
  mediantest<-c(mediantest,highlow[2])
  hightest<-c(hightest,highlow[3])
}

real<-consecutiveObsmatrix[(number-1),]
df<-data.frame(lowtest,mediantest,hightest,real)
df$dummy<-c(1:30)
library(reshape2)
dfm <- melt(df,id.vars = "dummy")
library(ggplot2)
p3<- ggplot(dfm,aes(x = dummy,y = value, width=0.75)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ xlab('location')+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))




testdf<-data.frame(test)
realdf<-data.frame(real)
colnames(testdf)<-c(1:30)
names<-colnames(testdf)
library(reshape)
meltData <- melt(testdf)
p3<- ggplot() + geom_boxplot(data = meltData, aes(x=variable,y=value),outlier.shape = NA)+theme_classic()+xlab('location')+ylab('occurance')+
  geom_point(data = realdf, aes(x=names,y=real),colour="red")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+scale_x_discrete(breaks=seq(0, 30, 5))








number = 5
test=array(0,c(1000,30))
for(i in 1:1000){
  for( location in 1:30){
    thresholdObs <-which(simulation[,i,location]< lowreal [location] )
    result <- rle(diff(thresholdObs))
    consecutiveObs<-sum(result$lengths>=number & result$values==1)
    test[i,location]<-consecutiveObs
  }}


lowtest<-c()
mediantest<-c()
hightest<-c()
for( location in 1:30){
  highlow<-quantile(test[,location], probs = c(0.05,0.5, 0.95))
  lowtest<-c(lowtest,highlow[1])
  mediantest<-c(mediantest,highlow[2])
  hightest<-c(hightest,highlow[3])
}

real<-consecutiveObsmatrix[(number-1),]
df<-data.frame(lowtest,mediantest,hightest,real)
df$dummy<-c(1:30)
library(reshape2)
dfm <- melt(df,id.vars = "dummy")
library(ggplot2)
p4<- ggplot(dfm,aes(x = dummy,y = value, width=0.75)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ xlab('location')+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))




testdf<-data.frame(test)
realdf<-data.frame(real)
colnames(testdf)<-c(1:30)
names<-colnames(testdf)
library(reshape)
meltData <- melt(testdf)
p4<- ggplot() + geom_boxplot(data = meltData, aes(x=variable,y=value),outlier.shape = NA)+theme_classic()+xlab('location')+ylab('occurance')+
  geom_point(data = realdf, aes(x=names,y=real),colour="red")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+scale_x_discrete(breaks=seq(0, 30, 5))








require(gridExtra)
grid.arrange(p1, p2,p3,p4,nrow=2, ncol=2)

