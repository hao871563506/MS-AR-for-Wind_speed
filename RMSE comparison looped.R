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

#one year is 122 days, we take out one year for test
n<-122
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
#1 year 122 days
T = 122
#T = dim(mydata)[1]
d = dim(mydata)[2]
#N.samples = 1 
#70 years of samples
N.samples = 70
new.arr <- array( unlist(mydata), dim=c(T, N.samples,d) ) 
#number of regimes
M = 4
#order of AR processes
order = 2
theta.init = init.theta.MSAR(new.arr,M=M,order=order,label="HH")
mod.hh = fit.MSAR(new.arr,theta.init,verbose=TRUE,MaxIter=100)
#mod.hh = fit.MSAR(new.arr,theta.init,verbose=TRUE,MaxIter=100,penalty = "SCAD",lambda1=.5,lambda2=.5)

N.samples2 = 1
test.arr <- array( unlist(mjdftest[,2:31]), dim=c(n, N.samples2,d) ) 

RMSEwithoutsum = function(m, o){
  sqrt((m - o)^2)
}

###########Use cross-validation to calcualte RMSE matrix
RMSEMSARmatrix <- matrix(0, 5,30)#2 to 21 days prediction
for(lengthforsimule in 4:8){
for(timeforsimule in 1:100){
  Bsim = 100
  N.samples = 1
  Ksim = Bsim*N.samples
  timeforsimuleplus1 = timeforsimule + 1
  Y0 <- array(rep(0, 2*Ksim*d), dim=c(2, Ksim, d))
  for(i in 1:Ksim){
    Y0[,i,]<-test.arr[timeforsimule:timeforsimuleplus1,1,]
  }
  #Y0 = array(test.arr[timeforsimule:timeforsimuleplus1,1,],c(2,Ksim,d))
  T = lengthforsimule
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

  lengthforsimuleminus2 = lengthforsimule -2
  predictdf = data.frame(V1 = medianYsim[3:lengthforsimule,])
  RMSEMSARsum = rep(0, 30)
  for(i in 1:lengthforsimuleminus2){
    predict<-predictdf[i,]
    predict<-as.numeric(as.vector(predict))
    actual<-mjdftest[timeforsimuleplus1+i,][2:31]
    actual<-as.numeric(as.vector(actual))
    RMSEMSAR = RMSEwithoutsum(actual,predict)
    RMSEMSARsum = RMSEMSARsum + RMSEMSAR
  }
  RMSEMSARmatrix[lengthforsimule-3,] <- RMSEMSARmatrix[lengthforsimule-3,]+RMSEMSARsum
}}
RMSEMSARmatrix<-RMSEMSARmatrix/100
#########This is because we take sum, so we need to average them
for(i in 2:6){
RMSEMSARmatrix[i-1,]<-RMSEMSARmatrix[i-1,]/i
}




#write.csv(RMSEMSARmatrix, file = "D:/columbia/research/californiawind/RMSEMSARmatrix.csv",row.names=FALSE)

library(Rfast)
Mymatrix<-as.matrix(mjdftest[,2:31])
Varianceofdata <- colVars(Mymatrix)
Varianceofdata
#######Calculate R square
VairanceofMSAR<-RMSEMSARmatrix^2
R2MSARdf<-as.data.frame(1 - sweep(VairanceofMSAR, 2, Varianceofdata, `/`))
R2MSARdf

RMSEMSARdf<-as.data.frame(RMSEMSARmatrix)

R2dfheatmap <- as.matrix(R2MSARdf[1,])
R2matrixheatmap <- matrix(0,6,5)
for(i in 0:5){
  for( j in 1:5){
    R2matrixheatmap[6-i,j]<-R2dfheatmap[5*i+j]
  }
}

library(reshape2)
library(ggplot2)
ggplot(melt(R2matrixheatmap), aes(Var1,Var2, fill=value)) + geom_raster()+scale_fill_gradient(low = "white", high = "#00AFBB")






library(ggplot2)
library(reshape2)
RMSEMSARdf$group <- row.names(RMSEMSARdf)
#RMSEMSARdf$group <- factor(RMSEMSARdf$group, levels = c('1', '2', '3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'))
RMSEMSARdf$group <- factor(RMSEMSARdf$group, levels = c('1', '2', '3','4','5'))
RMSEMSARdf.m <- melt(RMSEMSARdf, id.vars = "group")

ggplot(RMSEMSARdf.m, aes(group, value)) + geom_boxplot()+theme_classic()+xlab('days in prediction')+ylab('Root mean squarae error of prediction')+
  #ylim(0.6,3)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))




#################
##This one is hmm with order = 0 
###############





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
n<-122
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
order = 0
theta.init = init.theta.MSAR(new.arr,M=M,order=order,label="HH")
mod2.hh = fit.MSAR(new.arr,theta.init,verbose=TRUE,MaxIter=100)


test.arr <- array(unlist(mjdftest[,2:31]), dim=c(n, N.samples,d) ) 

RMSEwithoutsum = function(m, o){
  sqrt((m - o)^2)
}

RMSEhmmmatrix <- matrix(0, 5,30)#2 to 21 days prediction
for(lengthforsimule in 4:8){
  for(timeforsimule in 1:100){
    Bsim = 100
    Ksim = Bsim*N.samples
    timeforsimuleplus1 = timeforsimule + 1
    Y0 <- array(rep(0, 2*Ksim*d), dim=c(2, Ksim, d))
    for(i in 1:Ksim){
      Y0[,i,]<-test.arr[timeforsimule:timeforsimuleplus1,1,]
    }
    #Y0 = array(test.arr[timeforsimule:timeforsimuleplus1,1,],c(2,Ksim,d))
    T = lengthforsimule
    Y.sim = simule.nh.MSAR(mod2.hh$theta,Y0 = Y0,T,N.samples = Ksim)
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
    lengthforsimuleminus2 = lengthforsimule -2
    predictdf = data.frame(V1 = medianYsim[3:lengthforsimule,])
    RMSEhmmsum = rep(0, 30)
    for(i in 1:lengthforsimuleminus2){
      predict<-predictdf[i,]
      predict<-as.numeric(as.vector(predict))
      actual<-mjdftest[timeforsimuleplus1+i,][2:31]
      actual<-as.numeric(as.vector(actual))
      RMSEhmm = RMSEwithoutsum(actual,predict)
      RMSEhmmsum = RMSEhmmsum + RMSEhmm
    }
    RMSEhmmmatrix[lengthforsimule-3,] <- RMSEhmmmatrix[lengthforsimule-3,]+RMSEhmmsum
  }}
RMSEhmmmatrix<-RMSEhmmmatrix/100

for(i in 2:6){
  RMSEhmmmatrix[i-1,]<-RMSEhmmmatrix[i-1,]/i
}

RMSEhmmmatrix

#write.csv(RMSEhmmmatrix, file = "D:/columbia/research/californiawind/RMSEhmmmatrix.csv",row.names=FALSE)

RMSEhmmdf<-as.data.frame(RMSEhmmmatrix)


Vairanceofhmm<-RMSEhmmmatrix^2
R2hmmdf<-as.data.frame(1 - sweep(Vairanceofhmm, 2, Varianceofdata, `/`))
R2hmmdf




library(ggplot2)
library(reshape2)
RMSEhmmdf$group <- row.names(RMSEhmmdf)
#RMSEhmmdf$group <- factor(RMSEhmmdf$group, levels = c('1', '2', '3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'))
RMSEhmmdf$group <- factor(RMSEhmmdf$group, levels = c('1', '2', '3','4','5'))
RMSEhmmdf.m <- melt(RMSEhmmdf, id.vars = "group")

ggplot(RMSEhmmdf.m, aes(group, value)) + geom_boxplot()+theme_classic()+xlab('days in prediction')+ylab('average RMSE per day')+
  ylim(0.6,3)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))








############
#persistence






RMSEpersistencematrix <- matrix(0, 5,30)#2 to 21 days prediction

for(lengthforsimule in 4:8){
  for(timeforsimule in 1:100){
    lengthforsimuleminus2 = lengthforsimule -2
    predict<-mjdftest[timeforsimule,][2:31]
    predict<-as.numeric(as.vector(predict))
    RMSEpersistencesum = rep(0, 30)
    for(i in 2:lengthforsimuleminus2){
      actual<-mjdftest[timeforsimuleplus1+i,][2:31]
      actual<-as.numeric(as.vector(actual))
      RMSEpersistence = RMSEwithoutsum(actual,predict)
      RMSEpersistencesum = RMSEpersistencesum + RMSEpersistence
    }
    RMSEpersistencematrix[lengthforsimule-3,] <- RMSEpersistencematrix[lengthforsimule-3,]+RMSEpersistencesum
  }}
RMSEpersistencematrix<-RMSEpersistencematrix/100

for(i in 2:6){
  RMSEpersistencematrix[i-1,]<-RMSEpersistencematrix[i-1,]/i
}

RMSEpersistencematrix

#write.csv(RMSEpersistencematrix, file = "D:/columbia/research/californiawind/RMSEpersistencematrix.csv",row.names=FALSE)
R2persistencedf<-as.data.frame(1-RMSEpersistencematrix^2/Varianceofdata)
#R2persistencedf[R2persistencedf < 0] <- 0






#############
###Comparison 





MSARmedian<-as.numeric(as.vector(apply(R2MSARdf,1, mean, na.rm = TRUE)))
#hmmmedian<-as.numeric(as.vector(apply(R2hmmdf,1, mean, na.rm = TRUE)))
persistencemedian<-as.numeric(as.vector(apply(R2persistencedf,1, mean, na.rm = TRUE)))


Iter_name <- "MSAR"
dummy_name <- "persistence"
#viter_name <- "HMM"
#df <- data.frame(MSARmedian,persistencemedian,hmmmedian)
df <- data.frame(MSARmedian,persistencemedian)
#colnames(df) <- c(Iter_name, dummy_name,viter_name)
colnames(df) <- c(Iter_name, dummy_name)
df$number<- 1:5
library(ggplot2)
ggplot(df, aes(number)) + 
  geom_line(aes(y = MSAR, colour = "MSAR"),size=2) + 
  geom_line(aes(y = persistence, colour = "persistence"),size=2) + 
  #geom_line(aes(y = HMM, colour = "HMM"),size=2)+
  xlab("prediction of days ahead")+ylab("R square of prediction")+ 
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+ scale_color_discrete(name = "Method")
