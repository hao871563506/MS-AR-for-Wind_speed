RMSEwithoutsum = function(m, o){
  sqrt((m - o)^2)
}

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
mydata<-mjdftotal[,2:31]



aa<-as.matrix(mydata)%*%eigen(cor(mydata))$vectors ##as.matrix(mydata)%*%ae1
ae<-eigen(cor(mydata))$vectors

pca <- prcomp(mydata, scale = FALSE,center=FALSE)
summary(pca)
#aa1<-prcomp(mydata, scale = FALSE,center=FALSE)$x
#ae1<-prcomp(mydata, scale = FALSE,center=FALSE)$rotation
#b<-as.matrix(mydata)%*%ae1
#b[1:10,1:3]
b<-aa%*%solve(ae)
b[1:10,1:3]
  
library(forecast)
armafit1 <- auto.arima(aa[,1])
arma1<-armafit1$x
armafit2 <- auto.arima(aa[,2])
arma2<-armafit2$x
armafit3 <- auto.arima(aa[,3])
arma3<-armafit3$x
fitaa <- data.frame(arma1, arma2, arma3)
fitaa <- as.matrix(fitaa)
library(MASS)
inv <- ginv(ae[,1:3])
b<-fitaa%*%inv
b[1:10,1:3]

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}


totalRMSE = rep(0,30)
Varianceofdata = rep(0,30)
for(j in 1:30){
  totalRMSE[j] = totalRMSE[j] + RMSE(b[,j],mydata[,j])
  Varianceofdata[j] = var(mydata[,j])
}
R2 = totalRMSE/Varianceofdata
R2 <- R2[R2<1]
mean(R2)



RMSEtotal<-0
for( i in 1:20){
arfit<-ar(aa[1:8550,1],aic = TRUE)
b<-predict(arfit, aa[(8551+i):(8600+i),1], n.ahead = 10)$'pred'
b<-as.numeric(b)
actual<-aa[(8601+i):(8610+i),1]
RMSEtotal = RMSEtotal + sum(RMSEwithoutsum(actual,b))
}
RMSEtotal/(20*10)



#Rsquarelist<-c()
#for(ahead in 2:10){
#Rsquarelist<-c()
#for(j in 1:3){
ahead = 10
j=1
RMSEtotal<-0
library(forecast)
for( i in 1:100){
#armafit <- auto.arima(aa[(0+i):(8540+i),j])
armafit<-arima(aa[(0+i):(8540+i),j],order=c(1,0,0))
b<-forecast(armafit,h=ahead)
predictb<-as.numeric(b$mean)
actual<-aa[(8541+i):(8541+ahead-1+i),j]
RMSEtotal = RMSEtotal + RMSEwithoutsum(predictb,actual)
}
Variancetotal <- (RMSEtotal/(100))^2
Varmean <- var(aa[8541:8640,j])
Rsquare = 1 - Variancetotal/Varmean
Rsquare
#Rsquarelist<-c(Rsquarelist,Rsquare)
#}
#}








R1list<-c()
R2list<-c()
R3list<-c()
for(i in 1:9){
R1list<-c(R1list,Rsquarelist[1+3*(i-1)])
R2list<-c(R2list,Rsquarelist[2+3*(i-1)])
R3list<-c(R3list,Rsquarelist[3+3*(i-1)])
}

df<-data.frame(R1list,R2list,R3list)
colnames(df) <- c('PC1', 'PC2','PC3')
df$number<- 1:9
library(ggplot2)
ggplot(df, aes(number)) + 
  geom_line(aes(y = PC1, colour = "PC1"),size=2) + 
  geom_line(aes(y = PC2, colour = "PC2"),size=2) + 
  geom_line(aes(y = PC3, colour = "PC3"),size=2)+xlab("prediction of days ahead")+ylab("average R2 per day")+ 
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+ scale_color_discrete(name = "Method")

#################################        not average now
#######################################################################################
library(forecast)
Rsquarelist<-c()
ahead = 10
cvdays = 30
for(j in 1:3){
  RMSEtotal<-0
  RMSEt<-rep(0, 10)
  for( i in 1:cvdays){
    armafit <- auto.arima(aa[(0+i):(8540+i),j])
    b<-forecast(armafit,h=ahead)
    predictb<-as.numeric(b$mean)
    actual<-aa[(8541+i):(8541+ahead-1+i),j]
    for(t in 1:ahead){
      RMSEt[t] = RMSEt[t] + abs(actual[t]-predictb[t])
    }}
  Variancetotal <- (RMSEt/(cvdays))^2
  Varmean <- var(aa[8541:8580,j])
  Rsquare = 1 - Variancetotal/Varmean
  Rsquarelist<-c(Rsquarelist,Rsquare)
}

R1list<-Rsquarelist[1:10]
R2list<-Rsquarelist[11:20]
R3list<-Rsquarelist[21:30]
df<-data.frame(R1list,R2list,R3list)
colnames(df) <- c('PC1', 'PC2','PC3')
df$number<- 1:10
library(ggplot2)
ggplot(df, aes(number)) + 
  geom_line(aes(y = PC1, colour = "PC1"),size=2) + 
  geom_line(aes(y = PC2, colour = "PC2"),size=2) + 
  geom_line(aes(y = PC3, colour = "PC3"),size=2)+xlab("prediction of days ahead")+ylab("average R2 per day")+ 
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ 
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+ scale_color_discrete(name = "Method")



#######################################################################################


library(NHMSAR)
Rsquarelist2 <- c()
for( j in 1:5){
mydata2<-aa[1:8550,j]
T = length(mydata2)
d = 1
N.samples = 1 #dim(mydata)[3]
new.arr <- array( unlist(mydata2), dim=c(T, N.samples,d) ) 
#number of regimes
M = 4
#order of AR processes
order = 2
theta.init = init.theta.MSAR(new.arr,M=M,order=order,label="HH")
mod.hh = fit.MSAR(new.arr,theta.init,verbose=TRUE,MaxIter=30)

RMSEtotal<-0
for(timeforsimule in 1:20){
  Bsim = 100
  Ksim = Bsim*N.samples
  Y0 = array(aa[(8548+timeforsimule):(8549+timeforsimule),j],c(2,Ksim,d))
  T = 12
  Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = Ksim)
  medianYsim<-apply(Y.sim$Y[,,],1, median, na.rm = TRUE)
  medianYsim<-medianYsim[3:12]
  actual<-aa[(8551+timeforsimule):(8560+timeforsimule),j]
  RMSEtotal = RMSEtotal + sum(RMSEwithoutsum(actual,medianYsim))
}
Variancetotal <- (RMSEtotal/(20*10))^2
Varmean <- var(aa[8551:8580,j])
Rsquare = 1 - Variancetotal/Varmean
Rsquarelist2<-c(Rsquarelist2,Rsquare)
}



plot_loadings = data.frame(y = rep(c(1:6), each = 5),
                           x = rep(c(1:5), 6),
                           load_1 = ae[,1],
                           load_2 = ae[,2],
                           load_3 = ae[,3])

ggplot(plot_loadings)+
  geom_tile(aes(x,y,fill = load_1)) + xlab('lon') + ylab('lat')+
  scale_fill_continuous(limits = c(-0.35, 0.35))

ggplot(plot_loadings)+
  geom_tile(aes(x,y,fill = load_2)) + xlab('lon') + ylab('lat')+
  scale_fill_continuous(limits = c(-0.35, 0.35))

ggplot(plot_loadings)+
  geom_tile(aes(x,y,fill = load_3)) + xlab('lon') + ylab('lat')+
  scale_fill_continuous(limits = c(-0.35, 0.35))


