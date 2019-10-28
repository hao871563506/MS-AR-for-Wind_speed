
##load data
finaldf<-read.csv(file = "D:/columbia/research/californiawind/finalHMM_biggerwide.csv")
# 0 can not do boxcox, so we make it to 0.01
#finaldf[finaldf == 0] <- 0.01

##take the month 3,4,5,6
library(lubridate)
finaldf$date<-as.Date(finaldf$date, format = "%m/%d/%y")
mjdftotal<- subset(finaldf, month(date) == '3' |  month(date) == '4' |  month(date) == '5' |  month(date) == '6')

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



# extract the transition matrix for the k = num_states case
library(dplyr)

mod_trans = data.frame(mod.hh$theta$transmat) %>% dplyr::mutate(transition_from = row.names(mod.hh$theta$transmat))

mod_trans_long = reshape2::melt(mod_trans, id.vars = 'transition_from', variable.name = 'transition_to', value.name = 'prob')

library(ggplot2)
trans_matrix_plot = 
  ggplot(data = mod_trans_long, aes(y=transition_from,x=transition_to)) + 
  geom_tile(aes(fill=prob)) +
  geom_text(aes(label=sprintf("%.2f", round(prob,2)))) +
  scale_fill_gradient(limits = c(0,1),low="white",high="red") + 
  labs(x="t", y="t-1") + 
  # scale_x_continuous(breaks=seq(1,num_states)) +
  # scale_y_continuous(breaks=seq(1,num_states)) +
  theme(axis.title = element_text(size = 20)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ggtitle("Transition Matrix")

trans_matrix_plot

##simualtion
##Do a simulation to see if the model is good or not.
Bsim = 1
Ksim = Bsim*N.samples
Y0 = array(new.arr[1:2,sample(1:dim(new.arr)[2],Ksim,replace=T),],c(order,Ksim,d))
Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = Ksim)
#valid_all.MSAR(new.arr,Y.sim$Y)

###for predict, use Y0 the true 12 value list for the last day, and T = 90, for a season.
### Y0 = array(new.arr[8661:8662,1,],c(2,Ksim,d))
### T=120  #one season
### Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = Ksim)


###Box-cox back transfrom
for (i in 1:30){
  T_box = (dummy ^ lambda - 1)/lambda 
Y.sim$Y[,,i]<-(Y.sim$Y[,,i]*lambdalist[i] + 1)^(1/lambdalist[i])
}
Y.sim$Y[is.nan(Y.sim$Y)]<-0.1



###average speed for each location for each state
Y.sim_df = data.frame(V1 = Y.sim$Y[,,1])
##ii in 2:30 (12 just write 12)
for(ii in 2:30){
  Y.sim_df = cbind(Y.sim_df, Y.sim$Y[,,ii])  
}
Y.sim_df = Y.sim_df %>% setNames(c(names(mjdftotal2)[-1])) %>% dplyr::mutate(state = Y.sim$S)

Y.sim_df_long = Y.sim_df %>% reshape2::melt(, id.vars = 'state') %>%
  tidyr::separate(col = variable, into = c('lon', 'lat'), sep = '_') %>%
  tidyr::separate(col = lon, into = c('blank','lon'), sep = 'X.')

Y.sim_df_long_mean = Y.sim_df_long %>% dplyr::group_by(lat, lon, state) %>%
  dplyr::summarise(mean = mean(value)) %>%
  ungroup() %>%
  dplyr::mutate(lon = -1 * as.numeric(lon),
                lat = 1 * as.numeric(lat)) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::mutate(mean_anom = mean - mean(mean))




##map the pervious result
library(maps)
world <- data.frame(map("world", plot=FALSE)[c("x","y")])
state <- data.frame(map("state", plot=FALSE)[c("x","y")])
library(ggplot2)
#ggplot() +

#  geom_tile(data = Y.sim_df_long_mean, aes(x = (lon),y = lat,fill=(mean_anom))) + 
#  scale_fill_gradient2(name=expression(paste("Wind anom")),
#                       limits = c(-2,2),
#                       midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey") +
#  geom_path(data=world, aes(x,y,z=NULL), size = 0.25) + 
#  geom_path(data=state, aes(x,y,z=NULL), size = 0.25) + 
#  scale_y_continuous(limits = c(30,47.5)) +
#  # labs(x = "Longitude", y = "Latitude") +
#  scale_x_continuous(limits = c(-127.5,-116)) +
#  theme_bw() +
#  theme(legend.position = "bottom") +
#  facet_wrap(~state)

p1<-ggplot() +
  geom_tile(data = Y.sim_df_long_mean, aes(x = lon,y = lat,fill=(mean))) + 
  scale_fill_gradient2(name=expression(paste("Wind mean")),
                       limits = c(0,12),
                       midpoint = 6,low="blue", mid = "white",  high = "red",na.value = "grey") +
  geom_path(data=world, aes(x,y,z=NULL), size = 0.25) + 
  geom_path(data=state, aes(x,y,z=NULL), size = 0.25) + 
  scale_y_continuous(limits = c(27.5,50)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(limits = c(-127.5,-110)) +
  theme(legend.position = "bottom") +
  facet_wrap(~state)+
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))




Y.sim_df_long_mean = Y.sim_df_long %>% dplyr::group_by(lat, lon, state) %>%
  dplyr::summarise(mean = mean(value)) %>%
  ungroup() %>%
  dplyr::mutate(lon = -1 * as.numeric(lon),
                lat = 1 * as.numeric(lat)) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::mutate(std = sd(mean))

p2<-ggplot() +
  geom_tile(data = Y.sim_df_long_mean, aes(x = lon,y = lat,fill=(std))) + 
  scale_fill_gradient2(name=expression(paste("Wind standard deviation")),
                       limits = c(0,1),
                       midpoint = 0.5,low="blue", mid = "white",  high = "red",na.value = "grey") +
  geom_path(data=world, aes(x,y,z=NULL), size = 0.25) + 
  geom_path(data=state, aes(x,y,z=NULL), size = 0.25) + 
  scale_y_continuous(limits = c(27.5,50)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(limits = c(-127.5,-110)) +
  theme(legend.position = "bottom") +
  facet_wrap(~state)+
  theme_classic()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(legend.title=element_text(size=14),legend.text=element_text(size=14))



require(gridExtra)
grid.arrange(p1, p2,nrow=1, ncol=2)


state1A1<-mod.hh$theta$A$Regime1$A1
state1A2<-mod.hh$theta$A$Regime1$A2
state2A1<-mod.hh$theta$A$Regime2$A1
state2A2<-mod.hh$theta$A$Regime2$A2
state3A1<-mod.hh$theta$A$Regime3$A1
state3A2<-mod.hh$theta$A$Regime3$A2
state4A1<-mod.hh$theta$A$Regime4$A1
state4A2<-mod.hh$theta$A$Regime4$A2


state1A1<-apply(state1A1, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state1A2<-apply(state1A2, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state2A1<-apply(state2A1, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state2A2<-apply(state2A2, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state3A1<-apply(state3A1, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state3A2<-apply(state3A2, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state4A1<-apply(state4A1, 2, function(x) ifelse (abs(x) >=0.6,x,0))
state4A2<-apply(state4A2, 2, function(x) ifelse (abs(x) >=0.6,x,0))

rownames(state1A1)<-c(1:30)
colnames(state1A1)<-c(1:30)
library(reshape2)
library(ggplot2)
p1<-ggplot(melt(state1A1), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state1A1")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
rownames(state1A2)<-c(1:30)
colnames(state1A2)<-c(1:30)
rownames(state2A1)<-c(1:30)
colnames(state2A1)<-c(1:30)
rownames(state2A2)<-c(1:30)
colnames(state2A2)<-c(1:30)
rownames(state3A1)<-c(1:30)
colnames(state3A1)<-c(1:30)
rownames(state3A2)<-c(1:30)
colnames(state3A2)<-c(1:30)
rownames(state4A1)<-c(1:30)
colnames(state4A1)<-c(1:30)
rownames(state4A2)<-c(1:30)
colnames(state4A2)<-c(1:30)
p1<-ggplot(melt(state1A1), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state1A1")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p2<-ggplot(melt(state1A2), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state1A2")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p3<-ggplot(melt(state2A1), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state2A1")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p4<-ggplot(melt(state2A2), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state2A2")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p5<-ggplot(melt(state3A1), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state3A1")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p6<-ggplot(melt(state3A2), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state3A2")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p7<-ggplot(melt(state4A1), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state4A1")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")
p8<-ggplot(melt(state4A2), aes(Var1,Var2, fill=value)) + geom_raster()+  ggtitle("state4A2")+
  scale_fill_gradient2(name=expression(paste("ar value")),limits = c(-1,1),midpoint = 0,low="blue", mid = "white",  high = "red",na.value = "grey")


grid.arrange(p1, p2,p3,p4,p5,p6,p7,p8,nrow=4, ncol=2)



library(tidyverse)  
library(corrr)













#wind power
WPD<-Y.sim_df[,1:ncol(Y.sim_df)-1]
for (i in 1:nrow(WPD)){
 for (j in 1:ncol(WPD)){
   WPD[i,j] <- 0.5 * 1.225 * (WPD[i,j])^3 #1.225 is the density of air
 } 
}

#write.csv(WPD, file = "D:/columbia/research/californiawind/WPD_30_locations(4,2).csv", row.names=FALSE)





#get the viterbi
viterbi = regimes.plot.MSAR(mod.hh,new.arr)

#length of viterbi is 8662-n+2, the first 2 are NA, we ignore
viterbi<-viterbi[-c(1, 2)]

state1 <- which(viterbi %in% 1)
state2 <- which(viterbi %in% 2)
state3 <- which(viterbi %in% 3)
state4 <- which(viterbi %in% 4)


#t1=state1[1:200]
#y1=sin(t1)
#plot(t1,y1,type="l", xlab="time", ylab="Sine wave")
#t2 = state2[1:200]
#y2=sin(t2)
#lines(t2,y2,col='red')
#do a pdf plot, for 43,34,33,44. meaning the time interval of 4 and a 3, get all the time interval, plot a pdf



state1<-mydata[viterbi == 1,]

plot(state1[1:1700,1], state1[2:1701,1], main="State1 lag1",
     xlab="windspeed for t", ylab="windspeed for t+1")







sum(viterbi == 1)
sum(viterbi == 2)
sum(viterbi == 3)
sum(viterbi == 4)


sum1list<-c()
sum2list<-c()
sum3list<-c()
sum4list<-c()
sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
for(i in 0:70){
  sum1 = sum1 + sum(viterbi[(1+i*122):(40+i*122)] == 1)
  sum2 = sum2 + sum(viterbi[(1+i*122):(40+i*122)] == 2)
  sum3 = sum3 + sum(viterbi[(1+i*122):(40+i*122)] == 3)
  sum4 = sum4 + sum(viterbi[(1+i*122):(40+i*122)] == 4)
}
sum1list<-c(sum1list,sum1)
sum2list<-c(sum2list,sum2)
sum3list<-c(sum3list,sum3)
sum4list<-c(sum4list,sum4)
sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
for(i in 0:70){
  sum1 = sum1 + sum(viterbi[(41+i*122):(80+i*122)] == 1)
  sum2 = sum2 + sum(viterbi[(41+i*122):(80+i*122)] == 2)
  sum3 = sum3 + sum(viterbi[(41+i*122):(80+i*122)] == 3)
  sum4 = sum4 + sum(viterbi[(41+i*122):(80+i*122)] == 4)
}
sum1list<-c(sum1list,sum1)
sum2list<-c(sum2list,sum2)
sum3list<-c(sum3list,sum3)
sum4list<-c(sum4list,sum4)
sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
for(i in 0:70){
  sum1 = sum1 + sum(viterbi[(81+i*122):(122+i*122)] == 1)
  sum2 = sum2 + sum(viterbi[(81+i*122):(122+i*122)] == 2)
  sum3 = sum3 + sum(viterbi[(81+i*122):(122+i*122)] == 3)
  sum4 = sum4 + sum(viterbi[(81+i*122):(122+i*122)] == 4)
}
sum1list<-c(sum1list,sum1)
sum2list<-c(sum2list,sum2)
sum3list<-c(sum3list,sum3)
sum4list<-c(sum4list,sum4)


aa<-data.frame(c(1:8662),viterbi)
colnames(aa) <- c('time', 'viterbi')
write.csv(aa, file = "Viterbi.csv",row.names=FALSE)