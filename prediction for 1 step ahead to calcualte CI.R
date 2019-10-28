
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
n<-1
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


##iterated algorithm for n step

predictiontotal<- data.frame(matrix(ncol = 30, nrow = 0))
#number of regimes
M = 4
#order of AR processes
order = 2
theta.init = init.theta.MSAR(new.arr,M=M,order=order,label="HH")
mod.hh = fit.MSAR(new.arr,theta.init,verbose=TRUE,MaxIter=50)


y.p.2 = prediction.MSAR(new.arr,mod.hh$theta,ex=1:N.samples)


predictiony<-y.p.2$y.p
predictionvar<-y.p.2$var.p














sum = 0

for( i in 11:8640){
  for(location in 1:30){
  lowerbound = predictiony[i,1,location]-2*sqrt((predictionvar[i,1,location,location]))
  upperbound = predictiony[i,1,location]+2*sqrt((predictionvar[i,1,location,location]))
  #upperbound = predictiony[i,1,location]+2*sqrt(predictionvar[i,1,location,location])
  real = mydata[i+1,location]
  if (real>=upperbound||real<=lowerbound){ sum = sum + 1}
  }}

100-sum/(8630*30)*100




variancetotal<-sapply(mydata, sd)
averagewind<-colMeans(mydata)

sum = 0

for( i in 11:8640){
  for(location in 1:30){
    lowerbound = mydata[i,location]-2*variancetotal[location]
    upperbound = mydata[i,location]+2*variancetotal[location]
    real = averagewind[location]
    if (real>=upperbound||real<=lowerbound){ sum = sum + 1}
  }}

100-sum/(8630*30)*100





