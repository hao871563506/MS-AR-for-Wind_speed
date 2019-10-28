library(ncdf4)
ncpath <- "D:/columbia/research/nywinddata/"
ncname <- "uwind"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
ncin <- nc_open(ncfname)

#u10
v1 <- ncin$var[[1]]
#data1[lon,lat,time]
data1 <- ncvar_get(ncin,v1)

#lon 576 88:108
lontable<-v1$dim[[1]]$vals
#lat 361 245:267
lattable<-v1$dim[[2]]$vals
#time "minutes since 1980-01-01 00:30:00"
timetable<-v1$dim[[4]]$vals

nc_close(ncin)



u10df <- data.frame(lon=character(),
                 lat=character(), 
                 date=character(), 
                 windspeed=character(),
                 stringsAsFactors=FALSE) 

for (i in 1:length(lontable)){
  for (j in 1:length(lattable)){
    for (k in 1:length(timetable)){
      clon=lontable[i]
      clat=lattable[j]
      ctime=timetable[k]
      cwindspeed=data1[i,j,k]
      u10df[nrow(u10df) + 1,] = list(clon,clat,ctime,cwindspeed)
    }}}

u10df$lon <- as.numeric(as.character(u10df$lon))
u10df$lat <- as.numeric(as.character(u10df$lat))
u10df$windspeed <- as.numeric(as.character(u10df$windspeed))

#write.csv(u10df, file = "D:/columbia/research/nywinddata/uwind.csv")



gpclibPermit()
library(maps)
library(maptools)
# index the cells that are/are not over CA
land = map("state",region = "ca", plot = FALSE, fill = TRUE)
land_IDs <- sapply(strsplit(land$names, ":"), function(x) x[1])
land <- map2SpatialPolygons(land, IDs = land_IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))
land_pts <- SpatialPoints(u10df[,c("lon","lat")],proj4string = CRS(proj4string(land))) # JUST NEED TO CHANGE THIS LINE
land_ii = !is.na(over(land_pts,land))
u10df_CA = u10df[land_ii,]

#take a backup
#write.csv(u10df_CA, file = "D:/columbia/research/californiawind/u10df_CA.csv")
#u10df_CA<-read.csv(file = "D:/columbia/research/californiawind/u10df_CA.csv")
#v10df_CA<-read.csv(file = "D:/columbia/research/californiawind/v10df_CA.csv")
#finaldf<-u10df_CA
#for(i in 1:nrow(u10df_CA)){
#  finaldf[i,]$windspeed = sqrt(u10df_CA[i,]$windspeed^2+v10df_CA[i,]$windspeed^2)
#}
#write.csv(finaldf, file = "D:/columbia/research/californiawind/finaldf.csv")
#finaldf<-read.csv(file = "D:/columbia/research/californiawind/finaldf.csv")

#59052/444, there is 133 points, 444/12month = 37years.
Jantable <- data.frame(lon=numeric(),
                    lat=numeric(), 
                    date=numeric(), 
                    windspeed=numeric(),
                    stringsAsFactors=FALSE) 

for(i in seq(1, nrow(finaldf), 12)){
  Jantable[nrow(Jantable) + 1,]=finaldf[i,]
}
Janmean<-c()
Janvar<-c()
for(i in seq(1, nrow(Jantable)-37, 37)){
  j=37+i
  Janmean <- c(Janmean,abs(mean(Jantable[c(i:j),]$windspeed)))
  Janvar <- c(Janvar,var(Jantable[c(i:j),]$windspeed))
}

library(ggplot2)
require(gridExtra)
xaxis<-c(1:132)
Jandf<-data.frame(Janmean,Janvar,xaxis)
Jan1<-ggplot(Jandf, aes(x=xaxis, y=Janmean)) + geom_point(size=2, shape=23)+labs(title="Jan Mean",x ="points", y = "Mean")
Jan2<-ggplot(Jandf, aes(x=xaxis, y=Janvar)) + geom_point(size=2, shape=23)+labs(title="Jan Var",x ="points", y = "Var")
grid.arrange(Jan1, Jan2, nrow = 1)

Febtable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Martable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Aprtable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Maytable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Juntable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Jultable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Augtable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Septable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Octtable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Novtable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 
Dectable <- data.frame(lon=numeric(),lat=numeric(), date=numeric(), windspeed=numeric(),stringsAsFactors=FALSE) 

for(i in seq(2, nrow(u10df_CA), 12)){Febtable[nrow(Febtable) + 1,]=u10df_CA[i,]}
Febmean<-c();Febvar<-c()
for(i in seq(1, nrow(Febtable)-37, 37)){j=37+i
  Febmean <- c(Febmean,abs(mean(Febtable[c(i:j),]$windspeed)))
  Febvar <- c(Febvar,var(Febtable[c(i:j),]$windspeed))}

for(i in seq(3, nrow(u10df_CA), 12)){Martable[nrow(Martable) + 1,]=u10df_CA[i,]}
Marmean<-c();Marvar<-c()
for(i in seq(1, nrow(Martable)-37, 37)){j=37+i
  Marmean <- c(Marmean,abs(mean(Martable[c(i:j),]$windspeed)))
  Marvar <- c(Marvar,var(Martable[c(i:j),]$windspeed))}

for(i in seq(4, nrow(u10df_CA), 12)){Aprtable[nrow(Aprtable) + 1,]=u10df_CA[i,]}
Aprmean<-c();Aprvar<-c()
for(i in seq(1, nrow(Aprtable)-37, 37)){j=37+i
Aprmean <- c(Aprmean,abs(mean(Aprtable[c(i:j),]$windspeed)))
Aprvar <- c(Aprvar,var(Aprtable[c(i:j),]$windspeed))}

for(i in seq(5, nrow(u10df_CA), 12)){Maytable[nrow(Maytable) + 1,]=u10df_CA[i,]}
Maymean<-c();Mayvar<-c()
for(i in seq(1, nrow(Maytable)-37, 37)){j=37+i
Maymean <- c(Maymean,abs(mean(Maytable[c(i:j),]$windspeed)))
Mayvar <- c(Mayvar,var(Maytable[c(i:j),]$windspeed))}

for(i in seq(6, nrow(u10df_CA), 12)){Juntable[nrow(Juntable) + 1,]=u10df_CA[i,]}
Junmean<-c();Junvar<-c()
for(i in seq(1, nrow(Juntable)-37, 37)){j=37+i
Junmean <- c(Junmean,abs(mean(Juntable[c(i:j),]$windspeed)))
Junvar <- c(Junvar,var(Juntable[c(i:j),]$windspeed))}

for(i in seq(7, nrow(u10df_CA), 12)){Jultable[nrow(Jultable) + 1,]=u10df_CA[i,]}
Julmean<-c();Julvar<-c()
for(i in seq(1, nrow(Jultable)-37, 37)){j=37+i
Julmean <- c(Julmean,abs(mean(Jultable[c(i:j),]$windspeed)))
Julvar <- c(Julvar,var(Jultable[c(i:j),]$windspeed))}

for(i in seq(8, nrow(u10df_CA), 12)){Augtable[nrow(Augtable) + 1,]=u10df_CA[i,]}
Augmean<-c();Augvar<-c()
for(i in seq(1, nrow(Augtable)-37, 37)){j=37+i
Augmean <- c(Augmean,abs(mean(Augtable[c(i:j),]$windspeed)))
Augvar <- c(Augvar,var(Augtable[c(i:j),]$windspeed))}

for(i in seq(9, nrow(u10df_CA), 12)){Septable[nrow(Septable) + 1,]=u10df_CA[i,]}
Sepmean<-c();Sepvar<-c()
for(i in seq(1, nrow(Septable)-37, 37)){j=37+i
Sepmean <- c(Sepmean,abs(mean(Septable[c(i:j),]$windspeed)))
Sepvar <- c(Sepvar,var(Septable[c(i:j),]$windspeed))}

for(i in seq(10, nrow(u10df_CA), 12)){Octtable[nrow(Octtable) + 1,]=u10df_CA[i,]}
Octmean<-c();Octvar<-c()
for(i in seq(1, nrow(Octtable)-37, 37)){j=37+i
Octmean <- c(Octmean,abs(mean(Octtable[c(i:j),]$windspeed)))
Octvar <- c(Octvar,var(Octtable[c(i:j),]$windspeed))}

for(i in seq(11, nrow(u10df_CA), 12)){Novtable[nrow(Novtable) + 1,]=u10df_CA[i,]}
Novmean<-c();Novvar<-c()
for(i in seq(1, nrow(Novtable)-37, 37)){j=37+i
Novmean <- c(Novmean,abs(mean(Novtable[c(i:j),]$windspeed)))
Novvar <- c(Novvar,var(Novtable[c(i:j),]$windspeed))}

for(i in seq(12, nrow(u10df_CA), 12)){Dectable[nrow(Dectable) + 1,]=u10df_CA[i,]}
Decmean<-c();Decvar<-c()
for(i in seq(1, nrow(Dectable)-37, 37)){j=37+i
Decmean <- c(Decmean,abs(mean(Dectable[c(i:j),]$windspeed)))
Decvar <- c(Decvar,var(Dectable[c(i:j),]$windspeed))}

#Take max for each point
Maxforlocation<-c()
for(i in seq(1,132,1)){
  Maxforlocation<-c(Maxforlocation,max(Janmean[i],Febmean[i],Marmean[i],Aprmean[i],Maymean[i],Junmean[i],Julmean[i],Augmean[i],Sepmean[i],Octmean[i],Novmean[i],Decmean[i]))
}

#define a ratio of the maximum
ratio=0.4
#define a threshold for variance
CV=0.3
Janlocation<-c();Feblocation<-c();Marlocation<-c();Aprlocation<-c();Maylocation<-c();Junlocation<-c();Jullocation<-c();Auglocation<-c();Seplocation<-c();Octlocation<-c();Novlocation<-c();Declocation<-c();
for(i in seq(1,132,1)){
  if(Janmean[i]>ratio*Maxforlocation[i]){
    if(Janvar[i]>CV){
      Janlocation<-c(Janlocation,i)}}}
for(i in seq(1,132,1)){
  if(Febmean[i]>ratio*Maxforlocation[i]){
    if(Febvar[i]>CV){
      Feblocation<-c(Feblocation,i)}}}
for(i in seq(1,132,1)){
  if(Marmean[i]>ratio*Maxforlocation[i]){
    if(Marvar[i]>CV){
      Marlocation<-c(Marlocation,i)}}}
for(i in seq(1,132,1)){
  if(Aprmean[i]>ratio*Maxforlocation[i]){
    if(Aprvar[i]>CV){
      Aprlocation<-c(Aprlocation,i)}}}
for(i in seq(1,132,1)){
  if(Maymean[i]>ratio*Maxforlocation[i]){
    if(Mayvar[i]>CV){
      Maylocation<-c(Maylocation,i)}}}
for(i in seq(1,132,1)){
  if(Junmean[i]>ratio*Maxforlocation[i]){
    if(Junvar[i]>CV){
      Junlocation<-c(Junlocation,i)}}}
for(i in seq(1,132,1)){
  if(Julmean[i]>ratio*Maxforlocation[i]){
    if(Julvar[i]>CV){
      Jullocation<-c(Jullocation,i)}}}
for(i in seq(1,132,1)){
  if(Augmean[i]>ratio*Maxforlocation[i]){
    if(Augvar[i]>CV){
      Auglocation<-c(Auglocation,i)}}}
for(i in seq(1,132,1)){
  if(Sepmean[i]>ratio*Maxforlocation[i]){
    if(Sepvar[i]>CV){
      Seplocation<-c(Seplocation,i)}}}
for(i in seq(1,132,1)){
  if(Octmean[i]>ratio*Maxforlocation[i]){
    if(Octvar[i]>CV){
      Octlocation<-c(Octlocation,i)}}}
for(i in seq(1,132,1)){
  if(Novmean[i]>ratio*Maxforlocation[i]){
    if(Novvar[i]>CV){
      Novlocation<-c(Novlocation,i)}}}
for(i in seq(1,132,1)){
  if(Decmean[i]>ratio*Maxforlocation[i]){
    if(Decvar[i]>CV){
      Declocation<-c(Declocation,i)}}}
