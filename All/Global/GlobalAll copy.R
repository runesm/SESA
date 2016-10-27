#4/17/15.
library(lubridate)
library(lattice)
library(sensorgnome)
library(plyr)
library(ggplot2)
library(plotrix)
Sys.timezone()
Sys.setenv(TZ = "UTC")
setwd("C:/Users/RedHawk/Google Drive/UMaine/R/Sgnome/All/Global")
##Read Global detection files and subset important variables
Holbt2013 <- readRDS("GlobalHolbt2013.rds")
names(Holbt2013)
Holbt2013 <- Holbt2013[,c(1,2,6,17,18,20,21,23,24,25)]
names(Holbt2013)
Holbt2014 <- readRDS("GlobalHolbt2014.rds")
names(Holbt2014)
Holbt2014 <- Holbt2014[,c(1,2,6,18,19,21,22,24,25,26)]
names(Holbt2014)
Obrien2014 <- readRDS("GlobalObrien2014.rds")
names(Obrien2014)
Obrien2014 <- Obrien2014[,c(1,2,6,18,19,21,22,24,25,26)]
names(Obrien2014)
Pau2013 <- readRDS("GlobalPau2013.rds")
names(Pau2013)
Pau2013 <- Pau2013[,c(1,2,3,7,18,20,21,23,24,25)]
names(Pau2013)
Pau2014 <- readRDS("GlobalPau2014.rds")
names(Pau2014)
summary(Pau2014$proj.id)
Pau2014 <- Pau2014[,c(1,2,3,4,5,6,7,8)]
Pau2014$depYear <- "2014"
Pau2014$tagProj <- "Pau"
names(Pau2014)
df <- rbind(Holbt2013,Holbt2014,Obrien2014,Pau2013,Pau2014)
xtabs(~depYear + tagProj, data=df)
proj <- as.character(df$tagProj) #reorder to have Obrien as #3
unique(df$tagProj1)
###############################################################
#need to remove detections prior to deployment
#Import banding data, configure date in POSIXct format
Tag <- read.csv("Tag all.csv")
str(Tag)
Tag$ts1 <- with(Tag, ymd_hm(paste(Year, Month, Day, Hour, Minute),tz="America/New_York" ))
summary(Tag$ts1)
unique(Tag[Tag$Proj=="Holbt"&Tag$Year==2013,]$id)
unique(df[df$tagProj=="Holbt"&df$depYear==2013,]$id)
unique(df[!(df$id %in% Tag$id)&df$tagProj=="Holbt",]$id)
unique(df[df$id %in% c("367","365"),]$tagProj)
summary(df[df$tagProj=="Holbt" & df$id=="365",]$ts) #check ts for those tags vs 65&67
df <- df[!(df$id %in% c("367","365")),] #remove 2013 RUTU tags at PMI
Tag$ts1 <-   with_tz(Tag$ts1, "UTC") #Convert capture time to UTC
summary(Tag$ts1) 
Tag$depYear <- as.factor(Tag$Year)
Tag$tagProj<- Tag$Proj
Tag <- subset( Tag, select = -c(Year, Month, Day, Hour, Minute) )
str(Tag)
library(plyr)
print(test <- ddply(df, .(id), summarize,      #search for duplicate IDs
                    Year=length(unique(depYear)),
                    Proj=length(unique(tagProj)) )  )  
#merge the detection df with the tags df. include detections with no known BBL.ID
Global <- merge(df, Tag, by=c("id","depYear","tagProj"), all.x=TRUE)
print(test <- ddply(Global, .(id,depYear,tagProj), summarize,
                    Year=length(unique(depYear)),
                    Proj=length(unique(tagProj)),
                    BBL.ID=unique(BBL.ID) )  ) 
str(Global)
xtabs(~depYear + tagProj, data=Global)
summary(Global$id %in% Tag$id)#Check to see if there are any other ID's in Global besides those accountd for in Tags
xtabs(~depYear + tagProj, data=Global[!(Global$id %in% Tag$id),])
unique(Global[!(Global$id %in% Tag$id),]$id) #432:435 are sparrows
Global <- Global[Global$id %in% Tag$id,] #remove the foreign nanotags not found in Tag database (sparrows)
unique(Global$id)
summary(Global$id %in% Tag$id)
summary(Global$Species)
Global <- Global[Global$Species=="SESA",] #remove all non-SESA
summary(Global$Species)
Global <- droplevels(Global)
summary(Global$Species)
Pre <- subset(Global, ts < ts1) #identify detection of tags before tag deployment
unique(Pre$id)
xtabs(~Pre$tagProj+Pre$depYear)
Global <- Global[Global$ts > Global$ts1,] #remove pre-deployment detections
################################################################
library(plyr)
print(test <- ddply(Global, .(tagProj, depYear), summarize,
              ID=length(unique(id)),
              Detcs=length(id),
              First.ts = min(ts),
              Last.ts  = max(ts) )
              )
#df$date <- as.Date(df$ts)
#summary(df$date)
#Calculating DDays for All Detections
range(Global$ts)
TCalc <- ddply(Global, .(id, tagProj, depYear), summarize,
               First.ts=min(ts), Last.ts = max(ts))
str(TCalc)
TCalc$DDays <- as.numeric((TCalc$Last.ts - TCalc$First.ts), units="days")
summary(TCalc[TCalc$tagProj=="Obrien",]$DDays)
std.error(TCalc$DDays)
summary(TCalc$First.ts)
summary(TCalc$Last.ts)
summary(TCalc[TCalc$DDays > 3,]$DDays) ##DDays with those under 3 days 
#saveRDS(TCalc, "TCalc all copy.rds")
##################################################################
##############################################################################
require(rworldmap)
require(rworldxtra)
library(maps)
library(mapdata)
library(ggplot2)
# import line map
LineMap <- fortify(getMap(resolution = "high"))
## subset the NE USA from the world map
summary(Global$lat)
summary(Global$lon)
LineMap <- subset(LineMap, long > -73 & long < -61 & lat > 40 & lat < 47)
coast_map <- fortify(map("worldHires",
                         fill=TRUE,
                         plot=FALSE,
                         ylim=c(34,48),
                         xlim=c(-79,-61) ) )
## reorder proj by decreasing latitude
unique(Global$proj)
tapply(Global$lat, Global$proj, mean, na.rm=TRUE)
Global <- Global[order(Global$id, Global$ts),]#order the dataframe by id then time, for correct path plotting
##plot
Global$proj1 <- factor(Global$proj,
                       levels = unique(Global$proj[order(Global$lat,
                                                         decreasing=TRUE)] )  )
unique(Global$proj1) #plot lat ~ lon on map for all detections
p <- ggplot(data=LineMap, aes(long, lat))
p + 
  ylim(41,46)+
  xlim(-72,-62)+
  geom_map(data=coast_map, map=coast_map,
           aes(x=long, y=lat, map_id=region),
           fill="gray", alpha=0.5,colour="black") +
  geom_path(aes(group=group)) + 
  coord_map() + 
  geom_path(data=Global, aes(lon, lat,group=id),colour="blue",alpha=0.5)+ 
  #group makes sure id's are plotted separately
  geom_point(data=Global, aes(lon, lat, fill=proj1),shape=21, size=2)+
  # theme_bw() +
  theme(text = element_text(size=14))+
  ggtitle("All SESA Detections")+
  labs(x = "Longitude", y = "Latitude", fill = "Receiver Project")
#######################################################################
library(lattice)
summary(df[df$tagProj=="Pau",]$lat)
summary(df[df$tagProj=="Holbt",]$lat)
summary(df[df$tagProj=="Obrien",]$lat)
xyplot(lat ~ ts, group=tagProj, data=df[df$depYear=="2013",],
       #auto.key = list(columns = 2, cex=1.5),
       col = c("blue","magenta"),
       type = c("g","p"),
       scales=list(cex=1.5,x=list(rot=45,tick.number=10),y=list(tick.number=6)),
       xlab=list("2013 Tag Detections", cex=1.5, font=2) ,
       ylim=c(36.5,47),
       ylab=list("Latitude", cex=1.5, font=2),)
xyplot(lat ~ ts, group=tagProj, data=df[df$depYear=="2014",],
       #auto.key = list(columns = 3, cex=1.5),
       col = c("blue","red","magenta"),
       type = c("g","p"),
       scales=list(cex=1.5,x=list(rot=45,tick.number=10),y=list(tick.number=6)),
       ylim=c(36.5,47),
       xlab=list("2014 Tag Detections", cex=1.5, font=2) ,
       ylab=list("Latitude", cex=1.5, font=2),)
xyplot(lat~ lon, data=df,
       group=tagProj,
       #auto.key = list(columns = 3, cex=1.5),
       col = c("blue","red","magenta"),
       type = c("g","p"),
       aspect="iso",
       scales=list(cex=1.5),
       #pch=20,
       #cex=2,
       xlab=list("Longitude", cex=1.5, font=2) ,
       ylab=list("Latitude", cex=1.5, font=2),)
names(df)
latlong <- df[,c(5,6,7)]
#saveRDS(df, "GlobalAll.rds")
unique(Pau2014$id)
unique(Pau2013[Pau2013$lon > -66,]$id)
unique(Pau2013[Pau2013$lat > 44,]$id)
unique(Pau2014[Pau2014$lat < 42,]$id)
unique(Obrien2014$id)
unique(Obrien2014[Obrien2014$lat <40,]$id)
unique(Pau2014[Pau2014$lat < 42,]$id)

