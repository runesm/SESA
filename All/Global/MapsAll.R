#4/17/15.
setwd("C:/Users/RedHawk/Google Drive/UMaine/R/Sgnome/All/Global")
require(ggplot2)
require(lubridate)
require(plyr)
require(rworldmap)
require(rworldxtra)
require(mapproj)

## get the tower data
all.rc.2013 <- read.csv("Rec2013.csv") ## add in the appropriate directory structure here
#readRDS(df, "GlobalAll.rds")
# example plot
newmap <- getMap(resolution = "high")  # different resolutions available

newmap.df <- fortify(newmap)

## subset the Gulf of Maine from the world map
summary(df$lat)
summary(df$lon)
GOM.shoreline.df <- subset(newmap.df, long > -78 & long < -59 & lat > 35 & lat < 48)

## plot it
unique(df$tagProj)
p <- ggplot(data=GOM.shoreline.df, aes(long, lat))
p + geom_path(aes(group=group)) + 
  coord_map() + 
  geom_point(data=all.rc.2013, aes(long, lat),shape=24, col="black",fill="white", size=5)+ 
  #geom_point(data=df[df$tagProj=="Holbt",], aes(lon, lat), col="blue", size=5)+
  #geom_point(data=df[df$tagProj=="Obrien",], aes(lon, lat), col="red", size=5)+
  #geom_point(data=df[df$tagProj=="Pau",], aes(lon, lat), col="magenta",  size=5)+
  #theme_bw() +
  theme(text = element_text(size=24))+
  ylab("Latitude")+
  xlab("Longitude")


## map using gmap
mymap <- gmap("NB, Canada", zoom=6, scale=2, type="hybrid")
plot(mymap)

select.area <- drawExtent()
# now click 2 times on the map to select your region
mymap <- gmap(select.area)
plot(mymap)
RC2013 <- all.rc.2013[,(3:4)]
coordinates(RC2013) <- c("long", "lat")  # set spatial coordinates
plot(RC2013)

Merc.RC2013 <- Mercator(RC2013)  # Google Maps are in Mercator projection. 
# This function projects the points to that projection to enable mapping
plot(RC2013, pch = 20, col = "blue")
points(Merc.RC2013, pch = 20, col = "green")



