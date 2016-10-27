#June 15 2016, adding variables to master list
setwd("C:/Users/ASUS/R/Sgnome/All")
Sys.setenv(TZ = "UTC")
library(lattice)
library(car)
library(lubridate)
library(plotrix)
DDays <- read.csv(file="DDaysAll.csv", header = TRUE, sep = ",", quote = "\"",
                  fill = TRUE, comment.char = "", )
#FM <- read.csv(file="FMAllcopy.csv", header = TRUE, sep = ",", quote = "\"",
#               fill = TRUE, comment.char = "", )
Master <- read.csv(file="Mastercopy.csv", header = TRUE, sep = ",", quote = "\"",
                   fill = TRUE, comment.char = "", )
summary(is.na(Master$BBL.ID))
summary(is.na(Master$Mass))
#Master1 <- merge(Master, FM, by = "BBL.ID", all.x = TRUE)
#Master2 <- merge(Master1, DDays, by = "BBL.ID", all.x = TRUE)
########################################################################
Master$FM <- NA
summary(Master$Species)
summary(Master$Mass)
summary(Master$WC)
summary(Master[Master$Species!="SESA",]$Species)
Master$FM <- Master$Mass - (-9.0513 + 0.3134*Master$WC) #Calculate estimate of Fat Mass (for all captures)
summary(Master$FM)
length(Master[Master$Species!="SESA"&!is.na(Master$FM),]$FM) #identify FM values for non-SESA captures that are not NA
Master[Master$Species!="SESA"&!is.na(Master$FM),]$FM <- NA #make non-SESA FM values NA
summary(Master$FM)
length(Master[Master$FM < 0 & !(is.na(Master$FM)),]$FM)
Master[Master$FM < 0 & !(is.na(Master$FM)),]$FM <- 0 #Convert FM values of negative numbers (n=4) to 0
summary(Master$FM)
#write.csv(Master, file = "MasterOutput.csv")
