#4/8/15. PCA based on Everitt Book.  
setwd("C:/Users/RedHawk/Google Drive/UMaine/R/SGnome/All")
Sys.setenv(TZ = "UTC")
library(lattice)
library(car)
library(lubridate)
library(plotrix)
library(xtable)
library(plyr)
Morph <- read.csv(file="MorphAll copy.csv", header = TRUE, sep = ",", quote = "\"",
                  fill = TRUE, comment.char = "" )
summary(Morph$Species)
str(Morph)
summary(is.na(Morph$BBL.ID))
summary(is.na(Morph$Age))
summary(is.na(Morph$Nanotag))
xtabs(~Proj+Year,data=Morph)
#Output<-xtabs(~Proj+Year,data=Morph)
#Create an xtable for the object that you want to export:
Table1 <- ddply(Morph, .(Proj,Year), summarize,
                Mass= round( mean(Mass, na.rm = TRUE),1 ),
                Culmen = round( mean(Culmen, na.rm = TRUE),2 ),
                EC = round( mean(EC, na.rm = TRUE),2 ) )
print(Table1)
#Table1<-xtable(Table1)
#Export your xtabled object to LaTeX:
#print.xtable(Table1, type="html", file="Mass & Culmen Table copy.html")
summary(Morph[Morph$Culmen > Morph$EC,]$Proj) #check to see that Culmen and EC look good
#xtabs(~ Proj + Year, data=Morph[Morph$Culmen > Morph$EC,])
summary(aov(Culmen ~ Proj + Year + Age, data=Morph[Morph$BBL.ID!="1151-40355",]))
Culmen.AOV <- aov(Culmen ~ Proj, data=Morph[Morph$BBL.ID!="1151-40355",])
summary(Culmen.AOV)
tapply(Morph$Culmen, Morph$Proj, mean, na.rm=TRUE)
library(agricolae)
Culmen.test <- (scheffe.test(Culmen.AOV, "Proj") )
print(Culmen.test$group)
outlierTest(lm(Culmen~Proj,data=Morph))
#Morph[379,]$Culmen <- NA   #remove outliers
#Morph[379,]
library(ggplot2)
qplot(Proj, Culmen, fill=factor(Proj), data=Morph[Morph$Year==2014,],
      geom="boxplot",
      ylab="Culmen Length (mm)") + #label y axis
  theme_bw() +
  theme(text = element_text(size=20))
p <- ggplot(data=Morph, aes(Proj, Culmen))
p + geom_boxplot(data=Morph, aes(x=Proj, y=Culmen, fill=Proj ),notch=TRUE)+ 
  facet_wrap(~ Year)+
  #scale_fill_grey()+ #convert color fill to grayscale
  scale_fill_brewer(palette="Greys")+ #better grayscale
  theme_bw()+
  #theme(text = element_text(size=24))+
  ylab("Culmen Length (mm)")+
  xlab("Project")
## examining EC values
summary(aov(EC ~ Proj + Year + Age, data=Morph))
EC.AOV <- (aov(EC ~ Proj, data=Morph))
summary(EC.AOV)
Table2<-xtable(EC.AOV)
#Export your xtabled object to LaTeX:
#print.xtable(Table2, type="html", file="table EC ~ Proj.html")
EC.test <- (scheffe.test(EC.AOV, "Proj") )
print(EC.test$group)
outlierTest(aov(EC~Proj,data=Morph))
#Morph[c(259,260,369,379),]$EC <- NA   #remove outliers
#Morph[379,]
p <- ggplot(data=Morph, aes(Proj, EC))
p + geom_boxplot(data=Morph, aes(x=Proj, y=EC, fill=Year ),notch=TRUE)+ 
  #facet_wrap(~ Year)+
  #scale_fill_grey()+ #convert color fill to grayscale
  scale_fill_brewer(palette="Greys")+ #better grayscale
  theme_bw()+
  #theme(text = element_text(size=24))+
  ylab("Exposed Culmen Length (mm)")+
  xlab("Project")
########################################################################################
#Mass and Date
summary(Morph$Date)
Morph$Date1 <- mdy(Morph$Date,tz="America/New_York")
unique(Morph$Date1)
summary(Morph$Age)
Morph[!is.na(Morph$Age) & Morph$Age=="ASY",]$Age <- "AHY"
Morph <- droplevels(Morph)
summary(Morph$Age)
var.test(Mass ~ Age, data=Morph)
t.test(Mass ~ Age, data=Morph,var.equal=TRUE)
Mass.aov <- aov(Mass ~ Julian + Proj + Year + Age, data=Morph)
summary(Mass.aov)
#Mass.aov<-xtable(aov(Mass ~ Julian + Proj + Year + Age, data=Morph))
#Export your xtabled object to LaTeX:
#print.xtable(Mass.aov, type="html", file="table Mass.aov.html")
shapiro.test(residuals(Mass.aov)) #null hypothesis is normal
outlierTest(Mass.aov) #no outliers
Mass.test <- (scheffe.test(Mass.aov, "Proj") )
print(Mass.test$group)
p <- ggplot(data=Morph, aes(Proj, Mass))
p + geom_boxplot(data=Morph, aes(x=Proj, y=Mass),fill="light gray",notch=TRUE)+ 
  #facet_wrap(~ Year)+
  theme_bw()+
  theme(text = element_text(size=16))+
  ylab("Mass (g)")+
  xlab("Project")
xyplot(Proj ~ Date1, groups=Age, data=Morph[Morph$Year==2014,],
       auto.key = list(columns = 2, cex=1.5),
       #col = c("blue","green"),
       type = c("g","p"),
       jitter.x = TRUE, jitter.y = TRUE,
       scales=list(cex=1),
       xlab=list("Capture Date", cex=1.5, font=2) ,
       ylab=list(cex=1.5, font=2))
summary(Morph$Location)
summary(aov(Mass ~ Julian+Proj+Age, data=Morph[Morph$Year==2014,]))
xyplot(Mass ~ Date1, groups=Age, data=Morph[Morph$Year==2014,],
       auto.key = list(columns = 2),
       type = c("g","p"),
       pch = 20,
       #col = c("blue","green"),
       cex = 2,
       scales=list(cex=1.5,x=list(rot=45) ),
       xlab=list("Date",cex=1.5),
       ylab=list("Mass (g)", cex=1.5) )

#######################################################################
#### Calculating Predicted Fat Mass
summary(is.na(Morph$WC))
Morph$FM <- Morph$Mass - (-9.0513 + 0.3134*Morph$WC) 
summary(lm(FM ~ Fat, data=Morph))
unique(Morph$Fat)
Morph$Fat <- as.factor(Morph$Fat)
xyplot(FM ~ Fat, group=Age, data=Morph,
       type=c("p","r"))
p <- ggplot(data=Morph, aes(Proj, FM))
p + geom_boxplot(data=Morph, aes(x=Proj, y=FM),fill="light gray",notch=TRUE)+ 
     #facet_wrap(~ Year)+
      theme_bw()+
     #theme(text = element_text(size=24))+
     theme(text = element_text(size=16))+
     ylab("Predicted Fat Mass (g)")+
     xlab("Project")
p <- ggplot(data=Morph[!is.na(Morph$Fat),], aes(Fat, FM))
p + geom_boxplot(data=Morph[!is.na(Morph$Fat),], aes(x=Fat, y=FM),fill="light gray")+  
      #facet_wrap(~ Year)+
      theme_bw()+
     #scale_fill_grey()+ #convert color fill to grayscale
     theme(text = element_text(size=16))+
     ylab("Predicted Fat Mass (g)")+
      xlab("Fat Score")
summary(Morph$FM)
outlierTest(lm(Mass~WC,data=Morph))
#Morph[,]$FM <- NA #no outliers
length(Morph[Morph$FM < 0 & !(is.na(Morph$FM)),]$FM)
Morph[Morph$FM < 0 & !(is.na(Morph$FM)),]$FM <- 0 #Convert FM values of negative numbers (n=5) to 0
var.test(FM ~ Age, data=Morph)
t.test(FM ~ Age, data=Morph,var.equal=TRUE)
densityplot(Morph$FM,group=Morph$Age)
FM.aov<-(aov(FM ~ Julian+Proj+Year+Age, data=Morph))
summary(FM.aov)
#FM.aov<-xtable(FM.aov)
#Export your xtabled object to LaTeX:
#print.xtable(FM.aov, type="html", file="table FM.aov.html")
outlierTest(FM.aov) #no outliers
shapiro.test(resid(FM.aov))#null hypothesis is normal
xyplot(FM ~ Date1 | Age,groups=Proj, data=Morph[Morph$Year==2014,],
       auto.key = list(columns = 2,cex=1.5),
       type = c("g","p"), 
       pch = 20,
       #col = c("blue","green"),
       cex = 2,
       between=list(x=0.5),
       scales=list(x=list(rot=45,tick.number=8), alternating=FALSE, cex=1.5 ),
       xlab=list("2014 Capture Date",cex=1.5),
       ylab=list("Predicted Fat Mass (g)", cex=1.5) )

#Read Morph data for subset of tagged birds
Morph$Nanotag <- as.character(Morph$Nanotag)
str(Morph$Nanotag)
unique(Morph$Nanotag)
MorphTag <- Morph[!is.na(Morph$Nanotag),]
#write.csv(MorphTag, file = "MorphTag.csv")
str(MorphTag)
var.test(Morph[is.na(Morph$Nanotag),]$Mass, MorphTag$Mass)
t.test(Morph[is.na(Morph$Nanotag),]$Mass, MorphTag$Mass,var.equal=FALSE)
var.test(Morph[is.na(Morph$Nanotag),]$FM, MorphTag$FM)
t.test(Morph[is.na(Morph$Nanotag),]$FM, MorphTag$FM, var.equal=FALSE)

#Read telemetry detection period data
DDays <- readRDS(file="TCalc all copy.rds")
str(DDays)
DDays$Nanotag <- DDays$id
DDays$Proj <- DDays$tagProj
DDays$Year <- DDays$depYear
str(DDays)
unique(DDays$Proj)
DDays[DDays$Proj == "Obrien",]$Proj <- "OBrien"
DDays <- DDays[,6:9]
MorphTag <- merge(MorphTag, DDays, by=c("Year","Proj","Nanotag"), all.x=TRUE)
summary(MorphTag$DDays) #NAs are tags that were deployed but not detected
MorphTag$DDays3 <- NA
MorphTag[!is.na(MorphTag$DDays) & MorphTag$DDays > 3,]$DDays3 <- MorphTag[!is.na(MorphTag$DDays) & MorphTag$DDays > 3,]$DDays
summary(MorphTag$DDays3)
ddply(MorphTag, .(Proj), summarize,
     DDays3= round( mean(DDays3, na.rm = TRUE),2 ) )
std.error(MorphTag[MorphTag$Age =="AHY",]$DDays3, na.rm=TRUE)
std.error(MorphTag[MorphTag$Age =="HY",]$DDays3, na.rm=TRUE)
var.test(DDays3 ~ Age, data=MorphTag)
t.test(DDays3 ~ Age, data=MorphTag, var.equal=TRUE)
DDays3.aov <- aov(DDays3 ~ Julian + Proj + Year + Age, data=MorphTag)
summary(DDays3.aov)
summary(aov(DDays3 ~ Julian*Proj+Age, data=MorphTag))
#DDays3.aov<-xtable(aov(DDays3 ~ Julian + Proj + Year + Age, data=MorphTag))
#Export your xtabled object to LaTeX:
#print.xtable(DDays3.aov, type="html", file="table DDays3.aov.html")
shapiro.test(residuals(DDays3.aov)) #null hypothesis is normal
outlierTest(DDays3.aov) #no outliers
p <- ggplot(data=MorphTag, aes(Proj, DDays3))
p + geom_boxplot(data=MorphTag, aes(x=Proj, y=DDays3),fill="light gray",notch=FALSE)+ 
  geom_point(data=MorphTag, aes(x=Proj, y=DDays3), alpha=0.2,size=2)+
  #facet_wrap(~ Year)+
  theme_bw()+
  theme(text = element_text(size=16))+
  ylab("Detection Period (Days")+
  xlab("Project")
p <- ggplot(data=MorphTag, aes(Age, DDays3))
p + geom_boxplot(data=MorphTag, aes(x=Age, y=DDays3),fill="light gray",notch=FALSE)+ 
  facet_wrap(~ Proj)+
  theme_bw()+
  theme(text = element_text(size=16))+
  ylab("Detection Period (Days")+
  xlab("Age Class")
p <- ggplot(data=MorphTag, aes(DDays3, fill = Proj, colour = Proj)) 
p +  geom_density(alpha = 0.2)+
  xlim(0, 35)
p <- ggplot(data=MorphTag, aes(Julian, DDays3))
p + geom_point(data=MorphTag, aes(x=Julian, y=DDays3,colour=Proj))+ 
  geom_smooth()+ 
 # facet_wrap(~ Proj)+
  theme_bw()+
  theme(text = element_text(size=16))+
  ylab("Detection Period (Days")+
  xlab("Ordinal Date")
library(splines)
library(MASS)
summary(lm(DDays3 ~ FM, data=MorphTag[MorphTag$Proj=="Holbt",]))
shapiro.test(resid(lm(DDays3 ~ FM, data=MorphTag[MorphTag$Proj=="Holbt",]) )  ) #null hypothesis is normal
p <- ggplot(data=MorphTag[MorphTag$Proj=="Holbt",], aes(FM, DDays3))
p + 
  geom_point(data=MorphTag[MorphTag$Proj=="Holbt",], aes(x=FM, y=DDays3))+ 
  geom_smooth(method = "lm", formula = y ~ x,level=0.95,col="black")+ 
  theme_bw()+
  theme(text = element_text(size=16))+
  ylab("Detection Period (Days)")+
  xlab("Predicted Fat Mass (g)")
xyplot(DDays3 ~ FM|Proj, group= Age, data=MorphTag,
       auto.key = list(columns = 2,cex=1.5),
       type = c("g","p"), 
       pch = 20,
       #col = c("blue","green"),
       cex = 2,
       scales=list(x=list(tick.number=8), cex=1.5 ),
       xlab=list("Predicted Fat Mass (g)",cex=1.5,font=2),
       ylab=list("Detection Period (days)", cex=1.5,font=2) )
#########################################################################
#Print Table and Validate DDays 
DDays.table <- ddply(MorphTag, .(Proj,Year), summarize,
                Mean= round( mean(DDays, na.rm = TRUE),1 ),
                Median= round( median(DDays, na.rm = TRUE),1 ),
                Min= round( min(DDays, na.rm = TRUE),1 ),
                Max= round( max(DDays, na.rm = TRUE),1 ) )
print(DDays.table)
#DDays.table<-xtable(DDays.table)
#Export your xtabled object to LaTeX:
#print.xtable(DDays.table, type="html", file="DDays Table.html")
######################################################################
# Use plyr and reshape2 to create summary statistics for each group
library(reshape2)
library(plyr)
str(MorphTag)
Sub <- MorphTag[,c(2,8,15,17,19)]
Melted <- melt(Sub, id.vars=c("Proj")) #reshape df so each response variable is its own row
test <- ddply(Melted, c("Proj","variable"), summarise,
              mean = mean(value), sd = sd(value), n=length(value),
              sem = sd(value)/sqrt(length(value)))
std.error(Morph[Morph$Age =="HY",]$FM, na.rm=TRUE)
summary(Morph[Morph$Age =="HY",]$FM)
#########################