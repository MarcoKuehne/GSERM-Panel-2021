#########################################################
# GSERM - St. Gallen (2021)
#
# Analyzing Panel Data
# Prof. Christopher Zorn
#
# Day Two: "Unit Effects" models.
#
########################################################
# Load packages (install as needed, using
# install.packages()), and set some options:

library(RCurl)
library(haven)
library(psych)
library(plyr)
library(lme4)
library(plm)
library(gtools)
library(boot)
library(plyr)
library(dplyr)
library(texreg)
library(statmod)
library(pscl)
library(stargazer)
library(prais)
library(nlme)
# install.packages("tseries")
library(tseries)
# install.packages("pcse")
library(pcse)
# install.packages("panelView")
library(panelView)
library(performance)

options(scipen = 6) # bias against scientific notation
options(digits = 4) # show fewer decimal places

# setwd("~/Dropbox (Personal)/GSERM/Panel-2021/Notes and Slides")

########################################
# Cross-country example data, 1945-2014:

TSCSURL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-Panel-2021/master/Data/DemonstrationsTSCS.csv"
temp<-getURL(TSCSURL)
Demos<-read.csv(textConnection(temp))
rm(temp)

Demos$X<-NULL
Demos$ColdWar <- with(Demos, 
                      ifelse(Year<1990,1,0))

# Create variables; kill duplicates...

Demos$lnGDP <- log(Demos$GDP)
Demos <- Demos[!duplicated(Demos[c("Year", "ccode")]),] # kill duplicates
PDF <- pdata.frame(Demos,index=c("ccode","Year")) # panel data frame

summary(Demos)

#################################
# Visualizing panel data...

pdf("PanelDemosViz.pdf",7,5)
panelView(lnDemons~POLITY+Monarch,data=Demos,theme.bw=TRUE,
          outcome.type="continuous",type="outcome",
          by.timing=TRUE,index=c("ccode","Year"),
          main=" ",ylab="log(Demonstrations)",
          legendOff=TRUE,)
dev.off()

pdf("PanelMonarchViz.pdf",7,5)
panelView(lnDemons~Monarch,data=Demos,theme.bw=TRUE,
          by.timing=FALSE,index=c("ccode","Year"),
          color=c("white","darkgreen","red"),
          legend.labs=c("Missing","Not Monarchy","Monarchy"),
          main=" ",ylab="Country Code",axis.lab.gap=c(5,5),
          background="white")
dev.off()

#########################################
# Variation: Total, within, and between 

with(Demos, describe(lnDemons)) # all variation

pdf("DemonsAll.pdf",6,5)
par(mar=c(4,4,2,2))
plot(density(Demos$lnDemons,na.rm=TRUE),
     main="",xlab="ln(Demonstrations)",
     lwd=2)
dev.off()

DemonsMeans <- ddply(Demos,.(ccode),summarise,
                     DemonsMean = mean(lnDemons,na.rm=TRUE))

with(DemonsMeans, describe(DemonsMean)) # "between" variation

pdf("DemonsBetween.pdf",6,5)
par(mar=c(4,4,2,2))
plot(density(DemonsMeans$DemonsMean,na.rm=TRUE),
     main="",xlab="Mean ln(Demonstrations)",
     lwd=2)
dev.off()

Demos <- ddply(Demos, .(ccode), mutate,
               DemonsMean = mean(lnDemons,na.rm=TRUE))
Demos$DemonsWithin <- with(Demos, lnDemons-DemonsMean)

with(Demos, describe(DemonsWithin))

pdf("DemonsWithin.pdf",6,5)
par(mar=c(4,4,2,2))
plot(density(Demos$DemonsWithin,na.rm=TRUE),
     main="",xlab="Demonstrations: Within-Country Variation",
     lwd=2)
abline(v=0,lty=2)
dev.off()

#####################################################
# One-Way Unit Effects models:

# Pooled OLS:

OLS<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
        data=PDF,model="pooling")

summary(OLS)

# "Fixed" / within effects:

FE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
        data=PDF, effect="individual",model="within")

summary(FE)

# Make a table:

Table1 <- stargazer(OLS,FE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="UFX1.tex")

# Time-period Fixed Effects:

FE.Time<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
        data=PDF,effect="time",model="within")

# summary(FET)

# A comparison table:

FE.Units <- FE
CompFETable <- stargazer(FE.Units,FE.Time,
                    title="FE Models of Demonstrations",
                    column.separate=c(1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="FEComp.tex")

# Test for \alpha_i = 0:

pFtest(FE,OLS)
plmtest(FE,effect=c("individual"),type=c("bp"))
plmtest(FE,effect=c("individual"),type=c("kw"))

pFtest(FE.Time,OLS)
plmtest(FE.Time,effect=c("time"),type=c("bp"))
plmtest(FE.Time,effect=c("time"),type=c("kw"))

#####################
# Interpretation...

with(Demos, sd(Monarch,na.rm=TRUE)) # all variation

Demos <- ddply(Demos, .(ccode), mutate,
               MonMean = mean(Monarch,na.rm=TRUE))
Demos$MonarchWithin <- with(Demos, Monarch-MonMean)

with(Demos, sd(MonarchWithin,na.rm=TRUE)) # "within" variation

###################
# Between effects:

BE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
        data=PDF, effect="individual",model="between")

summary(BE)

Table2 <- stargazer(OLS,FE,BE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="UFX2.tex")

###################
# Random effects:

RE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar,
        data=PDF, effect="individual", model="random")

summary(RE)

Table3 <- stargazer(OLS,FE,BE,RE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="UFX3.tex")

# Hausman test:

phtest(FE, RE)  # ugh...

#######################
# HLM bit...

AltRE<-lmer(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar+
              (1|ccode), data=Demos)

summary(AltRE)

# Are they the same?

TableHLM <- stargazer(RE,AltRE,
                    title="RE Models of Demonstrations",
                    column.separate=c(1,1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="HLMTable.tex")

# ... yes. Yes they are.
#
# More HLM fun...

HLM1<-lmer(lnDemons~POLITY+I(POLITY^2)+lnGDP+(lnGDP|ccode)+
             Monarch+ColdWar, data=Demos,
             control=lmerControl(optimizer="bobyqa"))

summary(HLM1)

# Testing:

anova(AltRE,HLM1)
VarCorr(HLM1)

# Get some of those sweeeet random slopes:

Bs<-data.frame(coef(HLM1)[1])

head(Bs)
mean(Bs$ccode..Intercept.)
mean(Bs$ccode.lnGDP)


pdf("Demo-RandomIntercepts.pdf",6,5)
par(mar=c(4,4,2,2))
with(Bs, plot(density(ccode..Intercept.),lwd=3,
              main="",xlab="Intercept Values"))
abline(v=mean(Bs$ccode..Intercept.),lty=2)
dev.off()

pdf("Demo-GDPRandomSlopes.pdf",6,5)
par(mar=c(4,4,2,2))
with(Bs, plot(density(ccode.lnGDP),lwd=3,
              main="",xlab="GDP Slopes"))
abline(v=mean(Bs$ccode.lnGDP),lty=2)
dev.off()

REcorr<-with(Bs,cor(ccode..Intercept.,ccode.lnGDP))

pdf("Demo-HLMScatter.pdf",6,5)
par(mar=c(4,4,2,2))
with(Bs, plot(ccode..Intercept.,ccode.lnGDP,
              pch=19,main="",xlab="Intercept Values",
              ylab="GDP Slopes"))
text(3,0.4,cex=1.2,
     labels=paste0("r = ",round(REcorr,3)))
dev.off()

###############################################
# Separating "within" and "between" effects:
# Wealth and demonstrations:

Demos<-ddply(Demos,.(ccode),mutate,
             GDP.Between=mean(lnGDP,na.rm=TRUE))
Demos$GDP.Within<- (Demos$lnGDP - Demos$GDP.Between) 

WEBE.OLS<-lm(lnDemons~POLITY+I(POLITY^2)+GDP.Within+GDP.Between+
           Monarch+ColdWar,data=Demos)

summary(WEBE.OLS)

Table4 <- stargazer(WEBE.OLS,
                    title="BE + WE Model of Demonstrations",
                    column.separate=c(1,1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "Within-Country ln(GDP)",
                                       "Between-Country ln(GDP)",
                                       "Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="WEBE.tex")

#######################################################
# Two-way effects:

TwoWayFE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar,
        data=PDF,effect="twoway",model="within")

summary(TwoWayFE)

# Testing...
#
# Two-way effects:

pFtest(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar,
       data=PDF,effect="twoway",model="within")
plmtest(TwoWayFE,c("twoways"),type=("kw"))

# One-way effects in the two-way model:

plmtest(TwoWayFE,c("individual"),type=("kw"))
plmtest(TwoWayFE,c("time"),type=("kw"))

# Equivalence to -lm-:

TwoWayFE.BF<-lm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+
                   factor(ccode)+factor(Year),data=PDF)

summary(TwoWayFE.BF)

# Two-way random effects:

TwoWayRE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar,
              data=PDF,effect="twoway",model="random")

summary(TwoWayRE)

# Here's a nicer table:

Table4 <- stargazer(OLS,FE,BE,RE,TwoWayFE,TwoWayRE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1,1,1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    column.sep.width="1pt",
                    omit.stat=c("f"),out="UFX4.tex")
