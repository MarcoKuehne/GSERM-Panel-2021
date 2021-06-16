#########################################################
# GSERM - St. Gallen (2021)
#
# Analyzing Panel Data
# Prof. Christopher Zorn
#
# Day Three: Panel Dynamics.
#
########################################################
# Load packages (install as needed, using
# install.packages()), and set some options:

library(RCurl)
library(haven)
library(psych)
library(gtools)
library(boot)
library(plyr)
library(dplyr)
library(texreg)
library(statmod)
library(pscl)
library(lme4)
library(plm)
library(stargazer)
library(prais)
library(nlme)
library(panelAR)
# install.packages("tseries")
library(tseries)
# install.packages("pcse")
library(pcse)
# install.packages("pgmm")
library(pgmm)
# install.packages("dynpanel")
library(dynpanel)
# install.packages("panelView")
library(panelView)
library(dotwhisker)
library(performance)

options(scipen = 6) # bias against scientific notation
options(digits = 4) # show fewer decimal places

setwd("~/Dropbox (Personal)/GSERM/Panel-2021/Notes and Slides")
###########################################
# GLS-ARMA...
#
# Cross-country example data on 
# political demonstrations, 1945-2014:

TSCSURL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-Panel-2021/master/Data/DemonstrationsTSCS.csv"
temp<-getURL(TSCSURL)
Demos<-read.csv(textConnection(temp))
rm(temp,TSCSURL)

Demos$X<-NULL
Demos$ColdWar <- with(Demos, 
                      ifelse(Year<1990,1,0))

# Create variables; kill duplicates...

Demos$lnGDP <- log(Demos$GDP)
Demos <- Demos[!duplicated(Demos[c("Year", "ccode")]),] # kill duplicates
PDF <- pdata.frame(Demos,index=c("ccode","Year")) # panel data frame

summary(Demos)
###########################################
# Models:

# Pooled OLS:

OLS<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
         data=PDF,model="pooling")

summary(OLS)

# Prais-Winsten:

PraisWinsten<-panelAR(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
                     data=Demos,panelVar="ccode",timeVar="Year",
                     autoCorr="ar1",panelCorrMethod="none",
                     rho.na.rm=TRUE)

summary(PraisWinsten)
PraisWinsten$panelStructure$rho

# GLS with homoscedasticity & AR(1) autocorrelation:

GLS<- gls(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar,
              Demos,correlation=corAR1(form=~1|ccode),
              na.action="na.omit")
summary(GLS)

# PCSEs:

PCSE<-panelAR(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
                     data=Demos,panelVar="ccode",timeVar="Year",
                     autoCorr="ar1",panelCorrMethod="pcse",
                     rho.na.rm=TRUE)

summary(PCSE)
PCSE$panelStructure$rho

###########################################
# Dynamics:
#
# # Simulating unit roots, etc.
#
# Leaving this commented out for now...
# 
# set.seed(7222009)
# 
# # T=250 N(0,1) errors:
# 
# T <- 250
# Time <- seq(1:T)
# u <- rnorm(T)
# 
# # I(1) process:
# 
# I1<-cumsum(u)
# 
# # AR(1) with rho = 0.9 (done "by hand):
# 
# rho9 <- numeric(T)
# rho9[1] <- u[1] # initialize T=1 to the first error
# for(i in 2:T) {
#   rho9[i] <- (0.9* rho9[i-1]) + u[i]
# }
# 
# # Plot:
# 
# pdf("TSIllustrated.pdf",6,5)
# par(mar=c(4,4,2,2))
# plot(Time,I1,t="l",lwd=2,lty=3,col="red",
#      xlab="Time",ylab="Y")
# lines(Time,rho9,lwd=2,lty=2,col="blue")
# lines(Time,u,lwd=2)
# legend("bottomleft",
#        legend=c("N(0,1) errors","Unit Root","AR(1) with rho=0.9"),
#        col=c("black","red","blue"),bty="n",
#        lwd=c(2,2,2),lty=c(1,3,2))
# dev.off()
########################################
# Panel unit root tests...

lnDemos<-cbind(Demos$ccode,Demos$Year,Demos$lnDemons)
lnDemos<-na.omit(lnDemos)
purtest(lnDemos,exo="trend",test=c("levinlin"))
purtest(lnDemos,exo="trend",test=c("hadri"))
purtest(lnDemos,exo="trend",test=c("ips"))

###############################################
# Some regression models, with dynamics...

# Lagged -dependent-variable model:

PDF$lnDemons.L <- plm::lag(PDF$lnDemons,k=1) # be sure to use the
# -plm- version of -lag-

LDV.fit <- lm(lnDemons~lnDemons.L+POLITY+I(POLITY^2)+
                lnGDP+Monarch+ColdWar,data=PDF)

FD.fit <- plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
              data=PDF, effect="individual",model="fd")

FE.fit <- plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
              data=PDF, effect="individual",model="within")

LDV.FE.fit <- plm(lnDemons~lnDemons.L+POLITY+I(POLITY^2)+lnGDP+
                    Monarch+ColdWar,data=PDF, effect="individual",
                  model="within")

AB.fit<-pgmm(lnDemons~lag(lnDemons,1:2)+POLITY+I(POLITY^2)+
               lag(lnGDP,0:1)+Monarch+ColdWar|lag(lnDemons,2:99),
               data=PDF,effect="individual",model="twosteps")

# Table:

texreg(list(LDV.fit,FD.fit,FE.fit,LDV.FE.fit),
       custom.model.names=c("LDV","First Difference","FE","LDV + FE"),
       custom.coef.names=c("Intercept","Lagged ln(Demonstrations)",
                           "POLITY","POLITY Squared","ln(GDP)",
                           "Monarch","Cold War"),
       digits=3,stars=0.05)

# OPMs:

PDF$POLITYSQ <- PDF$POLITY^2
set.seed(7222009)
OPM.fit <- opm(lnDemons~POLITY+POLITYSQ+lnGDP+Monarch+ColdWar,
               data=PDF,index=c("ccode","Year"),n.samp=1000)

# Ladder plot of estimates & CIs:

pdf("OPM-Ladder.pdf",8,6)
par(mar=c(4,8,2,2))
caterplot(OPM.fit,main=c(""),
          xlab="Parameter Estimate",
          labels=c("Rho","Variance","POLITY",
                   "POLITY Squared","ln(GDP)","Monarchy",
                   "Cold War"))
abline(v=0,lty=2)
dev.off()

# Short- and long-run effects:

SREs<-summary(OPM.fit)$quants[3:7,3]
LREs<-numeric(5)
for(i in 1:5){
  LREs[i]<-quantile(OPM.fit$samples$beta[,i]/(1-OPM.fit$samples$rho),
                    probs=c(0.50))
}

print(cbind(round(SREs,4),round(LREs,4)))

# fin
