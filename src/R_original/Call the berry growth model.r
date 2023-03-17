#28/11/2011### redo the figure used in ISHS
### The model######

########################################
########################################
#########Experimental data##############
########################################
########################################

#Define the working directory
getwd()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

wd<-"F:/egfv2022/2 recherche/modélisation/fruit growth model (zhanwu)/R fruit growth model"
setwd(wd)

#Read climate data for model inputs
climate<-read.csv("grapevine_climat.csv",header=T,sep=",",dec=".")

#Define climate variables for model inputs
Tem<-climate$Temp
Hour<-climate$hour
RH<-climate$RH
Name<-as.character(Hour)
Temp06<-structure(.Data=list(temp06=Tem,hour=Hour),class="data.frame",row.names=Name)
RH06<-structure(.Data=list(RH06=RH),class="data.frame",row.names=Name)

# read experimental data for comparing observed and simulated values
hardata<-read.csv("harvestdata.csv",header=TRUE,skip=2,sep=",",dec=".")


#load the model
source("The berry growth model.r")

#Run model for treatment CK
   harck<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.25,Cpmin=0.17,sigma=0.9)
   #Run model for the treatment
   harrl<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.16,Cpmin=0.08,sigma=0.9)
   #making figures
  par(mfrow=c(2,2))
 #fresh weight
 plot(hardata$hour,hardata$CKFW,col=1,ylim=c(0.35,1.3),lwd=2)
 points(hardata$hour,hardata$CKM,col=6,pch=8,lwd=2)
 lines(harck[,1],harck[,5],col=2)
 plot(hardata$hour,hardata$RLFW,col=3,ylim=c(0.35,1.3))
 points(hardata$hour,hardata$RLM,col=6,pch=8)
   lines(harrl[,1],harrl[,5],col=2)  
     
# dry weight
 plot(hardata$hour,hardata$ckdw,col=1,ylim=c(0.05,0.3))
  points(hardata$hour,hardata$ckm,col=6,pch=8)
  lines(harck[,1],harck[,2],col=2)
  plot(hardata$hour,hardata$rldw,col=3,ylim=c(0.05,0.3))
   points(hardata$hour,hardata$rlm,col=6,pch=8)
    lines(harrl[,1],harrl[,2],col=2)
  
##### 16/06/2014############################
#### simulating water stress effect#########
 harck<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.25,Cpmin=0.17,sigma=0.9)
harws1<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.25*0.95,Cpmin=0.17*0.95,sigma=0.9,PTLmax=-2,PTLmin=-12)
harws2<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.25*0.93,Cpmin=0.17*0.93,sigma=0.9,PTLmax=-4,PTLmin=-14)
 harws3<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.25*0.92,Cpmin=0.17*0.92,sigma=0.9,PTLmax=-6,PTLmin=-16)
harws<-berryLpaa(tau=30,tstar=490,Lh=200,Lk=0.011,Lpmin=0.0034,Lpmax=0.023,fh=200,fk=0.5,fmin=0.01,fmax=0.01,s0=0.06085,w0=0.4341,acnst=0.0035,Cpmax=0.25*0.9,Cpmin=0.17*0.9,sigma=0.9,PTLmax=-8,PTLmin=-20)
   
 par(mfrow=c(2,2))
 

plot(harck[1:24,1],harck[1:24,14],type="l",ylim=c(-2.01,-0.099),lwd=2,xlab="Hours of the day",ylab="Stem water potential(MPa)",cex.axis=1.8,cex.lab=1.8)
lines(harws[,1],harws[,14],col=2,lwd=2)
 lines(harws1[,1],harws1[,14],col=colors()[52],lwd=1)
lines(harws2[,1],harws2[,14],col=colors()[32],lwd=1)
lines(harws3[,1],harws3[,14],col=colors()[33],lwd=1)


plot(harck[1:24,1],harck[1:24,26],type="l",ylim=c(0.15,0.255),lwd=2,xlab="Hours of the day",ylab="Phloem sugar concentration (g/g)",cex.axis=1.8,cex.lab=1.8)
lines(harws[,1],harws[,26],col=2,lwd=2)
 lines(harws1[,1],harws1[,26],col=colors()[52],lwd=1)
lines(harws2[,1],harws2[,26],col=colors()[32],lwd=1)
lines(harws3[,1],harws3[,26],col=colors()[33],lwd=1)


plot(harck[,1],harck[,5],type="l",xaxt="n",lwd=2,xlab="Days after flowering",ylab="Berry Fresh weight (g)",cex.axis=1.8,cex.lab=1.8,ylim=c(0.45,1.1))
  axis(1,seq(1,1008,by=240),seq(60,100,by=10),cex.axis=2)
 lines(harws[,1],harws[,5],col=2,lwd=2)
 lines(harws1[,1],harws1[,5],col=colors()[52],lwd=1)
lines(harws2[,1],harws2[,5],col=colors()[32],lwd=1)
lines(harws3[,1],harws3[,5],col=colors()[33],lwd=1)


plot(harck[,1],harck[,2],type="l",xaxt="n",lwd=2,xlab="Days after flowering",ylab="Berry Fresh weight (g)",cex.axis=1.8,cex.lab=1.8,ylim=c(0.08,0.25))
  axis(1,seq(1,1008,by=240),seq(60,100,by=10),cex.axis=2)
 lines(harws[,1],harws[,2],col=2,lwd=2)
 lines(harws1[,1],harws1[,2],col=colors()[52],lwd=1)
lines(harws2[,1],harws2[,2],col=colors()[32],lwd=1)
lines(harws3[,1],harws3[,2],col=colors()[33],lwd=1)
