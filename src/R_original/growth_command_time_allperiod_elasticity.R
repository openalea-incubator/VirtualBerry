
## prepare the input data ---------------------
   
   Climat <- read.table("grapevine_climat.txt",dec='.',header=T, sep='\t')

   Climat$Doy <- as.numeric(format(as.Date(substr(Climat$TIMESTAMP,1,10),"%d/%m/%Y"),'%j')) 
   Climat$Date <- as.Date(substr(Climat$TIMESTAMP,1,10),"%d/%m/%Y")
   Climat$Hour <- as.numeric(substr(Climat$TIMESTAMP,12,13))
   Climat$min <- as.numeric(substr(Climat$TIMESTAMP,15,16))
   Climat$Temp<- Climat$AirTC_Avg
   Climat$Rh <- (Climat$RH_Max + Climat$RH_Min)/2
  
   DATA <- aggregate(cbind(Temp,Rh,PPFD) ~ Hour+Doy,data=Climat,mean)
   DATA$dap <- DATA$Doy - 126                            ## assume that the pollination date is 05/05/2016
   DATA$Tps <- DATA$dap *24 + DATA$Hour
   DATA$RH <- DATA$Rh
   DATA$vpd <- 0.611*exp(17.27*DATA$Temp/(DATA$Temp+237.3))*(1-DATA$RH/100)

   RGmax <- aggregate(PPFD ~ Doy,data=DATA,max);   RGmin <- aggregate(PPFD ~ Doy,data=DATA,min)
   RGmax$dap <- RGmax$Doy - 126;   RGmin$dap <- RGmin$Doy - 126

   DATA$RGmin <- RGmin$PPFD[match(DATA$dap,RGmin$dap)]
   DATA$RGmax <- RGmax$PPFD[match(DATA$dap,RGmax$dap)]
   DATA <- subset(DATA,Tps >=0,select=c(dap,Tps,Temp,RH,PPFD,vpd,RGmin,RGmax))    
 
   DATA$Cpmin <- 75                                      ## mmol/L
   DATA$Cpmax <- 275                                     ## mmol/L
   DATA$Wpmax <- -1.5                                    ##  bar 
   DATA$Wpmin <- -7.5                                    ##  bar  

   DATA$Cp <- DATA$Cpmin+ (DATA$PPFD- DATA$RGmin)/(DATA$RGmax - DATA$RGmin) *(DATA$Cpmax- DATA$Cpmin)
   DATA$PTLx <- DATA$Wpmax+ (DATA$PPFD- DATA$RGmin)/(DATA$RGmax - DATA$RGmin) * (DATA$Wpmin- DATA$Wpmax)

   if(F){
   DATA$daily <- ifelse(sin((DATA$Tps%%24-5)*pi/15)>0,sin((DATA$Tps%%24-5)*pi/15),0)
   DATA$Cp <- DATA$Cpmin + DATA$daily *(DATA$Cpmax - DATA$Cpmin )
   DATA$PTLx <- DATA$Wpmax + DATA$daily * (DATA$Wpmin - DATA$Wpmax)
         }
   
   DATA_input <- subset(DATA,dap >=4 & dap < 106,select=c(Tps,Temp,RH,Cp,PTLx))

   if(F)
   {
     write.table(DATA_input,"DATA_input.txt",sep='\t',row.names=F)
   }


### ----------------- Measurements ---------------------------------------
  Obs_data <- read.table("grapevine_data.txt",header=T,sep="\t") 

  Obs_data$DMC[Obs_data$DMC >0.35&!is.na(Obs_data$DMC)] <- Obs_data$SDMC[Obs_data$DMC >0.35&!is.na(Obs_data$DMC)]/100
## Obs_data$DMC[Obs_data$DAA==4&!is.na(Obs_data$DMC)] <- Obs_data$SDMC[Obs_data$DAA==4&!is.na(Obs_data$DMC)]/100
  Obs_data$SDMC[Obs_data$FruitID %in% c(9112,9113)] <- Obs_data$DMC[Obs_data$FruitID %in% c(9112,9113)]*100
  Obs_data$SDMC[Obs_data$FruitID == 9115] <- mean(Obs_data$DMC[Obs_data$FruitID %in% c(9111,9112,9113,9114)]*100)
  Obs_data$DW <- Obs_data$FW * Obs_data$DMC


## chose the measurement after veraison-----------

   data <- subset(Obs_data,DAA >= 4 & DAA < 106)

begin <- Sys.time()

## growth simulation
   library(deSolve)
   source("growth_model_desolve.R")

   res <- modgrowth(
     # Parametres relatifs ?l'influx de carbone et a la respiration ------------
     qm293  = 5.9e-5,        # [g sucre g-1MS h-1] Maintenace respiration Coeff
     qg    = 0.02,           # [adim] Growth respiration Coeff [Dai et al.,2010; Ollat and Gaudillere 2000]
     q10   = 1.65,           # [h-1]  Depenance temperature de la respiration

     nuM0  = 8e-3,           # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten) 
     nuM_max  = 8e-3,        # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
     nuM_min  = 0,        # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
     nuM_tau  = 0.50,        # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
     nuM_t  = 48,          # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)

     tStar  = 500,      # [h] Variation de Vmax avec l'age du fruit	
     tauA  = 350,	       # [h] Variation de Vmax avec l'age du fruit

     Cfstar =0.13,         #[g hexose /g H2O] fruit sugar concentration of the inflection point
     kcf = 35,              #[g H2O/ g hexose] proportional to slope at the inflextion point of Ua

     Km    = 0.08,          # [g sacc /g eau] # 
     ps  = 2.7e-3,          # [cm.h-1] Permeability of the composite membran 
     
     # Parametres concernant la plasticit?---------------------------------------

     phi0=0.01,           # [bar-1]
     phiMax  = 0.020,       # [bar-1] 
     phiMin  = 0.010,         # [bar-1]
     kPhi   = 0.015,      # [h-1]
     coefPhi = 30*24,      # [h]

     el0 =  153.2,        # [bar]   elasticity
     kel = 32,            # [cm-1]   elasticity 

     yC0 = 0.5,            # [bar]   threshold value of turgor pressure for growth
     ymax =1.5,
     ymin =0.5,
     Tpsmax = 40,
     Tpsmin = 50,
     
     # Parametres concernant les flux d'eau --------------------------------------
     Hf   = 0.996,	   # [%] Humidite relative dans le fruit

     rho0 = 35,
     rhoMax = 274.9,         # [cm/h]
     rhoMin = 34.65,       # [cm/h]
     kRho = 0.06288,         # [jour-1]
     coefRho = 34.20,        # [day]


     Lp0 = 0.015,
     Lpmax    = 0.045,      #0.015*4,     # [g.cm-2.h-1.bar-1] Conductivity for water transport
     Lpmin    = 3.5e-3*1.6, #1.5,          # [g.cm-2.h-1.bar-1] Conductivity for water transport
     kLp = 0.10,              # [g-1] proportional to the slope at inflextion point of Lp
     coefLp = 55,          # [g] fresh mass of the inflection point

     Lx0  =0,
     Lxmax    = 0.24,    # [g.cm-2.h-1.bar-1] Conductivity for water transport
     Lxmin    = 0,    # [g.cm-2.h-1.bar-1] Conductivity for water transport
     kLx = 0.32,           # [g-1] proportional to the slope at inflextion point of Lp
     coefLx = 33,         # [g] fresh mass of the inflection point

     axp = 3.5e-3,      # [adim] Coefficient for the xylem area of composite membrane
     sigmaX = 1,           # [adim] reflection coefficient xylem

     sigmaP0 = 0.9,         # [adim] reflection coefficient phloem
     tauS = 0.50,      # [h] variation de la reflexion membranaire avec le temps
     t1  = 35,

     piP0 = 8.8,           # [bar] phloem osmotic pressure from other solutions 
     piF0 = 3.5,           # [bar] fruit osmotic pressure from other solutions
     pF0 = 3.0,            # [bar] initial phloem tugor pressure

     apaa = 21.180, #5.91978,
     bpaa = 10.493, #5.0524,
     cpaa = 0.353,  #0.32911,

     # Variables allometriques diverses ------------------------------------------
     
     assrat= 0.62606,  # [adim] % of solubles solids / total : pente
     bssrat= 0.28124,  # [adim] % of solubles solids / total : OO  
     cssrat= 58.203,     # [adim] % of solubles solids / total : OO  
     dssrat= 0.043,     # [adim] % of solubles solids / total : OO  

     rwtr  = 0.6,     # [adim] ratio of water coming from respiration
     gamma = 4.152,  # [adim] emperical parameter
     eta   = 0.707,  # [adim] emperical parameter for the relation of fruit surface and fruit weight

    share1 = 1.721e-02,    # [adim] emperical parameter to calculate stone dry weight
    share2 = 27.99,    # [adim] emperical parameter to calculate stone dry weight

    noysf1 = 0.01772239,   # [adim] emperical parameter to calculate stone fresh weight
    noysf2 = 0.63535422,   # [adim] emperical parameter to calculate stone fresh weight

    aRatio_suc = 0.35698, # [adim] emperical parameter to calculate sucrose ratio
    bRatio_suc = 0.34726, # [adim] emperical parameter to calculate sucrose ratio
    cRatio_suc = 46.66770, # [adim] emperical parameter to calculate sucrose ratio
    dRatio_suc = 0.11757, # [adim] emperical parameter to calculate sucrose ratio
     
     # Variables initiales, duree simulation et contraintes ----------------------
     
     stot0    = 0.00415,      # [g]	Initial dry weight of the fruit
     wtot0    = 0.0118,      # [g] Initial water amount in the fruit

     # Variables control calculations ----------------------
   elasticity = T,      # [T/F] Consider elaticity or not
   calcunuM =  T,         # [T/F] Calculate stone weight or not
   calcuphi = F ,         # [T/F] Calculate phi or not
   calcupiF0 = T,        # [T/F] Calculate piF0 demande
   calcuyC = F,            # [T/F] Calculate yC demande
   calcustone = T,      # [T/F] Calculate stone weight or not
   calcuel = T,            # [T/F] Calculate el or not

   calcuSigmaP = 3,    # [numeric] Calculate sigmaP method
   calcuRho = 4,          # [numeric] Calculate rho demande (the input rho will be replaced by the calculating)
   calcussrat = 3,      # [numeric] Calculate ssrat method
   calcuVmax = 3,        # [numeric] Calculate Vmax method
   calcuLp = 3,            # [numeric] Calculate Lp or not
   calcuLx = 3,            # [numeric] Calculate Lx or not
     
     # Matrice d'entree des conditions aux limites (5 colonnes) ------------------
     DATA  = DATA_input

   )

end <- Sys.time()
 begin - end  


if(F)
{
## graph --------

pdf("grapeallperiod.pdf",height=6,width=8,pointsize=13.5)

   par(lwd=2.0,cex.lab=1.5,cex.axis=1.3)
   plot(res$Tps/24, res$Stot, type='l',xlab='DAA',ylab='Dry weight (g)',ylim=range(res$Stot,data$DW,na.rm=T),lwd=2.5)
   points(data$DAA,data$DW)
   abline(v=c(35,49),lty=2)
   text(20,0.32,"stage I")
   text(42,0.32,"stage II")
   text(65,0.32,"stage III")

   plot(res$Tps/24, res$FFW, type='l',xlab='DAA',ylab='Fresh weight (g)',ylim=range(res$FFW,data$FW,na.rm=T),lwd=2.5)
   points(data$DAA,data$FW)
   abline(v=c(35,49),lty=2)
   text(20,1.5,"stage I")
   text(42,1.5,"stage II")
   text(65,1.5,"stage III")

   plot(res$Tps[-(1:10)]/24, res$pF[-(1:10)],type='l',xlab="DAA",ylab="Tugor pressure (bar)")
   #abline(h=c(0,0.5,1))
   abline(v=c(35,49),lty=2)
   text(20,3.3,"stage I")
   text(42,3.3,"stage II")
   text(65,3.3,"stage III")

  plot(res$Tps/24, res$Up, type='l',xlab='DAA',ylim=range(res$Up,res$Ux),ylab='Water flux (g/h)',col=2) 
  points(res$Tps/24,res$Ux, type='l',col=4)
  abline(h=0)
  abline(v=c(35,49),lty=2)
  legend("topright",c("Phloem","Xylem"),bty='n',pch=NA,lty=1,col=c(2,4))

  plot(res$Tps/24,res$activS,type='l',col=1,ylim=range(res$activS,res$pasflS,res$masflS),xlab='DAA',ylab='Carbon flux (g/h)')
  points(res$Tps/24,res$masflS,type='l',col=2)
  abline(h=0)
  abline(v=c(35,49),lty=2)
  legend("topright",c("active uptake","mass flow"),bty='n',pch=NA,lty=1,col=c(1,2))





res$daa <- trunc(res$Tps/24)

res_day <- aggregate(cbind(Tf,Ux,Up,activS,pasflS,masflS,Rf)~daa,data=res,sum)
res_day$net <- res_day$Tf -res_day$Up -res_day$Ux

res_select <- subset(res_day,daa%in%c(10,15,20,25,35,50,65,80,95),select=c(daa,Tf,Ux,Up))
res_select$Tf <- res_select$Tf*-1000
res_select$Ux <- res_select$Ux*1000
res_select$Up <- res_select$Up*1000
rownames(res_select) <- c(10,15,20,25,35,50,65,80,95)
barplot(t(as.matrix(res_select[,-1])),beside=T,col= c("blue","red","yellow"),ylab='water flux (mg/d/berry)',ylim=c(-200,150),xlab='DAA',
legend.text=c('Transpiration','Xylem flux','Pholem flux'),args.legend=list(x="topleft",bty='n'),space=c(0,1.5))
abline(h=0)
box()

dev.off()
## variables -----
   par(mfrow=c(3,2))

   plot(res$Tps/24, res$Stot, type='l',xlab='DAA',ylab='Dry weight (g)',ylim=range(res$Stot,data$DW,na.rm=T))
   points(data$DAA,data$DW)
   abline(v=c(35,49))

   plot(res$Tps/24, res$FFW, type='l',xlab='DAA',ylab='Fresh weight (g)',ylim=range(res$FFW,data$FW,na.rm=T))
   points(data$DAA,data$FW)
   abline(v=c(35,49),lty=2)


   plot(res$Tps/24, res$Cssm, type='l',xlab='DAA',ylab='Soluble sugars concentration (mol/mL)',
        ylim=range(res$Cssm,data$Osmotic_SUGAR/1e6,na.rm=T))
   points(data$DAA,data$Osmotic_SUGAR/1e6)  
   abline(v=c(35,49))

   plot(res$Tps/24, res$pF,type='l')
   abline(h=c(0,0.5,1))
   abline(v=c(35,49))

  plot(res$Tps/24, res$Up, type='l',xlab='DAA',ylim=range(res$Up,res$Tf,res$Ux)) 
  points(res$Tps/24, res$Tf, type='l',xlab='DAA',col=2)
  points(res$Tps/24,res$Ux, type='l',,col=4)
  abline(h=0)
  abline(v=c(35,49))

  plot(res$Tps/24,res$activS,type='l',col=1,ylim=range(res$activS,res$pasflS,res$masflS))
  points(res$Tps/24,res$masflS,type='l',col=2)
  points(res$Tps/24,res$pasflS,type='l',col=3)
  points(res$Tps/24,res$Rf,type='l',col=4)
  abline(h=0)
  abline(v=c(35,49))

## parameters variation -----
   x11()
   pdf("parameters_variation.pdf",height=8,width=6,pointsize=13.5)
   par(mfrow=c(3,2))
   plot(res$Tps/24, res$el, type='l',xlab='DAA',ylab="elasticity modulus")
   abline(v=c(35,49),lty=2)

   plot(res$Tps/24,res$phi, type='l',xlab='DAA',ylab="extensibility coefficient")
   abline(v=c(35,49),lty=2)
   
   plot(res$Tps/24,res$yC,type='l',xlab='DAA',ylab="yield threshold")
   abline(v=c(35,49),lty=2)


   plot(res$Tps/24,res$nuM,type='l',xlab='DAA',ylab="maximum sugar uptake rate")
   abline(v=c(35,49),lty=2)

   plot(res$Tps/24,res$Lp,type='l',col=1,ylim=range(res$Lp,res$Lx),xlab='DAA',ylab="hydraulic conductivity")
   lines(res$Tps/24,res$Lx,col=2)
   abline(v=c(35,49),lty=2)
   legend("topright",c("Lp","Lx"),bty='n',col=c(1,2),lty=1)

   plot(res$Tps/24,res$sigmap,type='l',xlab='DAA',ylab="phloem reflection coefficient")
   abline(v=c(35,49),lty=2)

   dev.off()

   plot(res$Tps/24,res$Lp*res$Af*3.5e-3,type='l',col=1,ylim=range(res$Lp*res$Af*3.5e-3,res$Lx*res$Af*3.5e-3))
   lines(res$Tps/24,res$Lx*res$Af*3.5e-3,col=2)



  plot(res$Tps/24,res$Ratio_piF,type='l',col=1,ylim=range(res$Ratio_piF,data$Osmotic_other/data$Osmotic_SUGAR,na.rm=T))
  points(data$DAA,data$Osmotic_other/data$Osmotic_SUGAR,type='p')

  plot(res$Tps/24,res$rho,type='l',ylim=range(res$rho,data$rou,na.rm=T))
  points(data$DAA, data$rou)

  plot(res$Tps/24,res$ssrat,type='l',ylim=range(res$ssrat,data$SSMS/100,na.rm=T))
  points(data$DAA,data$SSMS/100)



  plot(res$Tps/24,res$piP,type='l',ylim=range(res$piP,res$pP))
  points(res$Tps/24,res$pP,type='l', col=2)



  plot(res$Tps/24,res$Stot,type='l')
  lines(res$Tps/24,res$DW,type='l',col=2)
  range(res$Stot-res$DW)

  plot(res$Stot,res$Stot-res$DW)

  plot(res$Tps/24,res$Stot-res$DW,type='l',ylim=range(res$Stot-res$DW,data$DWseed,na.rm=T))
  points(data$DAA,data$DWseed,type='p',col=2)

  plot(res$Tps/24,res$Stot-res$DW,type='l',ylim=range(res$Stot-res$DW,data$DWseed,na.rm=T))
  points(data$DAA,data$DWseed)

  plot(res$Tps/24,abs(res$Ux)/(abs(res$Ux)+abs(res$Up)),type='l')


  plot(res$Tps/24,res$difPTLxy,type='l',ylim=range(res$difPTLxy,res$difPTLph))
  lines(res$Tps/24,res$difPTLph,col=2)

  plot(res$Tps/24,res$piP,type='l',col=1,ylim=range(res$piP,res$piF,res$pF))
  lines(res$Tps/24,res$piF,col=2)
  lines(res$Tps/24,res$pP,col=3)
  lines(res$Tps/24,res$pF,col=4)

  plot(res$Tps/24,1-res$Stot/res$FFW,type='l',ylim=range(1-data$DMC,1-res$Stot/res$FFW,na.rm=T))
  points(data$DAA,1-data$DMC,pch=1)

}

