   #29 March 2007, change Lp, fi 
     
berryLpaa <- function(
       qg=0.019, 
       qm25=0.000059,
       Q10=1.6,
       VmaxI=0.0027,
       Km=0.08, 
       tau=167, 
       tstar=241,
       sigma=1,
       PermS=0.0027,
       acnst=0.0045,
       Lpmin=0.0018,
       Lpmax=0.00526,
       Lk=0.04,
       Lh=500,
       fmax=0.005,
       fmin=0.00015,
       fk=0.04,
       fh=180,
       Lx=0,
       Hf=0.996,
       roa=133.46,
       rob=-2.98, 
       Y=2.5,
       Yo=0,                                                        
       R=83,rwtr=0.6,SpVw=18,
       seed=0.055,	s0=0.1745,w0=0.68554,
       PTLmax=-1,PTLmin=-10,
       Cpmax=0.17,Cpmin=0.09,
       tstop=1008,                                              
       te=1008,	
       hourmax=seq(6,1008,by=24), 
       temp=Temp06,RH=RH06 )
			

{ 
 #### respiration#####
#qg      coefficient of growth respiration (non dimensional, gCarbon g-1 DW)
#qm      coefficient of maintenance respiration (gc-1 g-1DW h-1)
#Q10    (non dimensional)
##### water and sugar transport#########
#VmaxI   [g sucrose (g DW)-1 h-1] , 
#Km       Active uptake of dry material (Michaelis Menten)
#tau     (h) characteristic time for inhibitor accumulation
#tstar  (h) kinetic parameter in Michaelis Menten, describing the activity of an inhibitor
#sigma   reflexion coefficient for the composite membrane (non dimensional)
#PermS    permeability of the composite membran for diffusion (cm h-1)
#acnst    coefficient to compute the area of composite membrane (non dimensional)
#Lf      conductivity of the composite membrane for water transport
 ###### transpiration ######
#Hf        Humidity of air space in fruit (non dimensional)
#roa, rob   parameters for calculation of transpiration coefficient

###### Lockart Equation dV/dt=V*fi*(P-Y)#########
#fi Cell wall extensibility coefficient in lockart equation (bar-1 h-1)
#Y threshold of turgor pressure for growth (bar)
#Yo Minimal tugor in the fruit (bar)
########   constant#######
#R gaz constant (cm3 bar mol-1 K-1)
#rwtr ratio of water coming from respiration (non dimensional)
#SpVw Specific water volume (cm3/mol)
###########input values#########
#ssrat % of solubles solids / total.  ratio of soluble sugars in the total carbohydrate pool
#paa ratio of osmotic potentials : Others compounds/sugars (non dimensional)
#stone constant fresh weight of the stone (g)
#s0 initial dry weight of the flesh (g)
#w0 initial water amount in the flesh (g)
#tstop number of hours of the simulation (h)
#Cpmax, Cpmin, PTLmax, PTLmin Minimal and maxi Concentrations and Water potentials in stem
#temp, temperature, RH, relative humidity
   
  s <- vector(length=dim(temp)[1])
  w <- s
  pulpe <- s
  Tf <- s[-length(s)]
  osmf <- Tf
  turgf <- Tf
  permea <- Tf
  eauxy <- Tf
  eauphlo <- Tf
  respi <- Tf
  carin <- Tf	#carbon influx
  carbal <- Tf  #carbon balance
  
  Wx<-Tf
  Rm<-Tf
  Rg<-Tf
  Surface<-Tf
  rou<-Tf
  fai<-Tf
  osmp<-Tf ##phloem  osmotic potential
    Turgf2<-Tf
    Turgf3<-Tf
    dwater<-Tf
    pwf<-Tf ##fruit water potential
    LP<-Tf
    CPG<-Tf
    paa<-Tf # add on 31 July, 2008
    ssrat<-Tf# add on 31 July, 2008
    activS<-Tf
    pasflS<-Tf
    masflS<-Tf
 #INITIALISATIONS
  #ssrat=0.55
  ###assume the sugar inputed is accumulated in the form of sucrose,0.28 is the ssrat[1]
####concentration glc
	##transfer the unit of sugar concentration
	###from g sugar/g solution  to  mol sugar/mol solution
	##for example : Glc and Fru (180.2 g/mol) is the major sugars accumulated in berry , There is Cgg g sugar (glc) and 1-Cgg g water
	### in 1 g solution with a concentration Cgg gsugar/gsolution; then the calculation would be:
	###Cmm (mol sugar/mol solution)=(Cgg/180.2)/(Cgg/180.2+(1-Cgg)/18), get Cmm=Cgg/(10-9Cgg)
	s[1]<-s0
	w[1]<-w0
	pulpe[1] <- s0+w0
	 
 # calculation
 
 for (i in  temp[-dim(temp)[1],"hour"])
 {
 paa[i]<- 0.0842+1.6028*exp(-0.0085*i)+0.0007*i
 ssrat[i]<-0.2783+0.4116*(1-0.987^i)
 #INPUT VALUES
  
  DAILY <- sin(pi*(i-6)/12)	
  if (DAILY <0) DAILY <- 0
  DAYsn <- (1+sin(pi*(i-6)/12))*0.5
  CpINT <- Cpmax-Cpmin
  Cp <- Cpmin + CpINT*DAYsn
  Cpm <- Cp/(19.02-18.02*Cp)#### sucrose is the major sugar in phloem!!!
  Temp <- 273.15+temp[i,"temp06"] 
  OSMp <- 12.53+R*Temp*Cpm/SpVw
  osmp[i]<-OSMp
  PTLint <- PTLmax-PTLmin
  PTLx <- PTLmax-PTLint*DAILY
  Wx[i]<-PTLx/10
  CPG[i]<-Cp
  
  Twght <- w[i] + s[i] + seed
 Af <- 4.152*(Twght)^0.7071
#Af<-9.21*(s[i]+seed*0.5)^0.5108
  Surface[i]<-Af
  Aw <- Af*acnst

  Lp<-Lpmin+(Lpmax-Lpmin)/(1+exp(Lk*(i-Lh)))
  
 ##Lp<-Lpmin+(Lpmax-Lpmin)/(1+10*exp(-exp(-(i-h0)/ki)))
#Lp=0.00972
#if(i<=Lpt1) {Lp<-Lp1} 
#else if (Lpt1<i&i<=Lpt2) {Lp<-Lp2}
#else {Lp<-Lp3}

LP[i]<-Lp
fi<-fmin+(fmax-fmin)/(1+exp(fk*(i-fh)))

   fai[i]<-fi
  eps <- fi*(w[i] + s[i])
   
  #Concentrations
  #ssrat<-0.2783+ssrata*(1-ssratb^i)

	  Cf  <- s[i] / (w[i]+s[i])
	  Css  <- ssrat[i]*s[i] / w[i]
	  Cssm <- Css / (10-9.01*Css)
	  Cssplp <-  ssrat[i]*s[i]/(w[i]+s[i])

  #Transpiration

	  Ha <- RH[i,"RH06"]/100
    Pstar <- 0.00804817*exp(0.0546961*(Temp-273.15))
	  alf <- 18.0*Pstar/(R*Temp) 
	  ro<-(roa*exp(rob*(s[i]+seed*0.5)))  #decrease ro to decrese the diurnal fluctuation 0<fro<1
	  #ro=85
Tf[i] <- ro*alf*Af*(Hf - Ha)
#Tf[i]<0.0028
    rou[i]<-ro      
        
  #Water potentials 	  
	  OSMf <-  (1+paa[i])*(R*Temp*Cssm/SpVw)
	  osmf[i] <- OSMf
	 
  #Turgor regulation
	#BIDA <- Aw*Lf*(2.0*PTLx+OSMf+OSMp-sigma*(OSMp-OSMf))   ##original model  equation
	#BIDB <- eps*Y-Tf[i]                                     ##original model equation
BIDA <- Aw*(Lx*(PTLx+OSMf)+Lp*(PTLx+OSMp-sigma*(OSMp-OSMf)))   ##original model  equation
BIDB <- eps*Y-Tf[i] 
BIDC<-Aw*(Lx+Lp)                                   ##original model equation                                    
	TURGf <- (BIDA + BIDB)/(BIDC+eps)        ##EQUATION 15
	Turgf2[i]<-(BIDA-Tf[i])/BIDC   ### EQUATION 16
	Turgf3[i]<-TURGf            ###Equation15
	    if (TURGf<=Y) 
	    {
	    #BIDC <- PTLx+OSMf-Tf[i]/(2.0*Aw*Lf)
	    #BIDD <-(1.0-sigma)*(OSMp-OSMf)/2.0
	    #TURGf<-BIDC+BIDD
	    TURGf<-(BIDA-Tf[i])/BIDC
	    }
	 if (TURGf<=Yo) TURGf <- Yo
	turgf[i] <- TURGf

	  PTLf <- TURGf - OSMf
	  TURGp <-  PTLx + OSMp
	  difPTLx <- PTLx - PTLf
	  difPTLp <- TURGp - TURGf - sigma*(OSMp- OSMf) 
   pwf[i]<-PTLf
	   
  #Water uptake
 	  Uwx <- Aw*Lx * difPTLx
	  Uwp <- Aw*Lp * difPTLp
	eauxy[i] <- Uwx
	eauphlo[i] <- Uwp
	  	   
  #Sugars uptake
      	  activS[i] <- VmaxI * Cp / ((Km + Cp)*(1+exp((i-tstar)/tau)))
 	  pasflS[i]<- -Aw*PermS*(Cssplp-Cp)/s[i]
	  masflS[i]<- (1.0-sigma)*Uwp*(Cssplp+Cp)*0.5/s[i]
	  Us <- activS[i] + pasflS[i] + masflS[i]
carin[i] <- Us*s[i]


 	  qm <- qm25*(Q10^((Temp-298.15)/10.0))
      	  Ds <- (Us - qm)*s[i] / (1+qg)
	  Rf <- qg * Ds + qm * s[i]

respi[i] <- Rf
Rm[i]<-qm*s[i]
Rg[i]<-qg*Ds
 
  #Water balance
	 Dw <- Uwx + Uwp + rwtr * Rf - Tf[i]
	   #Dw <- Uwp + rwtr * Rf - Tf[i]
	  U <- Uwx + Uwp + Us
	  Up<- Uwp + Us
   dwater[i]<-Dw
 
  #Integration of state variables
	     s[i+1] <- s[i] + Ds
	     w[i+1] <- w[i] + Dw
	     pulpe[i+1] <- s[i+1] + w[i+1]
 }	    
	     
  #RESULTATS 
  	resul <- cbind(temp[,"hour"],s+seed*0.5,w,s+w,s+w+seed,
			c(Tf,NA),c(osmf,NA),c(turgf,NA),c(permea,NA),
			c(eauxy,NA),c(eauphlo,NA),c(respi,NA),
			c(carin,NA),c(Wx,NA),c(Surface,NA),c(Rm,NA),c(Rg,NA),c(rou,NA),c(fai,NA),c(osmp,NA),c(Turgf2,NA),
      c(Turgf3,NA),c(dwater,NA),c(pwf,NA),c(LP,NA),c(CPG,NA),c(activS*s[-1],NA),c(pasflS*s[-1],NA),c(masflS*s[-1],NA))
  	dimnames(resul)[[2]] <- c("1 hour","2 s","3 w","4 pulpe","5 poids frais",
				"6 Tf","7 osmf","8 turgf","9 ro","10 eauxy","11 eauphlo",
				"12 respi","13 carin","14 xylem water potential","15 berry surface area","16 maintenance respiration",
        "17 growth respiration","18 rou","19 cell wall extention","20 phloem osmotic potential",
        "21 Turgor equ16","22 Turgor equ15","23 water change rate","24 fruit water potential","25 Lp","26,phloem sugar concentration",
        "27,active sugar","28,passive sugar","29,mass flow sugar")
  	resul
}  	
  	  




   