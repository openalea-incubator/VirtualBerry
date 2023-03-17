library(deSolve)
## Use ode function to solve the differential equations

modgrowth <- function(

# Parameters related to carbon influx and respiration ------------
    qm293  = qm293,      # [g sucre g-1MS h-1] Maintenace respiration Coeff
    qg    = qg,          # [adim] Growth respiration Coeff
    q10   = q10,         # [h-1]  Dependence temperature of the respiration

    nuM0 = nuM0,           # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
    nuM_max  = nuM_max,          # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
    nuM_min  = nuM_min,          # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
    nuM_tau  = nuM_tau,          # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
    nuM_t  = nuM_t,              # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)

    tStar  = tStar,      # [h] Variation de Vmax avec l'age du fruit	
    tauA  = tauA,	       # [h] Variation de Vmax avec l'age du fruit

    Cfstar =Cfstar,      #[g hexose /g H2O] fruit sugar concentration of the inflection point
    kcf = kcf,           #[g H2O/ g hexose] proportional to slope at the inflextion point of Ua

    Km  = Km,            # [g sucrose /g eau] ## 
    ps  = ps,            # [cm.h-1] Permeability of the composite membrane 

# Parametres concernant la plasticité ---------------------------------------
    phi0 = phi0,
    phiMax  = phiMax,      # [bar-1] 
    phiMin  = phiMin,      # [bar-1]
    kPhi     = kPhi,           # [h-1]
    coefPhi    = coefPhi,      # [h]

    el0 =  el0,            # [bar]   elasticity
    kel = kel,             # [cm-1]   elasticity

    yC0 = yC0,                 # [bar]   threshold value of turgor pressure for growth
    ymax =ymax,              # [bar]  Maxi threshold value of turgor pressure for growth
    ymin =ymin,              # [bar]  Mini threshold value of turgor pressure for growth
    Tpsmax = Tpsmax,         # [jour] Tps en-deça duquel Y = Ymax
    Tpsmin = Tpsmin,         # [jour] Tps au-delà duquel Y = Ymin

# Parametres concernant les flux d'eau --------------------------------------
    Hf    = Hf,	        # [%] Humidite relative dans le fruit

    rho0 = rho0,          # [cm/h]
    rhoMax = rhoMax,      # [cm/h]
    rhoMin = rhoMin,      # [cm/h]
    kRho = kRho,          # [g-1]
    coefRho = coefRho,    # [g]
    
    Lp0      = Lp0,      # [g.cm-2.h-1.bar-1] Conductivity for water transport
    Lpmax    = Lpmax,    # [g.cm-2.h-1.bar-1] Conductivity for water transport
    Lpmin    = Lpmin,    # [g.cm-2.h-1.bar-1] Conductivity for water transport
    kLp = kLp,           # [g-1] proportional to the slope at inflextion point of Lp
    coefLp = coefLp,     # [g] fresh mass of the inflection point
 
    Lx0      = Lx0,      # [g.cm-2.h-1.bar-1] Conductivity for water transport
    Lxmax    = Lxmax,    # [g.cm-2.h-1.bar-1] Conductivity for water transport
    Lxmin    = Lxmin,    # [g.cm-2.h-1.bar-1] Conductivity for water transport
    kLx = kLx,           # [g-1] proportional to the slope at inflextion point of Lp
    coefLx = coefLx,         # [g] fresh mass of the inflection point

    axp = axp,           # [adim]	Coefficient for the xylem or phloem area of composite membrane
    sigmaX = sigmaX,     # [adim] reflection coefficient xylem

    sigmaP0 = sigmaP0,     # [adim] reflection coefficient phloem
    tauS = tauS,         # [h] variation de la reflexion membranaire avec le temps
    t1  = t1, 

    piP0 = piP0,      # [bar] phloem osmotic pressure from other solutions 
    piF0 = piF0,      # [bar] fruit osmotic pressure from other solutions
    pF0 = pF0,        # [bar] initial fruit tugor pressure

    apaa = apaa,      # [adim] ratio of osmotic : Others compounds/sugars : pente
    bpaa = bpaa,      # [adim] ratio of osmotic : Others compounds/sugars : OO
    cpaa = cpaa,      # [adim] ratio of osmotic : Others compounds/sugars : OO
   
# Variables allometriques diverses ------------------------------------------
    assrat= assrat,     # [adim] % of solubles solids / total : pente
    bssrat= bssrat,     # [adim] % of solubles solids / total : OO 
    cssrat= cssrat,     # [adim] % of solubles solids / total : OO  
    dssrat= dssrat,     # [adim] % of solubles solids / total : OO 
 
    rwtr  = rwtr,       # [adim] ratio of water coming from respiration
    gamma =  gamma,     # [adim] emperical parameter
    eta   = eta,        # [adim] emperical parameter for the relation of fruit surface and fruit weight
    
    share1 = share1,    # [adim] emperical parameter to calculate stone dry weight
    share2 = share2,    # [adim] emperical parameter to calculate stone dry weight

    noysf1 = noysf1,   # [adim] emperical parameter to calculate stone fresh weight
    noysf2 = noysf2,   # [adim] emperical parameter to calculate stone fresh weight
    
    aRatio_suc = aRatio_suc, # [adim] emperical parameter to calculate sucrose ratio
    bRatio_suc = bRatio_suc, # [adim] emperical parameter to calculate sucrose ratio
    cRatio_suc = cRatio_suc, # [adim] emperical parameter to calculate sucrose ratio
    dRatio_suc = dRatio_suc, # [adim] emperical parameter to calculate sucrose ratio
   
# Variables initiales, duree simulation et contraintes ----------------------
    stot0    = stot0,         # [g]	Initial dry weight of the fruit
    wtot0    = wtot0,         # [g]	Initial water amount in the fruit

# Variables control calculations ----------------------
   elasticity = elasticity,      # [T/F] Consider elaticity or not
   calcunuM =  calcunuM,         # [T/F] Calculate stone weight or not
   calcuphi = calcuphi ,         # [T/F] Calculate phi or not
   calcupiF0 = calcupiF0,        # [T/F] Calculate piF0 demande
   calcuyC = calcuyC,            # [T/F] Calculate yC demande
   calcustone = calcustone,      # [T/F] Calculate stone weight or not
   calcuel = calcuel,            # [T/F] Calculate el or not

   calcuSigmaP = calcuSigmaP,    # [numeric] Calculate sigmaP method
   calcuRho = calcuRho,          # [numeric] Calculate rho demande (the input rho will be replaced by the calculating)
   calcussrat = calcussrat,      # [numeric] Calculate ssrat method
   calcuVmax = calcuVmax,        # [numeric] Calculate Vmax method
   calcuLp = calcuLp,            # [numeric] Calculate Lp or not
   calcuLx = calcuLx,            # [numeric] Calculate Lx or not

# Matrice d'entree des conditions aux limites (5 colonnes) ------------------
    DATA  = DATA      # [   h       °C      adim    bar     mMol sucre/l seve] 
		          #   Heure Temperature  RH  Potentiel   Concentration
		          #                           Xylème        Sucre
		          #                                     dans le Phloeme
)


{

# describe the input funtions used for growth model #

 times <- DATA$Tps                 #[h]   
 Temp <- DATA$Temp                 #[°C] 
 RH <- DATA$RH                     #[-]
 Cpm <- DATA$Cp                    #[mmol/L]
 WP <- DATA$PTLx                   #[bar]
 inputTemp <- approxfun(times,Temp, rule = 2)
 inputRH <- approxfun(times, RH, rule = 2)
 inputCpm <- approxfun(times, Cpm, rule = 2)
 inputWP <- approxfun(times, WP, rule =2)
 inputfun <- c(inputTemp,inputRH,inputCpm,inputWP)

# definition of the initial values for ode solver #

 s0 = stot0
 w0 = wtot0

 if(calcustone == T)
   { 
     stone0 = share1 * (1 - exp(-share2 * stot0))
     SFstone0 <- noysf1 + noysf2 * stone0

     Sstone = stone0
     s0 = stot0 - stone0 
#    w0 = wtot0 - (SFstone0 - stone0)
   }

 Xstart = c(S =s0, W = w0)
 if(elasticity == T)  Xstart = c(S = s0, W = w0 , pF = pF0)

# definition of the params for ode solver, the parameters used in the growth model #

 params = c(qm293,qg,q10,nuM0,nuM_max,nuM_min,nuM_tau,nuM_t,tStar,tauA,Cfstar,kcf,Km,ps,phi0,phiMax,phiMin,kPhi,coefPhi,
el0,kel,yC0,ymax,ymin,Tpsmax,Tpsmin,Hf,rho0,rhoMax,rhoMin,kRho,coefRho,Lp0,Lpmax,Lpmin,kLp,coefLp,Lx0,Lxmax,Lxmin,kLx,coefLx,
axp,sigmaX,sigmaP0,tauS,t1,piP0,piF0,pF0,apaa,bpaa,cpaa,assrat,bssrat,cssrat,dssrat,rwtr,gamma,eta,share1,share2,noysf1,noysf2,
aRatio_suc,bRatio_suc,cRatio_suc,dRatio_suc,elasticity,calcunuM,calcuRho,calcuphi,calcuLp,calcuLx,calcupiF0,calcuyC,calcustone,
calcuel,calcuSigmaP,calcussrat,calcuVmax)


# definition of the fruit growth model, or the differential functions used for ode solver #

 Growth <- function(t,x,parms,inputfun)

 { with(as.list(c(x,parms)),{

  # input
   Temp <- inputfun[[1]](t)
   RH <- inputfun[[2]](t)
   Cpm <- inputfun[[3]](t)
   WP <- inputfun[[4]](t)

  # Constantes physiques -------------------------------------------------
   R           <- 83       # [cm3.bar.mol-1.K-1] Gaz constant
   SpVw        <- 18       # [cm3.mol-1] Specific water volume
   Ms.eau      <- 18       # [g.mol-1]	 Masse molaire eau
   Ms.sucrose  <- 342.3    # [g.mol-1] Masse molaire saccharose
   Ms.glucose  <- 180      # [g.mol-1] Masse molaire glucose
   Dw          <- 1        # [g.cm -3] water density
   Ds          <- 1.6        # [g.cm -3] dry mass density ???
   
  # Climate characteristics  ----------------------------
      Temp   <- 273.15 + Temp		                           # [K] Temperature
      RH     <- RH/100	                                        # [%]  Relative Humidity 
      Psat   <- 0.00804817 * exp(0.0546961*(Temp-273.15))          # [bar] saturated vapor pressure

  # calculate the weights of the stone -----------------

    Stot <- S
    if(calcustone == T)
      { 
        Stot <- S + Sstone
        Sstone <<- share1 * (1 - exp(-share2 * Stot))
        SFstone <- noysf1 + noysf2 * Sstone
      }

  # Fruit characteristics ---------------------------
      FFW  <- W + S 	                                       # [g] fruit fresh weight
      if(calcustone == T)
        { FFW <- W + S + SFstone }
      Af <- gamma * FFW^eta	                                 # [cm2] fruit surface
      Ax <- Af * axp		                                 # [cm2]  Surf.Memb. composite
      Ap <- Af * axp		                                 # [cm2]  Surf.Memb. composite

  # Xylem and phloem characteristics
      Cpm <- (Cpm*1e-3/1000)                                   # [mol sucre/cm3 solution]
      Cpg <- Cpm * Ms.sucrose                                  # [g sucre/g eau] # suppose sucre phloeme = saccharose
                                        
      piP.suc <-  R * Temp * Cpm                               # [bar] Pression Osmotic des sucres phloem
      piP <- piP0 + piP.suc                                    # [bar] Pression Osmotic phloem
      pP <- WP + piP		                                 # [bar] Turg. phloeme
   
  # Concentrations and Potential in Fruit ------------------------
     
     ssrat <- switch(calcussrat,
                      0.52,                                                # 1
                      assrat*(1-exp(-bssrat * t/24)) + cssrat,             # 2[g sucre solubles/g MS]  
                      assrat/(1+exp(-bssrat*(t/24-cssrat)))+dssrat )       # 3

     Cf <- S / (W + S)		                                # [g MS/gMF]   Teneur en MS
     Cssplp <- ssrat * S / (W + S)	                          # [g sucres solubles/g MF]
     Css <- ssrat * S / W	                                # [g sucres solubles/g eau]


     Ratio_suc <- aRatio_suc / (1+ exp(bRatio_suc*(t/24 - cRatio_suc))) + dRatio_suc  # calculate the ratio of sucrose to soluble sugars


     Cssm <-  Css*(1- Ratio_suc)/Ms.glucose + Css *Ratio_suc /Ms.sucrose 	# [mol sucres/g eau] # suppose sucre fruit = glucose

     piF.suc <- R * Temp * Cssm                               # [bar] PO fruit sucres
     piF0 <- piF0
     paa <-  apaa*exp(-bpaa*Cssm*1000)+ cpaa                  ## Cssm *1000  mol/mL  to mol/L

     if(calcupiF0 == T)
       { piF0 <- paa * piF.suc }

     piF <- piF.suc + piF0                                     # [bar] PO fruit

  # Caractéristiques plasticité -----------------------------
     phi <- phi0
     if(calcuphi == T)
        { phi <- phiMin + ( phiMax - phiMin )/( 1 + exp( kPhi * (t-coefPhi ) ) )  }    # [bar-1.h-1]   Plasticite
                 
  # Calculate fruit volume  -----------------------------
     vol <- W/Dw + Stot/Ds		                               # [cm3]  fruit volume

  # Caractéristiques Reflexion membranaire ------------------
     sigmap <- switch(calcuSigmaP,
                       sigmaP0,                            #1
                       sigmaP0 - exp(-tauS *(t-t1)^2),      #2
                       sigmaP0/(1+exp(-tauS*(t/24-t1)))     #3
                      )     

     if(sigmap < 0) sigmap <- 0
      	                                          # [adim]  ReflexionMembranaire
  
  # Caractéristiques conductance ----------------------------
     rho <- switch(calcuRho,
                       rho0,                                          #1
                       rhoMax*exp(-kRho*t/24)+rhoMin,                 #2
                       rhoMax/(1+exp(kRho*(FFW-coefRho)))+rhoMin,     #3 [cm.h-1]  Conductance Fruit
                       rhoMax/(1+exp(kRho*(t/24-coefRho)))+rhoMin     #4
                      )     

 
 # Calculate yC threshold pressure
     yC <- yC0

     if(calcuyC == T)
       {
         Pte     <- (ymax-ymin)/(Tpsmax-Tpsmin)
         OrdOri  <- ymax - Pte * Tpsmax
         yC      <- Pte * t/24 + OrdOri
         if ( (t/24) < Tpsmax) { yC <- ymax }
         if ( (t/24) > Tpsmin) { yC <- ymin }
        }

 # Transpiration -------------------------
     alf <- Ms.eau / (R * Temp) * Psat               # [adim]
     Tf <- rho * alf * Af * (Hf - RH)                # [cm3.h-1] Transp
     Tfstar <- Tf/Af			             # [cm.h-1]  TranspSurf

  # Maintenance respiration -------------------------

     qm  <- qm293*(q10^((Temp-293.15)/10.0))     # [g sucre g MS-1.h-1] RespMaintenance      

  # Calculate the concentration and osmotic and turgor pressure in the vascular of fruit-------------------------
      Cv <- Cpg
      piV <- piP
      pV <- pP
      Cs <- ifelse(max(Up,0),Cv,Css)   #(Css + Cv) * 0.5 

  # Caractéristiques influx sucre ---------------------------
     nuM <- nuM0
     if(calcunuM == T)
        { nuM <- nuM_max/(1+exp(-nuM_tau*(t/24- nuM_t)))+ nuM_min }


     Vmax  <- switch(calcuVmax,
                       nuM,                                      #1
                       nuM/(1+exp((t-tStar)/tauA)),              #2 [g sucre/g MS/h]
                       nuM/(1+exp((Css -Cfstar )*kcf))           #3 [g sucre/g MS/h]
                     ) 
      
     activS <- Stot * Vmax * Cv / (Km + Cv)            # [g sucre.h-1] 

  # Calculate the hydrolic conductivity of phloem and xylem ---------
 
     Lp <- switch(calcuLp,
                       Lp0,                                                    #1
                       Lpmin + (Lpmax - Lpmin)/(1+exp(kLp * (FFW -coefLp))),   #2[g.cm-2.h-1.bar-1] Conductivity for water transport
                       Lpmin + (Lpmax - Lpmin)/(1+exp(kLp * (t/24 -coefLp)))   #3[g.cm-2.h-1.bar-1] Conductivity for water transport

                   ) 

     
     Lx <- switch(calcuLx,
                       Lx0,                                                    #1
                       Lxmin + (Lxmax - Lxmin)/(1+exp(kLx * (FFW -coefLx))),   #2[g.cm-2.h-1.bar-1] Conductivity for water transport
                       Lxmin + (Lxmax - Lxmin)/(1+exp(kLx * (t/24 -coefLx)))   #3[g.cm-2.h-1.bar-1] Conductivity for water transport

                   ) 


  # Water potentials 	   -----------------
     # Turgor regulation
        prod1 =  (Ax * Lx + Ap * Lp) 
        prod2 =  (Ax * Lx + Ap * Lp * sigmap)
        prod3 =  activS +(Ap * ps * (Cv - Css))
        prod4 =  Cs *(1 - sigmap) * Ap * Lp 
     # calculate aCoef
       div1 = (rwtr * (qm*Stot  + qg * prod3) - Tf * (1 + qg))/Dw
       div2 = (prod3 - qm*Stot)/Ds
       aCoef = (div1 + div2)/(1+qg)
     # calculate bCoef
       div1 = (rwtr * qg * prod4 +prod1 * (1 + qg))/Dw
       div2 = prod4/Ds
       bCoef = (div1 + div2)/(1+qg)
     # calculate cCoef
       prod4 = prod4 *sigmap
       div1 = (rwtr * qg * prod4 +prod2 * (1 + qg))/Dw
       div2 = prod4/Ds
       cCoef = (div1 + div2)/(1+qg)

  # calculate fruit turgor 
     if(elasticity == F) {
     pF1 = ( aCoef + bCoef * pV - cCoef * (piV - piF) + vol * phi*yC )/(bCoef + vol * phi)
     pF2 = ifelse(max(pV + ( aCoef - cCoef * ( piV - piF))/ bCoef,0),pV + ( aCoef - cCoef * ( piV - piF))/ bCoef,0) 
     pF = ifelse(max(pF1-yC,0),pF1,pF2) }

     PTLf <- pF - piF	                          # [bar] PotHyd Total Fruit
     difPTLx <- WP - PTLf	                          # [bar] DiffPotHyd Xy-Fruit
     difPTLp <- pP - pF - sigmap* (piP- piF)	      # [bar] DiffPotHyd Ph-Fruit
	                    
  # Water uptake ----------------------------------------
     Ux <- Ax * Lx * (pV - pF - sigmaX * (piV - piF) )		# [g eau.h-1] flux d'eau xylemien
     Up <<- Ap * Lp * (pV - pF - sigmap * (piV - piF) )		# [g eau.h-1] flux d'eau phloemien
	  	   
  # Sugars uptake ---------------------------------------
    pasflS <- Ap * ps * ( Cv- Css)             # [g sucre.h-1]
    masflS <- (1.0-sigmap) * Up * Cs	     # [g sucre.h-1]
    Us <- activS + pasflS + masflS             # [g sucre.h-1]

  # Sugars respiration-------------------------------------
    Rm <- qm *Stot                                                   # [g sucre.h-1] RespMaintenance
    Rg <- ifelse(max(Us - Rm,0),qg *(Us - Rm)/(1+qg),0)	      # [g sucre.h-1] RespCroissance
    Rf <- Rm + Rg		
				             
  # Carbon and Water balance -----------------------------------------
    dS <- Us - Rf
    dW <- Ux + Up + rwtr * Rf - Tf            # [g eau.h-1]  Croissance Eau

    if(calcustone == T){
    dStot <- Us -Rf
    dSstone <- share1 * share2 * exp(-share2 * Stot) * dStot
    dS <- dStot - dSstone

    dWtot <- Ux + Up + rwtr * Rf - Tf
    dSFstone <- noysf2 * dSstone
    dW <- dWtot - (dSFstone -dSstone)

                       }

  # Caractéristiques elasticity  -----------------------------
     el  <- el0
     if(calcuel == T)
        {  
           el <- kel * pF *  (3*vol/4/pi)^(1/3)    # [bar]   elasticity
        }

   if(elasticity == T) {
       dpF = ( el/vol )* ((Ux + Up + rwtr * Rf - Tf)/Dw + (Us - Rf)/Ds)
       if(pF > yC) dpF = dpF - el * phi * (pF - yC) }    

   res <- list(c(dS,dW), Stot,FFW,Af,vol,piP,pP,Cssm,piF,pF,phi,sigmap,rho,Tf,Ux,Up,activS,pasflS,masflS,Rf,
               WP,difPTLx,difPTLp,ssrat,paa,yC,Vmax,Lp,Lx,nuM,el)

   if(elasticity == T)
   res <- list(c(dS,dW,dpF), Stot,FFW,Af,vol,piP,pP,Cssm,piF,phi,sigmap,rho,Tf,Ux,Up,activS,pasflS,masflS,Rf,
               WP,difPTLx,difPTLp,ssrat,paa,yC,Vmax,Lp,Lx,nuM,el)

   return(res)})
}

  Up <- 0.01
  output <- ode(y = Xstart, times = times, func = Growth, parms = params, inputfun= inputfun,method="rk4")
  output <- data.frame(output)
  names(output) <- c("Tps","DW","WW","Stot","FFW","Af","vol","piP","pP","Cssm","piF","pF","phi","sigmap","rho","Tf","Ux","Up","activS",
                     "pasflS","masflS","Rf","WP","difPTLxy","difPTLph","ssrat","Ratio_piF","yC","Vmax","Lp","Lx","nuM","el")
  if(elasticity == T)
  names(output) <- c("Tps","DW","WW","pF","Stot","FFW","Af","vol","piP","pP","Cssm","piF","phi","sigmap","rho","Tf","Ux","Up","activS",
                    "pasflS","masflS","Rf","WP","difPTLxy","difPTLph","ssrat","Ratio_piF","yC","Vmax","Lp","Lx","nuM","el")

  return(output)

}






