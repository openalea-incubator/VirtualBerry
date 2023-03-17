import os
import rpy2.robjects as robjects

os.environ['R_HOME']='C:\\Users\\delivore\\Documents\\R\\R-4.1.2'
r_source = robjects.r['source']
robjects.r("library(deSolve)")
r_source("growth_model_desolve.R")

def modgrowth_wrapper(
        # Parameters related to carbon influx and respiration ------------
        qm293,  # [g sucre g-1MS h-1] Maintenace respiration Coeff
        qg,  # [adim] Growth respiration Coeff
        q10,  # [h-1]  Dependence temperature of the respiration

        nuM0,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_max,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_min,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_tau,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_t,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)

        tStar,  # [h] Variation de Vmax avec l'age du fruit
        tauA,  # [h] Variation de Vmax avec l'age du fruit

        Cfstar,  # [g hexose /g H2O] fruit sugar concentration of the inflection point
        kcf,  # [g H2O/ g hexose] proportional to slope at the inflextion point of Ua

        Km,  # [g sucrose /g eau] ##
        ps,  # [cm.h-1] Permeability of the composite membrane

        # Parametres concernant la plasticité ---------------------------------------
        phi0,
        phiMax,  # [bar-1]
        phiMin,  # [bar-1]
        kPhi,  # [h-1]
        coefPhi,  # [h]

        el0,  # [bar]   elasticity
        kel,  # [cm-1]   elasticity

        yC0,  # [bar]   threshold value of turgor pressure for growth
        ymax,  # [bar]  Maxi threshold value of turgor pressure for growth
        ymin,  # [bar]  Mini threshold value of turgor pressure for growth
        Tpsmax,  # [jour] Tps en-deça duquel Y = Ymax
        Tpsmin,  # [jour] Tps au-delà duquel Y = Ymin

        # Parametres concernant les flux d'eau --------------------------------------
        Hf,  # [%] Humidite relative dans le fruit

        rho0,  # [cm/h]
        rhoMax,  # [cm/h]
        rhoMin,  # [cm/h]
        kRho,  # [g-1]
        coefRho,  # [g]

        Lp0,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lpmax,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lpmin,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        kLp,  # [g-1] proportional to the slope at inflextion point of Lp
        coefLp,  # [g] fresh mass of the inflection point

        Lx0,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lxmax,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lxmin,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        kLx,  # [g-1] proportional to the slope at inflextion point of Lp
        coefLx,  # [g] fresh mass of the inflection point

        axp,  # [adim]	Coefficient for the xylem or phloem area of composite membrane
        sigmaX,  # [adim] reflection coefficient xylem

        sigmaP0,  # [adim] reflection coefficient phloem
        tauS,  # [h] variation de la reflexion membranaire avec le temps
        t1,

        piP0,  # [bar] phloem osmotic pressure from other solutions
        piF0,  # [bar] fruit osmotic pressure from other solutions
        pF0,  # [bar] initial fruit tugor pressure

        apaa,  # [adim] ratio of osmotic : Others compounds/sugars : pente
        bpaa,  # [adim] ratio of osmotic : Others compounds/sugars : OO
        cpaa,  # [adim] ratio of osmotic : Others compounds/sugars : OO

        # Variables allometriques diverses ------------------------------------------
        assrat,  # [adim] % of solubles solids / total : pente
        bssrat,  # [adim] % of solubles solids / total : OO
        cssrat,  # [adim] % of solubles solids / total : OO
        dssrat,  # [adim] % of solubles solids / total : OO

        rwtr,  # [adim] ratio of water coming from respiration
        gamma,  # [adim] emperical parameter
        eta,  # [adim] emperical parameter for the relation of fruit surface and fruit weight

        share1,  # [adim] emperical parameter to calculate stone dry weight
        share2,  # [adim] emperical parameter to calculate stone dry weight

        noysf1,  # [adim] emperical parameter to calculate stone fresh weight
        noysf2,  # [adim] emperical parameter to calculate stone fresh weight

        aRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio
        bRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio
        cRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio
        dRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio

        # Variables initiales, duree simulation et contraintes ----------------------
        stot0,  # [g]	Initial dry weight of the fruit
        wtot0,  # [g]	Initial water amount in the fruit

        # Variables control calculations ----------------------
        elasticity,  # [T/F] Consider elaticity or not
        calcunuM,  # [T/F] Calculate stone weight or not
        calcuphi,  # [T/F] Calculate phi or not
        calcupiF0,  # [T/F] Calculate piF0 demande
        calcuyC,  # [T/F] Calculate yC demande
        calcustone,  # [T/F] Calculate stone weight or not
        calcuel,  # [T/F] Calculate el or not

        calcuSigmaP,  # [numeric] Calculate sigmaP method
        calcuRho,  # [numeric] Calculate rho demande (the input rho will be replaced by the calculating)
        calcussrat,  # [numeric] Calculate ssrat method
        calcuVmax,  # [numeric] Calculate Vmax method
        calcuLp,  # [numeric] Calculate Lp or not
        calcuLx,  # [numeric] Calculate Lx or not

        # Matrice d'entree des conditions aux limites (5 colonnes) ------------------
        DATA  # [   h       °C      adim    bar     mMol sucre/l seve]
        #   Heure Temperature  RH  Potentiel   Concentration
        #                           Xylème        Sucre
        #                                     dans le Phloeme
):
    res = robjects.r['modgrowth'](
        # Parameters related to carbon influx and respiration ------------
        qm293=qm293,  # [g sucre g-1MS h-1] Maintenace respiration Coeff
        qg=qg,  # [adim] Growth respiration Coeff
        q10=q10,  # [h-1]  Dependence temperature of the respiration

        nuM0=nuM0,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_max=nuM_max,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_min=nuM_min,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_tau=nuM_tau,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)
        nuM_t=nuM_t,  # [(g sucrose).(g DW)-1. h-1]  Active uptake of sucrose (Mich-Menten)

        tStar=tStar,  # [h] Variation de Vmax avec l'age du fruit
        tauA=tauA,  # [h] Variation de Vmax avec l'age du fruit

        Cfstar=Cfstar,  # [g hexose /g H2O] fruit sugar concentration of the inflection point
        kcf=kcf,  # [g H2O/ g hexose] proportional to slope at the inflextion point of Ua

        Km=Km,  # [g sucrose /g eau] ##
        ps=ps,  # [cm.h-1] Permeability of the composite membrane

        # Parametres concernant la plasticité ---------------------------------------
        phi0=phi0,
        phiMax=phiMax,  # [bar-1]
        phiMin=phiMin,  # [bar-1]
        kPhi=kPhi,  # [h-1]
        coefPhi=coefPhi,  # [h]

        el0=el0,  # [bar]   elasticity
        kel=kel,  # [cm-1]   elasticity

        yC0=yC0,  # [bar]   threshold value of turgor pressure for growth
        ymax=ymax,  # [bar]  Maxi threshold value of turgor pressure for growth
        ymin=ymin,  # [bar]  Mini threshold value of turgor pressure for growth
        Tpsmax=Tpsmax,  # [jour] Tps en-deça duquel Y = Ymax
        Tpsmin=Tpsmin,  # [jour] Tps au-delà duquel Y = Ymin

        # Parametres concernant les flux d'eau --------------------------------------
        Hf=Hf,  # [%] Humidite relative dans le fruit

        rho0=rho0,  # [cm/h]
        rhoMax=rhoMax,  # [cm/h]
        rhoMin=rhoMin,  # [cm/h]
        kRho=kRho,  # [g-1]
        coefRho=coefRho,  # [g]

        Lp0=Lp0,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lpmax=Lpmax,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lpmin=Lpmin,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        kLp=kLp,  # [g-1] proportional to the slope at inflextion point of Lp
        coefLp=coefLp,  # [g] fresh mass of the inflection point

        Lx0=Lx0,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lxmax=Lxmax,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        Lxmin=Lxmin,  # [g.cm-2.h-1.bar-1] Conductivity for water transport
        kLx=kLx,  # [g-1] proportional to the slope at inflextion point of Lp
        coefLx=coefLx,  # [g] fresh mass of the inflection point

        axp=axp,  # [adim]	Coefficient for the xylem or phloem area of composite membrane
        sigmaX=sigmaX,  # [adim] reflection coefficient xylem

        sigmaP0=sigmaP0,  # [adim] reflection coefficient phloem
        tauS=tauS,  # [h] variation de la reflexion membranaire avec le temps
        t1=t1,

        piP0=piP0,  # [bar] phloem osmotic pressure from other solutions
        piF0=piF0,  # [bar] fruit osmotic pressure from other solutions
        pF0=pF0,  # [bar] initial fruit tugor pressure

        apaa=apaa,  # [adim] ratio of osmotic : Others compounds/sugars : pente
        bpaa=bpaa,  # [adim] ratio of osmotic : Others compounds/sugars : OO
        cpaa=cpaa,  # [adim] ratio of osmotic : Others compounds/sugars : OO

        # Variables allometriques diverses ------------------------------------------
        assrat=assrat,  # [adim] % of solubles solids / total : pente
        bssrat=bssrat,  # [adim] % of solubles solids / total : OO
        cssrat=cssrat,  # [adim] % of solubles solids / total : OO
        dssrat=dssrat,  # [adim] % of solubles solids / total : OO

        rwtr=rwtr,  # [adim] ratio of water coming from respiration
        gamma=gamma,  # [adim] emperical parameter
        eta=eta,  # [adim] emperical parameter for the relation of fruit surface and fruit weight

        share1=share1,  # [adim] emperical parameter to calculate stone dry weight
        share2=share2,  # [adim] emperical parameter to calculate stone dry weight

        noysf1=noysf1,  # [adim] emperical parameter to calculate stone fresh weight
        noysf2=noysf2,  # [adim] emperical parameter to calculate stone fresh weight

        aRatio_suc=aRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio
        bRatio_suc=bRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio
        cRatio_suc=cRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio
        dRatio_suc=dRatio_suc,  # [adim] emperical parameter to calculate sucrose ratio

        # Variables initiales, duree simulation et contraintes ----------------------
        stot0=stot0,  # [g]	Initial dry weight of the fruit
        wtot0=wtot0,  # [g]	Initial water amount in the fruit

        # Variables control calculations ----------------------
        elasticity=elasticity,  # [T/F] Consider elaticity or not
        calcunuM=calcunuM,  # [T/F] Calculate stone weight or not
        calcuphi=calcuphi,  # [T/F] Calculate phi or not
        calcupiF0=calcupiF0,  # [T/F] Calculate piF0 demande
        calcuyC=calcuyC,  # [T/F] Calculate yC demande
        calcustone=calcustone,  # [T/F] Calculate stone weight or not
        calcuel=calcuel,  # [T/F] Calculate el or not

        calcuSigmaP=calcuSigmaP,  # [numeric] Calculate sigmaP method
        calcuRho=calcuRho,  # [numeric] Calculate rho demande (the input rho will be replaced by the calculating)
        calcussrat=calcussrat,  # [numeric] Calculate ssrat method
        calcuVmax=calcuVmax,  # [numeric] Calculate Vmax method
        calcuLp=calcuLp,  # [numeric] Calculate Lp or not
        calcuLx=calcuLx,  # [numeric] Calculate Lx or not

        # Matrice d'entree des conditions aux limites (5 colonnes) ------------------
        DATA=DATA  # [   h       °C      adim    bar     mMol sucre/l seve]
        #   Heure Temperature  RH  Potentiel   Concentration
        #                           Xylème        Sucre
        #

        )
    return(res)
