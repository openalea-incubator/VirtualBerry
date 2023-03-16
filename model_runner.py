import os
import rpy2.robjects as robjects
import model_wrapper as mw

robjects.r("DATA_input = read.table('DATA_input.txt',sep='\t',header= T) ")

DATA_input = robjects.r['DATA_input']

# Parametres relatifs ?l'influx de carbone et a la respiration ------------
qm293 = 5.9e-5  # [g sucre g-1MS h-1] Maintenace respiration Coeff
qg = 0.02  # [adim] Growth respiration Coeff [Dai et al.,2010; Ollat and Gaudillere 2000]
q10 = 1.65  # [h-1]  Dependance temperature de la respiration

nuM0 = 8e-3  # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
nuM_max = 8e-3 # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
nuM_min = 0  # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
nuM_tau = 0.50  # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)
nuM_t = 48  # [(g sacc.).(g DW)-1. h-1]  Active uptake of dry material (Mich-Menten)

tStar = 500  # [h] Variation de Vmax avec l'age du fruit
tauA = 350  # [h] Variation de Vmax avec l'age du fruit

Cfstar = 0.13  # [g hexose /g H2O] fruit sugar concentration of the inflection point
kcf = 35  # [g H2O/ g hexose] proportional to slope at the inflextion point of Ua

Km = 0.08  # [g sacc /g eau] #
ps = 2.7e-3  # [cm.h-1] Permeability of the composite membran

# Parametres concernant la plasticit?---------------------------------------

phi0 = 0.01  # [bar-1]
phiMax = 0.020  # [bar-1]
phiMin = 0.010  # [bar-1]
kPhi = 0.015  # [h-1]
coefPhi = 30 * 24  # [h]

el0 = 153.2  # [bar]   elasticity
kel = 32  # [cm-1]   elasticity

yC0 = 0.5  # [bar]   threshold value of turgor pressure for growth
ymax = 1.5
ymin = 0.5
Tpsmax = 40
Tpsmin = 50

# Parametres concernant les flux d'eau --------------------------------------
Hf = 0.996 # [%] Humidite relative dans le fruit

rho0 = 35
rhoMax = 274.9  # [cm/h]
rhoMin = 34.65  # [cm/h]
kRho = 0.06288  # [jour-1]
coefRho = 34.20  # [day]

Lp0 = 0.015
Lpmax = 0.045  # 0.015*4,     # [g.cm-2.h-1.bar-1] Conductivity for water transport
Lpmin = 3.5e-3 * 1.6 # 1.5,          # [g.cm-2.h-1.bar-1] Conductivity for water transport
kLp = 0.10  # [g-1] proportional to the slope at inflextion point of Lp
coefLp = 55 # [g] fresh mass of the inflection point

Lx0 = 0
Lxmax = 0.24  # [g.cm-2.h-1.bar-1] Conductivity for water transport
Lxmin = 0  # [g.cm-2.h-1.bar-1] Conductivity for water transport
kLx = 0.32  # [g-1] proportional to the slope at inflextion point of Lp
coefLx = 33  # [g] fresh mass of the inflection point

axp = 3.5e-3  # [adim] Coefficient for the xylem area of composite membrane
sigmaX = 1  # [adim] reflection coefficient xylem

sigmaP0 = 0.9  # [adim] reflection coefficient phloem
tauS = 0.50  # [h] variation de la reflexion membranaire avec le temps
t1 = 35

piP0 = 8.8  # [bar] phloem osmotic pressure from other solutions
piF0 = 3.5  # [bar] fruit osmotic pressure from other solutions
pF0 = 3.0  # [bar] initial phloem tugor pressure

apaa = 21.180  # 5.91978,
bpaa = 10.493 # 5.0524,
cpaa = 0.353  # 0.32911,

# Variables allometriques diverses ------------------------------------------

assrat = 0.62606  # [adim] % of solubles solids / total : pente
bssrat = 0.28124  # [adim] % of solubles solids / total : OO
cssrat = 58.203  # [adim] % of solubles solids / total : OO
dssrat = 0.043  # [adim] % of solubles solids / total : OO

rwtr = 0.6  # [adim] ratio of water coming from respiration
gamma = 4.152  # [adim] emperical parameter
eta = 0.707  # [adim] emperical parameter for the relation of fruit surface and fruit weight

share1 = 1.721e-02  # [adim] emperical parameter to calculate stone dry weight
share2 = 27.99  # [adim] emperical parameter to calculate stone dry weight

noysf1 = 0.01772239  # [adim] emperical parameter to calculate stone fresh weight
noysf2 = 0.63535422  # [adim] emperical parameter to calculate stone fresh weight

aRatio_suc = 0.35698  # [adim] emperical parameter to calculate sucrose ratio
bRatio_suc = 0.34726  # [adim] emperical parameter to calculate sucrose ratio
cRatio_suc = 46.66770  # [adim] emperical parameter to calculate sucrose ratio
dRatio_suc = 0.11757  # [adim] emperical parameter to calculate sucrose ratio

# Variables initiales, duree simulation et contraintes ----------------------

stot0 = 0.00415  # [g]	Initial dry weight of the fruit
wtot0 = 0.0118  # [g] Initial water amount in the fruit

# Variables control calculations ----------------------
elasticity = True  # [T/F] Consider elaticity or not
calcunuM = True  # [T/F] Calculate stone weight or not
calcuphi = False  # [T/F] Calculate phi or not
calcupiF0 = True  # [T/F] Calculate piF0 demande
calcuyC = False  # [T/F] Calculate yC demande
calcustone = True  # [T/F] Calculate stone weight or not
calcuel = True  # [T/F] Calculate el or not

calcuSigmaP = 3  # [numeric] Calculate sigmaP method
calcuRho = 4  # [numeric] Calculate rho demande (the input rho will be replaced by the calculating)
calcussrat = 3  # [numeric] Calculate ssrat method
calcuVmax = 3  # [numeric] Calculate Vmax method
calcuLp = 3  # [numeric] Calculate Lp or not
calcuLx = 3  # [numeric] Calculate Lx or not

# Matrice d'entree des conditions aux limites (5 colonnes) ------------------
DATA = DATA_input

mw.modgrowth_wrapper(
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