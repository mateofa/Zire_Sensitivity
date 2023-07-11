import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad
#Import Efficientcy Points

f_Sens = open('Sensitivity_Zire.dat', 'w')

#f_new = open('./EfficiencyTXT2/Efficiency_Angle_0.txt', 'r')
#f_eff = open('/Users/mateo/Documents/GSSI/Nuses/Zire/Geant4/ZIRE/ResultAnalysis_NewGeo/Efficiency/EfficiencyVSenergy_180.dat', 'r')
f_eff = open('/Users/mateo/Documents/GSSI/Nuses/Zire/Geant4/ZIRE/ResultAnalysis_NewGeo/Eff_Area/EffAreaVSenergy_60.dat', 'r')

lines_eff=f_eff.readlines()
eff=[]
eff_energy=[]
for x in lines_eff:
    if (float(x.split(' ')[3])==0):
        eff.append(float(x.split(' ')[3])+1e-1)
        eff_energy.append(float(x.split(' ')[0]))
    else:
        eff.append(float(x.split(' ')[3]))
        eff_energy.append(float(x.split(' ')[0]))
f_eff.close()



##Import simulated and weighted background count Points - FIRST multiply by 2pi (it's in sr-1 and we integrate over half sphere) to get counts

## Background
f_back = open('/Users/mateo/Documents/GSSI/Nuses/Zire/Geant4/ZIRE/ResultAnalysis_NewGeo/Backgrounds/PST_Center/TotalBackground_ZIRE_LOW.dat', 'r')
lines_back=f_back.readlines()
back_flux=[]
back_energy=[]
deltaE=[]

for x in lines_back:
    deltaE.append(float(x.split(' ')[2])-float(x.split(' ')[1]))
    #print float(x.split(' ')[2])-float(x.split(' ')[1])
    back_flux.append(2*np.pi*float(x.split(' ')[3]))
    back_energy.append(float(x.split(' ')[0]))


f_back.close()



# Define Sensitivity function and create Sensitivity points in the background energy values

Sens_value=[]

A_tot=8*8
T=3.156e+7 #1 year in seconds
#T=1e+6 #Time in seconds
deltaOmega=2*np.pi
def Sens(e,B,deltaE):
    #return  3/e * np.sqrt(B/(A_tot*T*deltaE))
    return  3/0.68 * np.sqrt(B*2*np.pi/(e*T*deltaE))

#def Sens(B,deltaE,Aeff,deltaOmega):
#    return  3/epsilon_68 * np.sqrt(B*deltaOmega/(Aeff*T*deltaE))

#Get interpolated efficiency values in background energy bins
int_eff = np.interp(back_energy,eff_energy,eff)


##Compute Sensitivityx
for i in range(np.size(back_flux)):
    Sens_value.append(back_energy[i]*back_energy[i]*Sens(int_eff[i],back_flux[i],deltaE[i]))
    print >> f_Sens, back_energy[i]," ",back_energy[i]*back_energy[i]*Sens(eff[i],back_flux[i],deltaE[i])
    #print back_energy[i]," ",eff[i]," ",back_flux[i]," ",back_energy[i]*back_energy[i]*Sens(eff[i],back_flux[i],deltaE[i])*1.6e-9


    #Sens_value.append(Sens(int_eff[i],back_int[i]))
f_Sens.close()
##Get Sensitivity values for particular wavelengths
int_back_flux_511 = np.interp(511,back_energy,back_flux) # Background @ 511 keV
int_eff_511 = np.interp(511,eff_energy,eff) # Eff @ 511 keV
int_deltaE_511 = np.interp(511,back_energy,deltaE) # Background @ 511 keV

Sens_511 =  Sens(int_eff_511,int_back_flux_511,int_deltaE_511)
#print Sens_511

##Add other experiments ###

## eAstrogam
f_astrogam = open('./eAstrogamData.dat', 'r')
lines_astrogam=f_astrogam.readlines()
Sens_astrogam=[]
energy_astrogam=[]
for x in lines_astrogam:
    Sens_astrogam.append(float(x.split(' ')[1]))
    energy_astrogam.append(1e3*float(x.split(' ')[0]))
f_astrogam.close()


## SPI
f_SPI = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/Sensitivity/NewSensitivity_Efficiency/SPIData.dat', 'r')
lines_SPI=f_SPI.readlines()
Sens_SPI=[]
energy_SPI=[]
for x in lines_SPI:
    Sens_SPI.append(float(x.split(' ')[1]))
    energy_SPI.append(1e3*float(x.split(' ')[0]))
f_SPI.close()


## IBIS_ISG
f_IBIS_ISG = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/Sensitivity/NewSensitivity_Efficiency/IBIS_ISGRI_erg_Mev.dat', 'r')
lines_IBIS_ISG=f_IBIS_ISG.readlines()
Sens_IBIS_ISG=[]
energy_IBIS_ISG=[]
for x in lines_IBIS_ISG:
    Sens_IBIS_ISG.append(float(x.split(' ')[1]))
    energy_IBIS_ISG.append(1e3*float(x.split(' ')[0]))
f_IBIS_ISG.close()

## IBIS_PIC
f_IBIS_PIC = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/Sensitivity/NewSensitivity_Efficiency/IBIS_PICs_erg_Mev.dat', 'r')
lines_IBIS_PIC=f_IBIS_PIC.readlines()
Sens_IBIS_PIC=[]
energy_IBIS_PIC=[]
for x in lines_IBIS_PIC:
    Sens_IBIS_PIC.append(float(x.split(' ')[1]))
    energy_IBIS_PIC.append(1e3*float(x.split(' ')[0]))
f_IBIS_PIC.close()

## Comptel
f_Comptel = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/Sensitivity/NewSensitivity_Efficiency/Comptel_erg_MeV.dat', 'r')
lines_Comptel=f_Comptel.readlines()
Sens_Comptel=[]
energy_Comptel=[]
for x in lines_Comptel:
    Sens_Comptel.append(float(x.split(' ')[1]))
    energy_Comptel.append(1e3*float(x.split(' ')[0]))
f_Comptel.close()

## Egret
f_Egret = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/Sensitivity/NewSensitivity_Efficiency/EGRET_erg_Mev.dat', 'r')
lines_Egret=f_Egret.readlines()
Sens_Egret=[]
energy_Egret=[]
for x in lines_Egret:
    Sens_Egret.append(float(x.split(' ')[1]))
    energy_Egret.append(1e3*float(x.split(' ')[0]))
f_Egret.close()

keVTOerg = 1.6e-9
keVTOGeV = 1e-6
##Crab Flux

def crab_flux(E):
    if(E<1e2):
        return(10.2*(E)**(-2.105)*E*E*keVTOerg)
    if(1e2<=E):
        return(10.2*(E)**(-2.22)*E*E*keVTOerg)
v_crab_flux=np.vectorize(crab_flux)


## GRB flux

A = 1.48e-2 # s-1 cm-2 between 10-1000 keV
alpha = -0.62
Epeak = 185

def GRB(x):
     return A*(x/100)**(alpha)*np.exp(-(alpha+2)*x/Epeak)*x*x*keVTOerg

ener_keV = np.linspace(10,3e4)
ener_GRB = np.linspace(10,1e3)

## Check ratio of values between CrysrtalEye and other experiments

CE_Sens_value_1MeV= np.interp(1e3,back_energy,Sens_value)*keVTOerg
AstroGam_Sens_value_1MeV= np.interp(1e3,energy_astrogam,Sens_astrogam)

ratio_CE_Astrogam = CE_Sens_value_1MeV/AstroGam_Sens_value_1MeV
#print ratio_CE_Astrogam



#Cut first point (19 to start after the drop)
Sens_value=Sens_value[5:]

back_energy=back_energy[5:]

######Sensitivity Plotting #############

fig, ax1 = plt.subplots()

#Change units
#back_energy_MeV=np.multiply(back_energy,1e-3)
#Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg=np.multiply(Sens_value,keVTOerg) #keV to erg

ax1.plot(back_energy,Sens_value_erg,label='Zire Sensitivity 1yr', linestyle='solid',linewidth=3,alpha=0.8)
#plt.plot(energy_astrogam,Sens_astrogam,label='eAstrogam Sensitivity 1yr ', linestyle='dashed')
#plt.plot(energy_SIP,Sens_SIP,label='SPI Sensitivity 1Ms ', linestyle='dashed')
plt.plot(energy_SPI,Sens_SPI,label='SPI Sensitivity 1yr ', linestyle='dashed',alpha=0.8)
plt.plot(energy_IBIS_ISG,Sens_IBIS_ISG,label='IBIS_ISGRI 1yr', linestyle='dashed',alpha=0.8)
plt.plot(energy_IBIS_PIC,Sens_IBIS_PIC,label='IBIS_PICsIT 1yr', linestyle='dashed',alpha=0.8)
plt.plot(energy_Comptel,Sens_Comptel,label='Comptel 9yr', linestyle='dashed',alpha=0.8)
#plt.plot(energy_Egret,Sens_Egret,label='Egret 9yr', linestyle='dashed',alpha=0.8)
#ax1.plot(back_albedo_energy,Sens_value_erg,label='CrystalEye Sensitivity 1yr', linestyle='dashed')

plt.plot(ener_keV,v_crab_flux(ener_keV),label='1 Crab', linestyle='dotted',alpha=0.8)
plt.plot(ener_keV,v_crab_flux(ener_keV)*1e-1,label='100 mCrab', linestyle='dotted',alpha=0.8)
plt.plot(ener_keV,v_crab_flux(ener_keV)*1e-3,label='1 mCrab', linestyle='dotted',alpha=0.8)

#plt.plot(ener_GRB,GRB(ener_GRB),label='GRB 170817A (Goldstein 2017)', linestyle='dotted')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('E [keV]')
ax1.set_ylabel(r' 3$\sigma$ Sensitivity [erg cm$^{-2}$ s$^{-1}$]')
ax1.legend(loc='lower right', prop={'size': 10})
plt.title('Zire Continuum Sensitivity (low background)')
plt.xlim([10,8e4 ])
plt.ylim([1e-12,2e-8 ])

ax1.text(0.05, 0.1, 'Preliminary', transform=ax1.transAxes,
        fontsize=20, color='gray', alpha=0.5,
        ha='left', va='top', rotation=0)
#ax2 = ax1.twinx()
#ax2.plot(eff_energy,eff,label='Efficiency', linestyle='dashed',color='red')
#ax2.set_yscale('log')
#ax2.set_xscale('log')
#ax2.set_xlabel('E [keV]')
#ax2.set_ylabel('Efficiency')
#ax2.legend()

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.savefig('Sensitivity_fromEffArea_LOW.png')
