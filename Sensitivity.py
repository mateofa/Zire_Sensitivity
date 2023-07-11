import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad

#Import Efficientcy Points

f = open('./Efficiency_Angle180.txt', 'r')
lines=f.readlines()
eff=[]
eff_energy=[]
for x in lines:
    eff.append(float(x.split('  ')[1]))
    eff_energy.append(float(x.split('  ')[0]))
f.close()

#Import eAstrogam background points

f_astroGam = open('Back_eAstrogam.dat', 'r')
lines_astroGam=f_astroGam.readlines()
energy_astroGam=[]
back_astroGam=[]
for x in lines_astroGam:
    energy_astroGam.append(1e3*float(x.split(' ')[0]))
    back_astroGam.append(1e-7*float(x.split(' ')[1])/(2*np.pi))
f_astroGam.close()

# Define Background functions

def albedo_back(E):
    C_albedo = 1.48e-2 #Taken from Ajello 2008
    return(C_albedo/((E/33.7)**(-5)+(E/33.7)**(1.72)))

def cosmic_back(E):
    C_cosmic = 10.15e-2 #Taken from Ajello 2008
    return(C_cosmic/((E/30)**(1.32)+(E/30)**(2.88)))


def back_total(E):
    return cosmic_back(E) + albedo_back(E)


#Define binning and integrate in each binning
back_int=[]
bin_center=[]

deltaE=10 #Bin width in keV
E_min=1
E_max=5e4
int_limits=np.linspace(E_min,E_max,(E_max-E_min)/deltaE)


for i in range(np.size(int_limits)-1):
    back_int.append(2*np.pi*quad(back_total,int_limits[i],int_limits[i+1])[0])
    bin_center.append(int_limits[i]+(int_limits[i+1]-int_limits[i])/2)


#Get interpolated efficiency values in bin_centers
int_eff = np.interp(bin_center,eff_energy,eff)

# Define Sensitivity function and create Sensitivity points in the bin_centers
Sens_value=[]
Sens_value2=[]

#A_geo=1321 #Area of Crystal CrystalEye
A_det1=12*12 #Area of bottom window Zire
A_det2=12*10 #Area of side window Zire
A_geo=A_det1+A_det2
T=5*3.156e+7 #1 year in seconds
#T=1e6 #1 M seconds

def Sens(e,B):
    return 3/(e*deltaE) * np.sqrt(B/(A_geo*T)) #in this case B is in ph cm-2 s-1

def Sens2(e,B):
    return 3/e * np.sqrt(B/(A_geo*T*deltaE)) #in this case B is in ph cm-2 s-1 keV-1


for i in range(np.size(bin_center)):
    Sens_value.append(bin_center[i]*bin_center[i]*Sens(int_eff[i],back_int[i]))
    Sens_value2.append(bin_center[i]*bin_center[i]*Sens2(int_eff[i],2*np.pi*back_total(bin_center[i])))


## eAstrogam Sensitivity
f_astrogam = open('./eAstrogamData.dat', 'r')
lines_astrogam=f_astrogam.readlines()
Sens_astrogam=[]
energy_astrogam=[]
for x in lines_astrogam:
        Sens_astrogam.append(float(x.split(' ')[1]))
        energy_astrogam.append(float(x.split(' ')[0]))
f_astrogam.close()


## SPI Sensitivity
f_SPI = open('./SPIData.dat', 'r')
lines_SPI=f_SPI.readlines()
Sens_SPI=[]
energy_SPI=[]
for x in lines_SPI:
        Sens_SPI.append(float(x.split(' ')[1]))
        energy_SPI.append(float(x.split(' ')[0]))
f_SPI.close()


## eCOMPTEL Sensitivity
f_COMPTEL = open('./COMPTELData.dat', 'r')
lines_COMPTEL=f_COMPTEL.readlines()
Sens_COMPTEL=[]
energy_COMPTEL=[]
for x in lines_COMPTEL:
        Sens_COMPTEL.append(float(x.split(' ')[1]))
        energy_COMPTEL.append(float(x.split(' ')[0]))
f_COMPTEL.close()
###Plotting######

######Sensitivity Plottting #############

#Change units
bin_center_MeV=np.multiply(bin_center,1e-3)
Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg=np.multiply(Sens_value,1.6e-9)
Sens2_value_erg=np.multiply(Sens_value2,1.6e-9)


#plt.title("Sensitivity 1 yr")
#plt.plot(bin_center_MeV,Sens_value_erg,label='Zire', linestyle='dashed')
plt.plot(bin_center_MeV,Sens2_value_erg,label='ZIRE (5 yr)', linestyle='solid', linewidth=1.5)

plt.plot(energy_astrogam,Sens_astrogam,label='eAstrogam (1 yr)', linestyle='dashed')
plt.plot(energy_SPI,Sens_SPI,label='SPI (1M s)', linestyle='dashed')
plt.plot(energy_COMPTEL,Sens_COMPTEL,label='COMPTEL (9 yr)', linestyle='dashed')


plt.yscale('log')
plt.xscale('log')
plt.xlabel('E [MeV]')
plt.ylabel(r' 3$\sigma$ Flux [erg cm$^{-2}$ s$^{-1}$]')
plt.legend()

plt.savefig('Sensitivity.png')

###Backgound plotting###########

#Define energy points and get log-flux

E_full = np.logspace(1,4,100)
albedo_flux = albedo_back(E_full)
cosmic_flux =cosmic_back(E_full)
total_flux =back_total(E_full)

#log_E = np.log10(E_full)
#log_albedo_flux = np.log10(albedo_back(E_full))
#log_cosmic_flux = np.log10(cosmic_back(E_full))

plt.plot(energy_astroGam,back_astroGam,label='Cosmic eAstrogam', linestyle='dashed')

plt.plot(E_full,albedo_flux,label='albedo', linestyle='dashed')
plt.plot(E_full,cosmic_flux,label='cosmic', linestyle='dashed')
plt.plot(E_full,total_flux,label='total')
#plt.plot(bin_center,integral,label='integral')
plt.legend()
plt.xlabel('E [keV]')
plt.ylabel(r'Flux [ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$ sr$^{-1}$]')

plt.yscale('log')
plt.xscale('log')

plt.savefig('Flux_Back.png')


##Efficiency Plotting
fig = plt.figure()
#plt.yscale('log')
plt.xscale('log')
plt.title('Efficiency')
plt.xlabel('E [keV]')
plt.ylabel(r'N$_{cuts}$/N$_{FoV}$')
plt.plot(eff_energy,eff,label='Efficiency (180$^o$)', linestyle='dashed')
plt.savefig('Efficiency_Zire.png')
