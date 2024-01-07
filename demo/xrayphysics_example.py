import sys
import numpy as np
import matplotlib.pyplot as plt

from xrayphysics import *
physics = xrayPhysics()

whichPlot = 4


#########################################################################################
# Example 1: Getting linear attenuation coefficients for a compound
#########################################################################################
# Set a numpy array of the x-ray energies
gammas = np.array(range(20,100),dtype=np.float32)+1.0

# Set the LAC (cm^-1) of water and its components
chemForm = 'H2O'
mu_water = physics.mu(chemForm,gammas,1.0)
mu_water_PE = physics.muPE(chemForm,gammas,1.0)
mu_water_CS = physics.muCS(chemForm,gammas,1.0)
mu_water_RS = physics.muRS(chemForm,gammas,1.0)

# Plot results
if whichPlot == 1:
    plt.plot(gammas,mu_water,gammas,mu_water_PE,gammas,mu_water_CS,gammas,mu_water_RS)
    plt.title('Water LAC')
    plt.xlabel('x-ray energy (keV)')
    plt.ylabel('linear attenuation coefficient (cm^-1)')
    plt.legend(['total','Photoelectric','Compton Scatter','Rayleigh Scatter'])
    plt.show()


#########################################################################################
# Example 2: Modeling a total system spectral response of the x-ray system
#########################################################################################
# Define the kV of the source voltage and the take-off angle (degrees)
kV = 100.0
takeOffAngle = 11.0

# First simulate the source spectrum (units are photons/(bin * mAs * sr))
Es, s = physics.simulateSpectra(kV,takeOffAngle)

# Then model the detector response as the product of the
# x-ray energy and the stopping power of the scintillator
detResp = physics.detectorResponse('O2SGd2', 7.32, 0.01, Es)

# Finally model the attenuation due to the filters
filtResp = physics.filterResponse('Al', 2.7, 0.1, Es)

# Take the product of all three factors
s_total = s*filtResp*detResp

if whichPlot == 2:
    plt.plot(Es, s_total/np.sum(s_total), 'k-')
    plt.title('Total System Spectral Response')
    plt.xlabel('x-ray energy (keV)')
    plt.ylabel('normalized response (unitless)')
    plt.show()


#########################################################################################
# Example 3: Effective Energy Through a Material of Variable Thickness
#########################################################################################
thicknesses = np.array(range(26))/10.0
e_effs = thicknesses.copy()
for n in range(thicknesses.size):
    e_effs[n] = physics.effectiveEnergy('H2O', 1.0, thicknesses[n], s_total, Es)
if whichPlot == 3:
    plt.plot(thicknesses, e_effs, 'k-')
    plt.title('Effective Energy (keV)')
    plt.xlabel('thickness (cm)')
    plt.ylabel('LAC (cm^-1)')
    plt.show()


#########################################################################################
# Example 4: Generate BH Lookup Table
#########################################################################################
LUT, T_lut = physics.setBHlookupTable('Al', s_total, Es)
monoAttens = np.array(range(LUT.size))*T_lut
if whichPlot == 4:
    plt.plot(monoAttens, LUT, 'k-')
    plt.title('Beam Hardening Transfer Function')
    plt.xlabel('monochromatic attenuation (unitless)')
    plt.ylabel('polychromatic attenuation (unitless)')
    plt.show()
    

#########################################################################################
# Example 5: Generate BHC Lookup Table
#########################################################################################
import time
LUT, T_lut = physics.setBHClookupTable('Al', s_total, Es)
polyAttens = np.array(range(LUT.size))*T_lut
if whichPlot == 5:
    plt.plot(polyAttens, LUT, 'k-')
    plt.title('Beam Hardening Correction Transfer Function')
    plt.xlabel('polychromatic attenuation (unitless)')
    plt.ylabel('monochromatic attenuation (unitless)')
    plt.show()


#########################################################################################
# Example 6: Generate Polynomial BHC Coefficients
#########################################################################################
coeff = physics.polynomialBHC('Al', 2.7, s_total, Es, referenceEnergy=0.0, maxThickness=10.0, order=2)
print(coeff)


#########################################################################################
# Example 7: Calculate the effective atomic number of a compound
#########################################################################################
chemForm = 'H2O'
Ze = physics.effectiveZ(chemForm)
print('effective energy of ' + str(chemForm) + ' is ' + str(Ze))
