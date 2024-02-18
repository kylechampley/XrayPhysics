import sys
import numpy as np
import matplotlib.pyplot as plt

from xrayphysics import *
physics = xrayPhysics()

# The next line prints out the version number and some other information about this package
#physics.about()

whichPlot = 2


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
# Example 2: Calculate cross section of bone by mass fractions
#########################################################################################
if whichPlot == 2:
    gammas = np.array(range(20,100),dtype=np.float32)+1.0

    # Mass Density of Bone (g/cm^3)    
    massDensity = 1.85
    
    # Here is how to calculate the cross section of bone explicitly
    sigma_bone = gammas.copy()
    sigma_bone[:] = 0.0
    sigma_bone += 0.063984*physics.sigma(1, gammas)
    sigma_bone += 0.278000*physics.sigma(6, gammas)
    sigma_bone += 0.027000*physics.sigma(7, gammas)
    sigma_bone += 0.410016*physics.sigma(8, gammas)
    sigma_bone += 0.002000*physics.sigma(12, gammas)
    sigma_bone += 0.070000*physics.sigma(15, gammas)
    sigma_bone += 0.002000*physics.sigma(16, gammas)
    sigma_bone += 0.147000*physics.sigma(20, gammas)

    # Here is how to calculate the cross section of bone using a string (which I think is more convenient)
    sigma_bone_2 = physics.sigma('0.063984*H+0.278000*C+0.027000*N+0.410016*O+0.002000*Mg+0.070000*P+0.002000*S+0.147000*Ca', gammas)
    
    # Here is another example for salt water where the salt is 17% of the total mass
    # Note that the scaling of the mass fractions does not matter, the XrayPhysics package always normalizes by the sum of all the mass fractions
    #sigma_salt_water = physics.sigma('17*NaCl + 83*H2O', gammas)
    
    plt.plot(gammas, sigma_bone*massDensity, gammas, sigma_bone_2*massDensity)
    plt.title('Bone LAC')
    plt.xlabel('x-ray energy (keV)')
    plt.ylabel('linear attenuation coefficient (cm^-1)')
    plt.show()


#########################################################################################
# Example 3: Modeling a total system spectral response of the x-ray system
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

if whichPlot == 3:
    plt.plot(Es, s_total/np.sum(s_total), 'k-')
    plt.title('Total System Spectral Response')
    plt.xlabel('x-ray energy (keV)')
    plt.ylabel('normalized response (unitless)')
    plt.show()


#########################################################################################
# Example 4: Effective Energy Through a Material of Variable Thickness
#########################################################################################
thicknesses = np.array(range(26))/10.0
e_effs = thicknesses.copy()
for n in range(thicknesses.size):
    e_effs[n] = physics.effectiveEnergy('H2O', 1.0, thicknesses[n], s_total, Es)
if whichPlot == 4:
    plt.plot(thicknesses, e_effs, 'k-')
    plt.title('Effective Energy (keV)')
    plt.xlabel('thickness (cm)')
    plt.ylabel('LAC (cm^-1)')
    plt.show()


#########################################################################################
# Example 5: Generate BH Lookup Table
#########################################################################################
LUT, T_lut = physics.setBHlookupTable(s_total, Es, 'Al')
monoAttens = np.array(range(LUT.size))*T_lut
if whichPlot == 5:
    plt.plot(monoAttens, LUT, 'k-')
    plt.title('Beam Hardening Transfer Function')
    plt.xlabel('monochromatic attenuation (unitless)')
    plt.ylabel('polychromatic attenuation (unitless)')
    plt.show()
    

#########################################################################################
# Example 6: Generate BHC Lookup Table
#########################################################################################
import time
LUT, T_lut = physics.setBHClookupTable(s_total, Es, 'Al')
polyAttens = np.array(range(LUT.size))*T_lut
if whichPlot == 6:
    plt.plot(polyAttens, LUT, 'k-')
    plt.title('Beam Hardening Correction Transfer Function')
    plt.xlabel('polychromatic attenuation (unitless)')
    plt.ylabel('monochromatic attenuation (unitless)')
    plt.show()


#########################################################################################
# Example 7: Generate Polynomial BHC Coefficients
#########################################################################################
coeff = physics.polynomialBHC(s_total, Es, 'Al', 2.7, referenceEnergy=0.0, maxThickness=10.0, order=2)
print(coeff)


#########################################################################################
# Example 8: Calculate the effective atomic number of a compound
#########################################################################################
chemForm = 'H2O'
Ze = physics.effectiveZ(chemForm)
print('effective atomic number of ' + str(chemForm) + ' is ' + str(Ze))


#########################################################################################
# Example 9: Calculate scattering distributions
#########################################################################################
thetas = np.array(range(180+1),dtype=np.float32)
dsigma_incoh = physics.incoherentScatterDistribution('H2O', 60.0, thetas, doNormalize=True)
dsigma_coh = physics.coherentScatterDistribution('H2O', 60.0, thetas, doNormalize=True)
# Setting doNormalize=True makes these so that they act like probability distributions over [0, 180]
    
if whichPlot == 9:
    plt.figure(figsize=(10,5))
    plt.subplot(1, 2, 1)
    plt.plot(thetas, dsigma_incoh)
    plt.title('Incoherent (Compton) Scatter Distributions')
    plt.xlabel('angle (degrees)')
    plt.ylabel('normalized differential cross section')
    plt.xlim((0.0,180.0))
    
    plt.subplot(1, 2, 2)
    plt.plot(thetas, dsigma_coh)
    plt.title('Coherent (Rayleigh) Scatter Distributions')
    plt.xlabel('angle (degrees)')
    plt.ylabel('normalized differential cross section')
    plt.xlim((0.0,180.0))
    plt.show()
