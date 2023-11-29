import sys
import numpy as np
import matplotlib.pyplot as plt

from xrayphysics import *
physics = xrayPhysics()


#########################################################################################
# Example 1: Getting linear attenuation coefficients for a compound
#########################################################################################
# Set a numpy array of the x-ray energies
gammas = np.array(range(20,100),dtype=np.float32)+1.0

# Set the LAC (mm^-1) of water and its components
mu_water = physics.mu('H2O',gammas,1.0)
mu_water_PE = physics.muPE('H2O',gammas,1.0)
mu_water_CS = physics.muCS('H2O',gammas,1.0)
mu_water_RS = physics.muRS('H2O',gammas,1.0)

# Plot results
plt.plot(gammas,mu_water,gammas,mu_water_PE,gammas,mu_water_CS,gammas,mu_water_RS)
plt.title('Water LAC')
plt.xlabel('x-ray energy (keV)')
plt.ylabel('linear attenuation coefficient (mm^-1)')
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

plt.plot(Es, s_total/np.sum(s_total), 'k-')
plt.title('Total System Spectral Response')
plt.xlabel('x-ray energy (keV)')
plt.ylabel('normalized response (unitless)')
plt.show()
