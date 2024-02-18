
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
from xrayphysics import *
physics = xrayPhysics()

'''
This script demonstrates how to generate a dual energy decomposition lookup table.
This table can be used to perform dual energy decomposition of dual energy data.
LEAP has examples of how these tables are used.  For examples, see here: https://github.com/LLNL/LEAP/tree/main/demo_leapctype
The way this works is one must create models of the full system spectral response
of a pair of spectra, one chooses two energy basis functions, and two energies we call "reference energies".
Let us call these spectra "low" and "high" energy.
Then the tables are generated that map the polychromatic low and high attenuation
data to monochromatic low and high attenuation data at the two reference energies.
If one wishes to recover the actual basis coefficients, one can do this with a simple
linear transformation.
These transformations are described in the following paper:
https://ieeexplore.ieee.org/abstract/document/8638824?casa_token=K_9cFGKJGvMAAAAA:EzTpZfY0qJHMvdxGniguZBS_dATpx-4vqhsDPZwB1VFh02loJFD0hvizr5RNKj5z5xgvU8Iq8g
'''

# Define the kV of the source voltage and the take-off angle (degrees)
kV_L = 100.0
kV_H = 160.0
takeOffAngle = 11.0

# First simulate the source spectrum (units are photons/(bin * mAs * sr))
# Make sure to define both spectra on the same energies
# Note that "Es" is an array of the energies at which the spectra are defined
Es, s_H = physics.simulateSpectra(kV_H, takeOffAngle)
Es, s_L = physics.simulateSpectra(kV_L, takeOffAngle, None, Es)

# Then model the detector response as the product of the
# x-ray energy and the stopping power of the scintillator
# Here we assume the detector is 0.01 cm thick GOS
detResp = physics.detectorResponse('O2SGd2', 7.32, 0.01, Es)

# Finally model the attenuation due to the filters
# The low energy uses 1 mm of aluminum and the high energy uses
# 1 mm of aluminum and 1 mm of copper
filtResp_L = physics.filterResponse('Al', 2.698899, 0.1, Es)
filtResp_H = physics.filterResponse('Cu',  8.959999, 0.1, Es)
filtResp_H[:] = filtResp_H[:]*filtResp_L[:]

# Take the product of all three factors to get the full system spectral response
s_L = s_L*filtResp_L*detResp
s_H = s_H*filtResp_H*detResp

# Now we define the basis functions.  The two most common energy basis functions
# are to either use the cross section of two materials or use the so-called
# "Compton-Photoelectric" basis.
# Use the next two lines for a material basis
sigma_water = physics.sigma('H2O', Es)
sigma_Al = physics.sigma('Al', Es)
# Or use the next two lines for the Compton-Photoelectric basis
PhotoBasisFcn = physics.PhotoelectricBasis(Es)
ComptonBasisFcn = physics.ComptonBasis(Es)
# Or use the next two lines for the PCA basis of the given materials
b_1, b_2 = physics.PCAbases(['C','N','O','Al'], Es)


# Now let's choose the reference energies of the decomposition.  The user is free to choose what they want,
# but we recommend that these energies be within the energy range of the spectra and the low energy
# reference energy be lower than the high energy reference energy.
# Here we will just use the mean energies of each spectra as the reference energy.
# This is what would be chosen if these parameters were not specified (i.e., the default value)
# Although not necessary, we round these values to the nearest whole number
referenceEnergy_L = np.round(physics.meanEnergy(s_L, Es))
referenceEnergy_H = np.round(physics.meanEnergy(s_H, Es))

# Now we generate the lookup tables.  This should take less than 4 seconds.
# Note that if s_L doesn't have a lower mean energy than s_H, this function may not work or may run slower.
# We recommend that you specify the lower energy spectra as s_L and the high energy spectra as s_H
# The return values LUT and T_atten are the lookuptables and the sampling rate of the attenuation, respectively.
# The lookup table is 3 dimensional.  The first dimension is the lookup table for the low energy channel,
# the second dimension is the lookup table for the high energy channel, and
# the third dimension is the L^2 error of the transformation (just to verify that everything worked correctly).
startTime = time.time()
LUT, T_atten = physics.setDEDlookupTable(s_L, s_H, Es, sigma_water, sigma_Al, [referenceEnergy_L, referenceEnergy_H])
#LUT,T_atten = physics.setDEDlookupTable(s_L, s_H, Es, PhotoBasisFcn, ComptonBasisFcn, [referenceEnergy_L, referenceEnergy_H])
#LUT,T_atten = physics.setDEDlookupTable(s_L, s_H, Es, b_1, b_2, [referenceEnergy_L, referenceEnergy_H])
print('DED LUT generation time: ' + str(time.time()-startTime) + ' seconds')

# Plot Results
plt.subplot(1, 3, 1)
plt.imshow(np.squeeze(LUT[0,:,:]), extent=[0.0, T_atten*(LUT.shape[1]-1), T_atten*(LUT.shape[1]-1), 0.0])
plt.xlabel('polychromatic high energy attenuation')
plt.ylabel('polychromatic low energy attenuation')
plt.title('Low Energy Synthesized Monochromatic Attenuation')

plt.subplot(1, 3, 2)
plt.imshow(np.squeeze(LUT[1,:,:]), extent=[0.0, T_atten*(LUT.shape[1]-1), T_atten*(LUT.shape[1]-1), 0.0])
plt.xlabel('polychromatic high energy attenuation')
plt.ylabel('polychromatic low energy attenuation')
plt.title('High Energy Synthesized Monochromatic Attenuation')

plt.subplot(1, 3, 3)
plt.imshow(np.squeeze(LUT[2,:,:]), extent=[0.0, T_atten*(LUT.shape[1]-1), T_atten*(LUT.shape[1]-1), 0.0])
plt.xlabel('polychromatic high energy attenuation')
plt.ylabel('polychromatic low energy attenuation')
plt.title('Decomposition L^2 Error')
plt.show()
