import ctypes
import os
import sys
from sys import platform as _platform
from numpy.ctypeslib import ndpointer
import numpy as np


class xrayPhysics:

    def __init__(self, lib_dir=""):
        if len(lib_dir) > 0:
            current_dir = lib_dir
        else:
            current_dir = os.path.abspath(os.path.dirname(__file__))

        if _platform == "linux" or _platform == "linux2":
            import readline
            from ctypes import cdll
            self.libxrayphysics = cdll.LoadLibrary(os.path.join(current_dir, "../build/lib/libxrayphysics.so"))
        elif _platform == "win32":
            from ctypes import windll
            self.libxrayphysics = windll.LoadLibrary(os.path.join(current_dir, r'..\win_build\bin\Release\libxrayphysics.dll'))
        elif _platform == "darwin":  # Darwin is the name for MacOS in Python's platform module
            from ctypes import cdll
            # Adjust the path to where your .dylib file is located
            self.libxrayphysics = cdll.LoadLibrary(os.path.join(current_dir, "../build/lib/libxrayphysics.dylib"))

    def mu(self, Z, gamma, massDensity):
        return massDensity * self.sigma(Z, gamma)
        
    def muPE(self, Z, gamma, massDensity):
        return massDensity * self.sigmaPE(Z, gamma)
        
    def muCS(self, Z, gamma, massDensity):
        return massDensity * self.sigmaCS(Z, gamma)
        
    def muRS(self, Z, gamma, massDensity):
        return massDensity * self.sigmaRS(Z, gamma)
        
    def muPP(self, Z, gamma, massDensity):
        return massDensity * self.sigmaPP(Z, gamma)
        
    def muTP(self, Z, gamma, massDensity):
        return massDensity * self.sigmaTP(Z, gamma)

    def sigma(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaCompound.restype = ctypes.c_float
            self.libxrayphysics.sigmaCompound.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCompound(chemicalFormula, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCompound(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigma.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigma.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigma(Z, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigma(Z, gamma)
        
    def sigmaPE(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaCompoundPE.restype = ctypes.c_float
            self.libxrayphysics.sigmaCompoundPE.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCompoundPE(chemicalFormula, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCompoundPE(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaPE.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaPE.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaPE(Z, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaPE(Z, gamma)
            
    def sigmaCS(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaCompoundCS.restype = ctypes.c_float
            self.libxrayphysics.sigmaCompoundCS.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCompoundCS(chemicalFormula, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCompoundCS(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaCS.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaCS.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCS(Z, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCS(Z, gamma)
            
    def sigmaRS(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaCompoundRS.restype = ctypes.c_float
            self.libxrayphysics.sigmaCompoundRS.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCompoundRS(chemicalFormula, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCompoundRS(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaRS.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaRS.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaRS(Z, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaRS(Z, gamma)
            
    def sigmaPP(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaCompoundPP.restype = ctypes.c_float
            self.libxrayphysics.sigmaCompoundPP.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCompoundPP(chemicalFormula, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCompoundPP(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaPP.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaPP.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaPP(Z, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaPP(Z, gamma)
            
    def sigmaTP(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaCompoundTP.restype = ctypes.c_float
            self.libxrayphysics.sigmaCompoundTP.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCompoundTP(chemicalFormula, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaCompoundTP(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaTP.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaTP.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaTP(Z, float(gamma[n]))
                return retVal
            else:
                return self.libxrayphysics.sigmaTP(Z, gamma)
    
    def simulateSpectra(self, kV, takeOffAngle=11.0, Z=74, gammas=None):
        #simulateSpectra(float kV, float takeOffAngle, int Z, float* gammas, int N, float* output)
        
        if gammas is None:
            maxEnergy = max(1,int(np.ceil(kV)))
            minEnergy = max(1,int(0.1*np.ceil(kV)))
            gammas = np.ascontiguousarray(np.array(range(minEnergy,maxEnergy+1)), dtype=np.float32)
        
        sourceSpectrum = np.ascontiguousarray(np.zeros(gammas.size), dtype=np.float32)
        
        self.libxrayphysics.simulateSpectra.restype = ctypes.c_bool
        self.libxrayphysics.simulateSpectra.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]
        self.libxrayphysics.simulateSpectra(kV, takeOffAngle, Z, gammas, gammas.size, sourceSpectrum)
        return gammas, sourceSpectrum
        
    def changeTakeOffAngle(self, kVp, takeOffAngle_cur, takeOffAngle_new, Z, gammas, spectrum_cur):
        #bool changeTakeOffAngle(float kVp, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s)
        self.libxrayphysics.changeTakeOffAngle.restype = ctypes.c_bool
        self.libxrayphysics.changeTakeOffAngle.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]
        spectrum_new = spectrum_cur.copy()
        self.libxrayphysics.changeTakeOffAngle(kVp, takeOffAngle_cur, takeOffAngle_new, Z, gammas, gammas.size, spectrum_new)
        return spectrum_new
        
    def detectorResponse(self, chemicalFormula, density, thickness, gammas):
        detResp = np.ascontiguousarray(np.zeros(gammas.size), dtype=np.float32)
        for n in range(gammas.size):
            detResp[n] = gammas[n]*(1.0-np.exp(-self.mu(chemicalFormula, gammas[n], density)*thickness))
        return detResp
        
    def filterResponse(self, chemicalFormula, density, thickness, gammas):
        filtResp = np.ascontiguousarray(np.zeros(gammas.size), dtype=np.float32)
        for n in range(gammas.size):
            filtResp[n] = np.exp(-self.mu(chemicalFormula, gammas[n], density)*thickness)
        return filtResp
        
    def meanEnergy(self, spectralResponse, gammas):
        self.libxrayphysics.meanEnergy.restype = ctypes.c_float
        self.libxrayphysics.meanEnergy.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
        return self.libxrayphysics.meanEnergy(spectralResponse, gammas, gammas.size)
        
    def effectiveAttenuation(self, Z, density, thickness, spectralResponse, gammas):
        #float effectiveAttenuation(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.effectiveAttenuation_compound.restype = ctypes.c_float
            self.libxrayphysics.effectiveAttenuation_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveAttenuation_compound(chemicalFormula, density, thickness, spectralResponse, gammas, gammas.size)
        else:
            self.libxrayphysics.effectiveAttenuation.restype = ctypes.c_float
            self.libxrayphysics.effectiveAttenuation.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveAttenuation(Z, density, thickness, spectralResponse, gammas, gammas.size)
    
    def effectiveEnergy(self, Z, density, thickness, spectralResponse, gammas):
        #float effectiveEnergy(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.effectiveEnergy_compound.restype = ctypes.c_float
            self.libxrayphysics.effectiveEnergy_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveEnergy_compound(chemicalFormula, density, thickness, spectralResponse, gammas, gammas.size)
        else:
            self.libxrayphysics.effectiveEnergy.restype = ctypes.c_float
            self.libxrayphysics.effectiveEnergy.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveEnergy(Z, density, thickness, spectralResponse, gammas, gammas.size)
            
    def transmission(self, Z, density, thickness, spectralResponse, gammas):
        #float transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.transmission_compound.restype = ctypes.c_float
            self.libxrayphysics.transmission_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.transmission_compound(chemicalFormula, density, thickness, spectralResponse, gammas, gammas.size)
        else:
            self.libxrayphysics.transmission.restype = ctypes.c_float
            self.libxrayphysics.transmission.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.transmission(Z, density, thickness, spectralResponse, gammas, gammas.size)
            
    def setBHClookupTable(self, Z, spectralResponse, gammas, T_lac=0.0, N_lac=0, referenceEnergy=0.0):
        #bool setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
        
        if N_lac <= 0 or T_lac <= 0.0:
            max_lac = 12.0
            #T_lac = 2.0e-4 # about 0.06MB; 14 counts on a 16-bit detector
            T_lac = 1.0e-4 # about 0.12MB
            #T_lac = 1.0e-5 # about 1.2MB
            #T_lac = 1.0e-6 # about 12MB
            N_lac = int(np.ceil(max_lac / T_lac)) + 1
            max_lac = float(N_lac-1)*T_lac
        
        N_gamma = gammas.size
        if referenceEnergy <= 0.0:
            #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
            referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spectralResponse, gammas)
            print('referenceEnergy = ' + str(referenceEnergy))
        LUT = np.zeros(N_lac, dtype=np.float32)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.setBHClookupTable_compound.restype = ctypes.c_bool
            self.libxrayphysics.setBHClookupTable_compound.argtypes = [ctypes.c_char_p, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int, ctypes.c_float]
            self.libxrayphysics.setBHClookupTable_compound(chemicalFormula, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy)
        else:
            self.libxrayphysics.setBHClookupTable.restype = ctypes.c_bool
            self.libxrayphysics.setBHClookupTable.argtypes = [ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int, ctypes.c_float]
            self.libxrayphysics.setBHClookupTable(Z, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy)
        return LUT, T_lac
        
    def polynomialBHC(self, Z, density, spectralResponse, gammas, referenceEnergy=0.0, maxThickness=10.0, order=2):
        order = max(1, min(order, 10))
        if maxThickness <= 0.0:
            return None
        if referenceEnergy <= 0.0:
            #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
            referenceEnergy = self.effectiveEnergy(Z, density, 0.0, spectralResponse, gammas)
            print('referenceEnergy = ' + str(referenceEnergy))
            
        L = (np.array(range(100))+1.0)*maxThickness/100.0
        A = np.zeros((order,order))
        b = np.zeros((order,1))
        polyAttenPowers = np.zeros((2*order+1,1))
        for n in range(L.size):
            curThickness = L[n]
            monoAtten = self.mu(Z, referenceEnergy, density)*curThickness
            polyAtten = -np.log(self.transmission(Z, density, curThickness, spectralResponse, gammas))
            
            #print(str(curThickness) + ': mono, poly = ' + str(monoAtten) + ', ' + str(polyAtten))
            
            polyAttenPowers = polyAtten**np.array(range(polyAttenPowers.size))
            
            for i in range(order):
                for j in range(order):
                    A[i,j] += polyAttenPowers[i+j+2]
                b[i] += monoAtten * polyAttenPowers[i+1]
        return np.concatenate((np.zeros((1,1)),np.matmul(np.linalg.inv(A), b)))
        