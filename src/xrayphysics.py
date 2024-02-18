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

        #'''
        if _platform == "linux" or _platform == "linux2":
            import readline
            from ctypes import cdll
            
            fullPath = os.path.join(current_dir, 'libxrayphysics.so')
            fullPath_backup = os.path.join(current_dir, '../build/lib/libxrayphysics.so')
            
            if os.path.isfile(fullPath):
                self.libxrayphysics = cdll.LoadLibrary(fullPath)
            elif os.path.isfile(fullPath_backup):
                self.libxrayphysics = cdll.LoadLibrary(fullPath_backup)
            else:
                print('Error: could not find XrayPhysics dynamic library at')
                print(fullPath)
                print('or')
                print(fullPath_backup)
                self.libxrayphysics = None
            
        elif _platform == "win32":
            from ctypes import windll
        
            fullPath = os.path.join(current_dir, 'libxrayphysics.dll')
            fullPath_backup = os.path.join(current_dir, r'..\win_build\bin\Release\libxrayphysics.dll')
        
            if os.path.isfile(fullPath):
                try:
                    self.libxrayphysics = windll.LoadLibrary(fullPath)
                except:
                    self.libxrayphysics = ctypes.CDLL(fullPath, winmode=0)
            elif os.path.isfile(fullPath_backup):
                try:
                    self.libxrayphysics = windll.LoadLibrary(fullPath_backup)
                except:
                    self.libxrayphysics = ctypes.CDLL(fullPath_backup, winmode=0)
            else:
                print('Error: could not find XrayPhysics dynamic library at')
                print(fullPath)
                print('or')
                print(fullPath_backup)
                self.libxrayphysics = None
        
        elif _platform == "darwin":  # Darwin is the name for MacOS in Python's platform module
            # there is current no support for XrayPhysics on Mac, but maybe someone can figure this out
            from ctypes import cdll
            
            fullPath = os.path.join(current_dir, 'libxrayphysics.dylib')
            fullPath_backup = os.path.join(current_dir, '../build/lib/libxrayphysics.dylib')
            
            if os.path.isfile(fullPath):
                self.libxrayphysics = cdll.LoadLibrary(fullPath)
            elif os.path.isfile(fullPath_backup):
                self.libxrayphysics = cdll.LoadLibrary(fullPath_backup)
            else:
                print('Error: could not find XrayPhysics dynamic library at')
                print(fullPath)
                print('or')
                print(fullPath_backup)
                self.libxrayphysics = None
        #'''

        '''
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
        #'''
            
        self.detectorBits = 16
        self.length_units = 'cm'

    def about(self):
        """prints info about this package, including the version number"""
        return self.libxrayphysics.about()

    def use_cm(self):
        """Call this function so that all units that use lengths are cm-based (this is the default setting)"""
        self.length_units = 'cm'
        
    def use_mm(self):
        """Call this function so that all units that use lengths are mm-based"""
        self.length_units = 'mm'

    def cross_section_scalar(self):
        """Returns the scalar to convert the units of cross sections
        Should only be used to alter a return values from self.libxrayphysics
        or to do inverse scaling to an argument to a self.libxrayphysics function
        """
        if self.length_units == 'mm':
            return 1.0e2
        else:
            return 1.0
            
    def LAC_scalar(self):
        """Returns the scalar to convert the units of Linear Attenuation Coefficient (LAC)
        Should only be used to alter a return values from self.libxrayphysics
        or to do inverse scaling to an argument to a self.libxrayphysics function
        """
        if self.length_units == 'mm':
            return 0.1
        else:
            return 1.0
            
    def density_scalar(self):
        """Returns the scalar to convert the units of density
        Should only be used to alter a return values from self.libxrayphysics
        or to do inverse scaling to an argument to a self.libxrayphysics function
        """
        if self.length_units == 'mm':
            return 1.0e-3
        else:
            return 1.0
            
    def areal_density_scalar(self):
        """Returns the scalar to convert the units of areal density
        Should only be used to alter a return values from self.libxrayphysics
        or to do inverse scaling to an argument to a self.libxrayphysics function
        """
        if self.length_units == 'mm':
            return 1.0e-2
        else:
            return 1.0
            
    def thickness_scalar(self):
        """Returns the scalar to convert the units of thickness
        Should only be used to alter a return values from self.libxrayphysics
        or to do inverse scaling to an argument to a self.libxrayphysics function
        """
        if self.length_units == 'mm':
            return 10.0
        else:
            return 1.0

    def atomicMass(self, Z):
        self.libxrayphysics.atomicMass.restype = ctypes.c_float
        self.libxrayphysics.atomicMass.argtypes = [ctypes.c_int]
        return self.libxrayphysics.atomicMass(Z)

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
            else:
                retVal = self.libxrayphysics.sigmaCompound(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigma.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigma.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigma(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigma(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal
                
    def sigma_e(self, Z, gamma):
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.sigmaeCompound.restype = ctypes.c_float
            self.libxrayphysics.sigmaeCompound.argtypes = [ctypes.c_char_p, ctypes.c_float]
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaeCompound(chemicalFormula, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmaeCompound(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmae.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmae.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmae(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmae(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal

    def rhoe(self, Z, massDensity):
        #massDensity * sigma = rhoe*sigma_e
        return massDensity * self.sigma(Z, 1.0) / self.sigma_e(Z, 1.0)
        
    def rho(self, Z, electronDensity):
        #rho * sigma = rhoe*sigma_e
        return electronDensity * self.sigma_e(Z, 1.0) / self.sigma(Z, 1.0)
        
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
            else:
                retVal = self.libxrayphysics.sigmaCompoundPE(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaPE.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaPE.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaPE(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmaPE(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal
            
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
            else:
                retVal = self.libxrayphysics.sigmaCompoundCS(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaCS.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaCS.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaCS(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmaCS(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal
        
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
            else:
                retVal = self.libxrayphysics.sigmaCompoundRS(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaRS.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaRS.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaRS(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmaRS(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal
            
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
            else:
                retVal = self.libxrayphysics.sigmaCompoundPP(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaPP.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaPP.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaPP(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmaPP(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal
            
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
            else:
                retVal = self.libxrayphysics.sigmaCompoundTP(chemicalFormula, float(gamma))
        else:
            self.libxrayphysics.sigmaTP.argtypes = [ctypes.c_float, ctypes.c_float]
            self.libxrayphysics.sigmaTP.restype = ctypes.c_float
            if type(gamma) is np.ndarray:
                retVal = gamma.copy()
                for n in range(gamma.size):
                    retVal[n] = self.libxrayphysics.sigmaTP(Z, float(gamma[n]))
            else:
                retVal = self.libxrayphysics.sigmaTP(Z, gamma)
        retVal *= self.cross_section_scalar()
        return retVal
        
    def incoherentScatterDistribution(self, Z, gamma, theta, doNormalize=False):
        #float incoherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.incoherentScatterDistributionCompound.restype = ctypes.c_float
            self.libxrayphysics.incoherentScatterDistributionCompound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float]
            if type(theta) is np.ndarray:
                retVal = theta.copy()
                for n in range(theta.size):
                    retVal[n] = self.libxrayphysics.incoherentScatterDistributionCompound(chemicalFormula, gamma, theta[n])
            else:
                retVal = self.libxrayphysics.incoherentScatterDistributionCompound(chemicalFormula, gamma, theta)
        else:
            self.libxrayphysics.incoherentScatterDistribution.restype = ctypes.c_float
            self.libxrayphysics.incoherentScatterDistribution.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float]
            if type(theta) is np.ndarray:
                retVal = theta.copy()
                for n in range(theta.size):
                    retVal[n] = self.libxrayphysics.incoherentScatterDistribution(float(Z), gamma, theta[n])
            else:
                retVal = self.libxrayphysics.incoherentScatterDistribution(float(Z), gamma, theta)
        retVal *= self.cross_section_scalar()
        if doNormalize:
            retVal = retVal / self.incoherentScatterDistribution_normalizationFactor(Z, gamma)
        return retVal
    
    def coherentScatterDistribution(self, Z, gamma, theta, doNormalize=False):
        #float coherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.coherentScatterDistributionCompound.restype = ctypes.c_float
            self.libxrayphysics.coherentScatterDistributionCompound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float]
            if type(theta) is np.ndarray:
                retVal = theta.copy()
                for n in range(theta.size):
                    retVal[n] = self.libxrayphysics.coherentScatterDistributionCompound(chemicalFormula, gamma, theta[n])
            else:
                retVal = self.libxrayphysics.coherentScatterDistributionCompound(chemicalFormula, gamma, theta)
        else:
            self.libxrayphysics.coherentScatterDistribution.restype = ctypes.c_float
            self.libxrayphysics.coherentScatterDistribution.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float]
            if type(theta) is np.ndarray:
                retVal = theta.copy()
                for n in range(theta.size):
                    retVal[n] = self.libxrayphysics.coherentScatterDistribution(float(Z), gamma, theta[n])
            else:
                retVal = self.libxrayphysics.coherentScatterDistribution(float(Z), gamma, theta)
        retVal *= self.cross_section_scalar()
        if doNormalize:
            retVal = retVal / self.coherentScatterDistribution_normalizationFactor(Z, gamma)
        return retVal
            
    ###
    def incoherentScatterDistribution_normalizationFactor(self, Z, gamma):
        #float incoherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.incoherentScatterDistributionCompound_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.incoherentScatterDistributionCompound_normalizationFactor.argtypes = [ctypes.c_char_p, ctypes.c_float]
            return self.libxrayphysics.incoherentScatterDistributionCompound_normalizationFactor(chemicalFormula, gamma) * self.cross_section_scalar()
        else:
            self.libxrayphysics.incoherentScatterDistribution_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.incoherentScatterDistribution_normalizationFactor.argtypes = [ctypes.c_float, ctypes.c_float]
            return self.libxrayphysics.incoherentScatterDistribution_normalizationFactor(float(Z), gamma) * self.cross_section_scalar()
    
    def coherentScatterDistribution_normalizationFactor(self, Z, gamma):
        #float coherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.coherentScatterDistributionCompound_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.coherentScatterDistributionCompound_normalizationFactor.argtypes = [ctypes.c_char_p, ctypes.c_float]
            return self.libxrayphysics.coherentScatterDistributionCompound_normalizationFactor(chemicalFormula, gamma) * self.cross_section_scalar()
        else:
            self.libxrayphysics.coherentScatterDistribution_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.coherentScatterDistribution_normalizationFactor.argtypes = [ctypes.c_float, ctypes.c_float]
            return self.libxrayphysics.coherentScatterDistribution_normalizationFactor(float(Z), gamma) * self.cross_section_scalar()
    ###
    
    def simulateSpectra(self, kV, takeOffAngle=11.0, Z=74, gammas=None):
        #simulateSpectra(float kV, float takeOffAngle, int Z, float* gammas, int N, float* output)
        
        if gammas is None:
            maxEnergy = max(1,int(np.ceil(kV)))
            minEnergy = max(1,int(0.1*np.ceil(kV)))
            T_E = max(1, int(np.floor(maxEnergy/100.0)))
            if T_E > 1:
                N = int(np.ceil(float(maxEnergy - minEnergy) / float(T_E)))
                minEnergy = maxEnergy - N*T_E
            gammas = np.ascontiguousarray(np.array(range(minEnergy,maxEnergy+1,T_E)), dtype=np.float32)
        if Z is None:
            Z = 74
        
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
        
    def normalizeSpectrum(self, spectralResponse, gammas):
        self.libxrayphysics.normalizeSpectrum.restype = ctypes.c_bool
        self.libxrayphysics.normalizeSpectrum.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
        return self.libxrayphysics.normalizeSpectrum(spectralResponse, gammas, gammas.size)
        
    def effectiveZ(self, chemicalFormula, min_energy=10.0, max_energy=100.0, arealDensity=0.0):
        if isinstance(chemicalFormula, str):
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            #float effectiveZ(const char* chemForm, float min_energy, float max_energy, float arealDensity);
            self.libxrayphysics.effectiveZ.restype = ctypes.c_float
            self.libxrayphysics.effectiveZ.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ctypes.c_float]
            return self.libxrayphysics.effectiveZ(chemicalFormula, min_energy, max_energy, arealDensity/self.areal_density_scalar())
        else:
            print('Error: first argument must be a chemical formula string')
            return 0.0
        
    def effectiveAttenuation(self, Z, density, thickness, spectralResponse, gammas):
        #float effectiveAttenuation(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.effectiveAttenuation_compound.restype = ctypes.c_float
            self.libxrayphysics.effectiveAttenuation_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveAttenuation_compound(chemicalFormula, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)
        else:
            self.libxrayphysics.effectiveAttenuation.restype = ctypes.c_float
            self.libxrayphysics.effectiveAttenuation.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveAttenuation(Z, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)
    
    def effectiveEnergy(self, Z, density, thickness, spectralResponse, gammas):
        #float effectiveEnergy(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.effectiveEnergy_compound.restype = ctypes.c_float
            self.libxrayphysics.effectiveEnergy_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveEnergy_compound(chemicalFormula, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)
        else:
            self.libxrayphysics.effectiveEnergy.restype = ctypes.c_float
            self.libxrayphysics.effectiveEnergy.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveEnergy(Z, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)
            
    def transmission(self, Z, density, thickness, spectralResponse, gammas):
        #float transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.transmission_compound.restype = ctypes.c_float
            self.libxrayphysics.transmission_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.transmission_compound(chemicalFormula, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)
        else:
            self.libxrayphysics.transmission.restype = ctypes.c_float
            self.libxrayphysics.transmission.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.transmission(Z, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)
            
    def setBHlookupTable(self, spectralResponse, gammas, Z, referenceEnergy=0.0, T_atten=0.0, N_atten=0):
        if N_atten <= 0 or T_atten <= 0.0:
            max_lac = 48.0
            T_atten = 1.0e-3 # about 66 counts max
            #T_atten = 2.0e-4 # about 0.06MB; 14 counts on a 16-bit detector
            #T_atten = 1.0e-4 # about 0.12MB
            #T_atten = 1.0e-5 # about 1.2MB
            #T_atten = 1.0e-6 # about 12MB
            N_atten = int(np.ceil(max_lac / T_atten)) + 1
            max_lac = float(N_atten-1)*T_atten
        
        N_gamma = gammas.size
        if referenceEnergy <= 0.0:
            #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
            referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spectralResponse, gammas)
            print('referenceEnergy = ' + str(referenceEnergy))
        LUT = np.zeros(N_atten, dtype=np.float32)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.setBHlookupTable_compound.restype = ctypes.c_bool
            self.libxrayphysics.setBHlookupTable_compound.argtypes = [ctypes.c_char_p, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int, ctypes.c_float]
            self.libxrayphysics.setBHlookupTable_compound(chemicalFormula, spectralResponse, gammas, N_gamma, LUT, T_atten, N_atten, referenceEnergy)
        else:
            self.libxrayphysics.setBHlookupTable.restype = ctypes.c_bool
            self.libxrayphysics.setBHlookupTable.argtypes = [ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int, ctypes.c_float]
            self.libxrayphysics.setBHlookupTable(Z, spectralResponse, gammas, N_gamma, LUT, T_atten, N_atten, referenceEnergy)
        return LUT, T_atten

    def setBHClookupTable(self, spectralResponse, gammas, Z, referenceEnergy=0.0, T_atten=0.0, N_atten=0):
        #bool setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_atten, int N_atten, float referenceEnergy);
        
        if N_atten <= 0 or T_atten <= 0.0:
            max_atten = 12.0
            T_atten = 1.0e-3 # about 66 counts max
            #T_atten = 2.0e-4 # about 0.06MB; 14 counts on a 16-bit detector
            #T_atten = 1.0e-4 # about 0.12MB
            #T_atten = 1.0e-5 # about 1.2MB
            #T_atten = 1.0e-6 # about 12MB
            N_atten = int(np.ceil(max_atten / T_atten)) + 1
            max_atten = float(N_atten-1)*T_atten
        
        N_gamma = gammas.size
        if referenceEnergy <= 0.0:
            #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
            referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spectralResponse, gammas)
            print('referenceEnergy = ' + str(referenceEnergy))
        LUT = np.zeros(N_atten, dtype=np.float32)
        if isinstance(Z, str):
            chemicalFormula = Z
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
            self.libxrayphysics.setBHClookupTable_compound.restype = ctypes.c_bool
            self.libxrayphysics.setBHClookupTable_compound.argtypes = [ctypes.c_char_p, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int, ctypes.c_float]
            self.libxrayphysics.setBHClookupTable_compound(chemicalFormula, spectralResponse, gammas, N_gamma, LUT, T_atten, N_atten, referenceEnergy)
        else:
            self.libxrayphysics.setBHClookupTable.restype = ctypes.c_bool
            self.libxrayphysics.setBHClookupTable.argtypes = [ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int, ctypes.c_float]
            self.libxrayphysics.setBHClookupTable(Z, spectralResponse, gammas, N_gamma, LUT, T_atten, N_atten, referenceEnergy)
        return LUT, T_atten
        
    def polynomialBHC(self, spectralResponse, gammas, Z, density, referenceEnergy=0.0, maxThickness=10.0, order=2):
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
        
    def setTwoMaterialBHClookupTable(self, spectralResponse, gammas, sigma_1, sigma_2, referenceEnergy=None, T_atten=0.0, N_atten=0):
        N = gammas.size
        if spectralResponse.size != N or sigma_1.size != N or sigma_2.size != N:
            print('Error: Input arrays must all be the same size!')
            return None, None
            
        if T_atten <= 0.0 or N_atten <= 0:
            T_atten = 0.01
            attenMax = np.ceil(-np.log(2.0**(-self.detectorBits)))
            N_atten = int(np.ceil(attenMax/T_atten))+1
        
        if referenceEnergy is None:
            referenceEnergy = np.floor(self.meanEnergy(spectralResponse, gammas))
            print('Using reference energy: ' + str(referenceEnergy) + ' keV')
        
        sigmas = np.zeros((2,gammas.size), dtype=np.float32)
        sigmas[0,:] = sigma_1[:]
        sigmas[1,:] = sigma_2[:]
        
        LUT = np.zeros((N_atten,N_atten), dtype=np.float32)
        self.libxrayphysics.setTwoMaterialBHClookupTable.restype = ctypes.c_bool
        self.libxrayphysics.setTwoMaterialBHClookupTable.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int]
        self.libxrayphysics.setTwoMaterialBHClookupTable(spectralResponse, gammas, gammas.size, referenceEnergy, sigmas, LUT, T_atten, N_atten)
        return LUT, T_atten
        
    def setDEDlookupTable(self, spectralResponse_L, spectralResponse_H, gammas, basisFunction_1, basisFunction_2, referenceEnergies=None, T_atten=0.0, N_atten=0):
        #bool generateDEDlookUpTables(float* spectralResponses, float* gammas, int N_gamma, float* referenceEnergies, float* basisFunctions, float* LUT, float T_lac, int N_lac);
        
        N = gammas.size
        if spectralResponse_L.size != N or spectralResponse_H.size != N or basisFunction_1.size != N or basisFunction_2.size != N:
            print('Error: Input arrays must all be the same size!')
            return None, None
        
        if T_atten <= 0.0 or N_atten <= 0:
            #T_atten = np.sqrt(8.0*1.0e-5) # = 0.008944; max interpolation error of 1e-5 max |f''(x)|
            #T_atten = sqrt(8.0*0.5e-5) # = 0.006325; max interpolation error of 0.5e-5 max |f''(x)|
            #T_atten = 0.005

            T_atten = 0.01
            attenMax = np.ceil(-np.log(2.0**(-self.detectorBits)))
            N_atten = int(np.ceil(attenMax/T_atten))+1
            
        if referenceEnergies is None:
            referenceEnergies = np.array([self.meanEnergy(spectralResponse_L, gammas), self.meanEnergy(spectralResponse_H, gammas)], dtype=np.float32)
            referenceEnergies = np.floor(referenceEnergies)
            print('Using reference energies: ' + str(referenceEnergies) + ' keV')
        else:
            referenceEnergies = np.array(referenceEnergies, dtype=np.float32)
            
        spectralResponses = np.zeros((2,spectralResponse_L.size), dtype=np.float32)
        spectralResponses[0,:] = spectralResponse_L[:]
        spectralResponses[1,:] = spectralResponse_H[:]
        
        basisFunctions = np.zeros((2,spectralResponse_L.size), dtype=np.float32)
        basisFunctions[0,:] = basisFunction_1[:]
        basisFunctions[1,:] = basisFunction_2[:]
            
        LUT = np.zeros((3,N_atten,N_atten), dtype=np.float32)
        self.libxrayphysics.generateDEDlookUpTables.restype = ctypes.c_bool
        self.libxrayphysics.generateDEDlookUpTables.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_float, ctypes.c_int]
        self.libxrayphysics.generateDEDlookUpTables(spectralResponses, gammas, gammas.size, referenceEnergies, basisFunctions, LUT, T_atten, N_atten)
        return LUT, T_atten
        
    def PhotoelectricBasis(self, gammas):
        basisFcn = gammas.copy()
        basisFcn[:] = gammas[:]**-3.0
        basisFcn *= self.cross_section_scalar()
        return basisFcn
        
    def ComptonBasis(self, gammas):
        ELECTRON_REST_MASS_ENERGY = 510.975
        CLASSICAL_ELECTRON_RADIUS = 2.8179403267e-13
        AVOGANDROS_NUMBER = 6.0221414107e23
        
        KNconstant = CLASSICAL_ELECTRON_RADIUS*CLASSICAL_ELECTRON_RADIUS*AVOGANDROS_NUMBER
        two_PI_KNconstant = 2.0*np.pi*KNconstant
    
        alpha = gammas.copy() / ELECTRON_REST_MASS_ENERGY
        one_plus_two_alpha = 1.0+2.0*alpha
        log_term = np.log(one_plus_two_alpha)
        basisFcn = (1.0+one_plus_two_alpha)/(2.0*one_plus_two_alpha*one_plus_two_alpha) + ((alpha*alpha-1.0-one_plus_two_alpha)*log_term + 4.0*alpha)/(2.0*alpha*alpha*alpha) * two_PI_KNconstant
        basisFcn *= self.cross_section_scalar()
        return basisFcn
        