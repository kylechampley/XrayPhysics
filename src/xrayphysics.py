import ctypes
import os
import sys
from sys import platform as _platform
from numpy.ctypeslib import ndpointer
import numpy as np


class xrayPhysics:
    """ Python class for tomographicModels bindings
    Usage Example:
    from xrayphysics import *
    physics = xrayPhysics()
    ...
    """

    def __init__(self, lib_dir=""):
        if len(lib_dir) > 0:
            current_dir = lib_dir
        else:
            current_dir = os.path.abspath(os.path.dirname(__file__))

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
            
        self.detectorBits = 16
        self.length_units = 'cm'

    def about(self):
        """Prints info about this package, including the version number"""
        return self.libxrayphysics.about()

    def use_cm(self):
        """Call this function so that all units that use lengths are cm-based (this is the default setting)"""
        self.length_units = 'cm'
        
    def use_mm(self):
        """Call this function so that all units that use lengths are mm-based
        
        The documentation of all functions in this class is written so that the arguments and return values are cm-based,
        but if this function is called one should assume that all function arguments and return values are mm-based.
        
        """
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
            
    def fixArray(self, x):
        if type(x) is np.ndarray:
            x = np.ascontiguousarray(x, dtype=np.float32)
        return x
        
    def resample(self, gammas, s, gammas_out):
        """Resamples the binning array on which a spectra is defined
        
        Args:
            gammas (numpy array): the energy samples (keV) of the given spectra
            s (numpy array): the given spectra to resample
            gammmas_out (numpy array): the energy samples (keV) to rebin to
        
        Returns:
            resampled spectra
        """
        s_out = gammas_out.copy()
        
        in_hi = 0.5*(gammas + np.roll(gammas,-1))
        in_hi[-1] = gammas[-1]+0.5*(gammas[-1]-gammas[-2])
        
        in_lo = 0.5*(np.roll(gammas,1) + gammas)
        in_lo[0] = max(0.0, gammas[0]-0.5*(gammas[1]-gammas[0]))
        
        out_hi = 0.5*(gammas_out + np.roll(gammas_out,-1))
        out_hi[-1] = gammas_out[-1]+0.5*(gammas_out[-1]-gammas_out[-2])
        
        out_lo = 0.5*(np.roll(gammas_out,1) + gammas_out)
        out_lo[0] = max(0.0, gammas_out[0]-0.5*(gammas_out[1]-gammas_out[0]))
        
        s_out[:] = 0.0
        for i in range(gammas_out.size):
            for j in range(gammas.size):
                intersection = max(0.0, min(out_hi[i], in_hi[j]) - max(out_lo[i], in_lo[j]))
                s_out[i] += intersection * s[j]
        
        return s_out

    def elementSymbolToAtomicNumber(self, elementStr):
        if sys.version_info[0] == 3:
            elementStr = bytes(str(elementStr), 'ascii')
        self.libxrayphysics.elementSymbolToAtomicNumber.restype = ctypes.c_int
        self.libxrayphysics.elementSymbolToAtomicNumber.argtypes = [ctypes.c_char_p]
        return self.libxrayphysics.elementSymbolToAtomicNumber(elementStr)

    def atomicMass(self, Z):
        """Returns the atomic mass of the given atomic number"""
        self.libxrayphysics.atomicMass.restype = ctypes.c_float
        self.libxrayphysics.atomicMass.argtypes = [ctypes.c_int]
        if isinstance(Z, str):
            Z = self.elementSymbolToAtomicNumber(Z)
        return self.libxrayphysics.atomicMass(Z)
        
    def density(self, Z):
        """Returns the mass density of the given atomic number or element"""
        return self.massDensity(Z)
        
    def massDensity(self, Z):
        """Returns the mass density of the given atomic number or element"""
        if isinstance(Z, str):
            atomicNumber = self.elementSymbolToAtomicNumber(Z)
            if atomicNumber < 1:
                if Z in materialDensities.keys():
                    return materialDensities[Z] * self.density_scalar()
                else:
                    return 0.0
            else:
                Z = atomicNumber
        self.libxrayphysics.massDensity.restype = ctypes.c_float
        self.libxrayphysics.massDensity.argtypes = [ctypes.c_int]
        return self.libxrayphysics.massDensity(Z) * self.density_scalar()

    def mu(self, Z, gamma, massDensity=None):
        """Returns the Linear Attenuation Coefficient (LAC) of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the LAC
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The LAC (cm^-1) of the material at the specified energy or energies
        
        """
        if massDensity is None:
            massDensity = self.massDensity(Z)
        return massDensity * self.sigma(Z, gamma)
        
    def muPE(self, Z, gamma, massDensity=None):
        """Returns the Photoelectric component of the Linear Attenuation Coefficient (LAC) of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the LAC
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The Photoelectric component of the LAC (cm^-1) of the material at the specified energy or energies
        
        """
        if massDensity is None:
            massDensity = self.massDensity(Z)
        return massDensity * self.sigmaPE(Z, gamma)
        
    def muCS(self, Z, gamma, massDensity=None):
        """Returns the Compton Scatter component of the Linear Attenuation Coefficient (LAC) of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the LAC
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The Compton Scatter component of the LAC (cm^-1) of the material at the specified energy or energies
        
        """
        if massDensity is None:
            massDensity = self.massDensity(Z)
        return massDensity * self.sigmaCS(Z, gamma)
        
    def muRS(self, Z, gamma, massDensity=None):
        """Returns the Rayleigh Scatter component of the Linear Attenuation Coefficient (LAC) of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the LAC
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The Rayleigh Scatter component of the LAC (cm^-1) of the material at the specified energy or energies
        
        """
        if massDensity is None:
            massDensity = self.massDensity(Z)
        return massDensity * self.sigmaRS(Z, gamma)
        
    def muPP(self, Z, gamma, massDensity=None):
        """Returns the Pair Production component of the Linear Attenuation Coefficient (LAC) of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the LAC
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The Pair Production component of the LAC (cm^-1) of the material at the specified energy or energies
        
        """
        if massDensity is None:
            massDensity = self.massDensity(Z)
        return massDensity * self.sigmaPP(Z, gamma)
        
    def muTP(self, Z, gamma, massDensity=None):
        """Returns the Triplet Production component of the Linear Attenuation Coefficient (LAC) of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the LAC
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The Triplet Production component of the LAC (cm^-1) of the material at the specified energy or energies
        
        """
        if massDensity is None:
            massDensity = self.massDensity(Z)
        return massDensity * self.sigmaTP(Z, gamma)

    def sigma(self, Z, gamma):
        """Returns the mass cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The mass cross section (cm^2/g) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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

    def sigmaPE(self, Z, gamma):
        """Returns the Photoelectric component of the mass cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The Photoelectric component of the mass cross section (cm^2/g) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        """Returns the Compton Scatter component of the mass cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The Compton Scatter component of the mass cross section (cm^2/g) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        """Returns the Rayleigh Scatter component of the mass cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The Rayleigh Scatter component of the mass cross section (cm^2/g) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        """Returns the Pair Production component of the mass cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The Pair Production component of the mass cross section (cm^2/g) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        """Returns the Triplet Production component of the mass cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The Triplet Production component of the mass cross section (cm^2/g) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
    
    def rhoe(self, Z, massDensity):
        """Returns the electron density of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, or mixture of compounds with mass fractions
            massDensity (scalar): the mass density (g/cm^3) of the material
            
        Returns:
            The electron density (mol * e / cm^3) of the material
        
        """
        #massDensity * sigma = rhoe*sigma_e
        return massDensity * self.sigma(Z, 1.0) / self.sigma_e(Z, 1.0)
        
    def rho(self, Z, electronDensity):
        """Returns the mass density of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, or mixture of compounds with mass fractions
            electronDensity (scalar): the electron density (mol * e / cm^3) of the material
            
        Returns:
            The mass density (g / cm^3) of the material
        
        """
        #rho * sigma = rhoe*sigma_e
        return electronDensity * self.sigma_e(Z, 1.0) / self.sigma(Z, 1.0)
                
    def sigma_e(self, Z, gamma):
        """Returns the electron cross section of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar or numpy array): the energy or energies (keV) for which to calculate the mass cross section
            
        Returns:
            The electron cross section (cm^2/(mol * e)) of the material at the specified energy or energies
        
        """
        gamma = self.fixArray(gamma)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
    
    def ComptonScatterDistribution(self, Z, gamma, theta, doNormalize=False):
        """Returns the Compton Scatter distribution of the given material
        
        Alias for incoherentScatterDistribution
        
        """
        return self.incoherentScatterDistribution(Z, gamma, theta, doNormalize)
    
    def incoherentScatterDistribution(self, Z, gamma, theta, doNormalize=False):
        """Returns the Incoherent Scatter distribution of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar): the energy (keV) for which to calculate the Scatter distribution
            theta (scalar or numpy array): the angle(s) (degrees) for which to calculate the Scatter distribution
            doNormalize (bool): if True, normalizes the values so that they can express as a probability density function
            
        Returns:
            The Incoherent Scatter distribution of the material at the specified energy and angles
        
        """
        #float incoherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
    
    def RayleighScatterDistribution(self, Z, gamma, theta, doNormalize=False):
        """Returns the Rayleigh Scatter distribution of the given material
        
        Alias for coherentScatterDistribution
        
        """
        return self.coherentScatterDistribution(Z, gamma, theta, doNormalize)
    
    def coherentScatterDistribution(self, Z, gamma, theta, doNormalize=False):
        """Returns the Coherent Scatter distribution of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar): the energy (keV) for which to calculate the Scatter distribution
            theta (scalar or numpy array): the angle(s) (degrees) for which to calculate the Scatter distribution
            doNormalize (bool): if True, normalizes the values so that they can express as a probability density function
            
        Returns:
            The Coherent Scatter distribution of the material at the specified energy and angles
        
        """
        #float coherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
    def ComptonScatterDistribution_normalizationFactor(self, Z, gamma):
        """Returns the Compton Scatter distribution normalization factor of the given material
        
        Alias for incoherentScatterDistribution_normalizationFactor
        
        """
        return self.incoherentScatterDistribution_normalizationFactor(Z, gamma)
    
    def incoherentScatterDistribution_normalizationFactor(self, Z, gamma):
        """Returns the Incoherent Scatter distribution normalization factor of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar): the energy (keV) for which to calculate the normalization factor
            
        Returns:
            The Incoherent Scatter distribution normalization factor of the material at the specified energy
        
        """
        #float incoherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
                
            self.libxrayphysics.incoherentScatterDistributionCompound_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.incoherentScatterDistributionCompound_normalizationFactor.argtypes = [ctypes.c_char_p, ctypes.c_float]
            return self.libxrayphysics.incoherentScatterDistributionCompound_normalizationFactor(chemicalFormula, gamma) * self.cross_section_scalar()
        else:
            self.libxrayphysics.incoherentScatterDistribution_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.incoherentScatterDistribution_normalizationFactor.argtypes = [ctypes.c_float, ctypes.c_float]
            return self.libxrayphysics.incoherentScatterDistribution_normalizationFactor(float(Z), gamma) * self.cross_section_scalar()
    
    def RayleighScatterDistribution_normalizationFactor(self, Z, gamma):
        """Returns the Rayleigh Scatter distribution normalization factor of the given material
        
        Alias for coherentScatterDistribution_normalizationFactor
        
        """
        return self.coherentScatterDistribution_normalizationFactor(Z, gamma)
    
    def coherentScatterDistribution_normalizationFactor(self, Z, gamma):
        """Returns the Coherent Scatter distribution normalization factor of the given material
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            gamma (scalar): the energy (keV) for which to calculate the normalization factor
            
        Returns:
            The Coherent Scatter distribution normalization factor of the material at the specified energy
        
        """
        #float coherentScatterDistribution(int Z, float gamma, float theta)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
                
            self.libxrayphysics.coherentScatterDistributionCompound_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.coherentScatterDistributionCompound_normalizationFactor.argtypes = [ctypes.c_char_p, ctypes.c_float]
            return self.libxrayphysics.coherentScatterDistributionCompound_normalizationFactor(chemicalFormula, gamma) * self.cross_section_scalar()
        else:
            self.libxrayphysics.coherentScatterDistribution_normalizationFactor.restype = ctypes.c_float
            self.libxrayphysics.coherentScatterDistribution_normalizationFactor.argtypes = [ctypes.c_float, ctypes.c_float]
            return self.libxrayphysics.coherentScatterDistribution_normalizationFactor(float(Z), gamma) * self.cross_section_scalar()
    
    def KleinNishinaScatterDistribution(self, gamma, theta, doNormalize=False):
        ELECTRON_REST_MASS_ENERGY = 510.975
        CLASSICAL_ELECTRON_RADIUS = 2.8179403267e-13
        AVOGANDROS_NUMBER = 6.0221414107e23
        KNconstant = CLASSICAL_ELECTRON_RADIUS*CLASSICAL_ELECTRON_RADIUS*AVOGANDROS_NUMBER
        cos_theta = np.cos(theta*np.pi/180.0)
        P = 1.0 / (1.0 + (gamma / ELECTRON_REST_MASS_ENERGY) * (1.0 - cos_theta))
        retVal = 0.5 * P * P * (P + 1.0 / P - 1.0 + cos_theta * cos_theta) / (self.ComptonBasis(gamma) / KNconstant)
        if doNormalize:
            return retVal / self.KleinNishinaScatterDistribution_normalizationFactor(gamma)
        else:
            return retVal
        
    def KleinNishinaScatterDistribution_normalizationFactor(self, gamma):
        thetas = np.linspace(0.0, 180.0, 1800)
        KN = self.KleinNishinaScatterDistribution(gamma, thetas, False)
        return np.sum(KN*np.sin(thetas*np.pi/180.0)*2.0*np.pi*np.pi/180.0*0.1)
        
    ###
    
    def simulateSpectra(self, kV, takeOffAngle=11.0, Z=74, gammas=None):
        r"""x-ray source spectra model
        
        This model is based on the following paper:
        
        Finkelshtein, A. L., and T. O. Pavlova.
        "Calculation of x-ray tube spectral distributions."
        X-Ray Spectrometry: An International Journal 28, no. 1 (1999): 27-32.
        
        The strength of the characteristic lines were determined using empirical methods.
        
        Args:
            kV (scalar): voltage of a bremsstrahlung spectrum
            takeOffAngle (scalar or list): the take-off angle (degrees) of the reflection anode
            Z (int): the atomic number of the anode material, must be 29, 42, 74, or 79 (Cu, Mo, W, Au)
            gammas (numpy array): the energies for which to model the spectra, if unspecified uses default values
            
        Returns:
            x-ray energy samples, x-ray source spectra model (will be 2D if multiple take-off angles are given)
        
        """
        #simulateSpectra(float kV, float takeOffAngle, int Z, float* gammas, int N, float* output)
        
        
        if gammas is None:
            maxEnergy = max(1,int(np.ceil(kV)))
            minEnergy = max(1,int(0.1*np.ceil(kV)))
            T_E = max(1, int(np.floor(maxEnergy/100.0)))
            if T_E > 1:
                N = int(np.ceil(float(maxEnergy - minEnergy) / float(T_E)))
                minEnergy = maxEnergy - N*T_E
            gammas = np.ascontiguousarray(np.array(range(minEnergy,maxEnergy+1,T_E)), dtype=np.float32)
        else:
            gammas = self.fixArray(gammas)
        if Z is None:
            Z = 74
        
        if isinstance(takeOffAngle, list) or isinstance(takeOffAngle, np.ndarray):
            # user specified several take-off angles, so create a spectra for each
            numSpectrum = len(takeOffAngle)
            s_all = np.zeros((numSpectrum, gammas.size), dtype=np.float32)
            for i in range(numSpectrum):
            
                gammas, sourceSpectrum = self.simulateSpectra(kV, takeOffAngle[i], Z, gammas)
                s_all[i,:] = sourceSpectrum[:]
            return gammas, s_all
        
        sourceSpectrum = np.ascontiguousarray(np.zeros(gammas.size), dtype=np.float32)
        
        self.libxrayphysics.simulateSpectra.restype = ctypes.c_bool
        self.libxrayphysics.simulateSpectra.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]
        self.libxrayphysics.simulateSpectra(kV, takeOffAngle, Z, gammas, gammas.size, sourceSpectrum)
        return gammas, sourceSpectrum
        
    def changeTakeOffAngle(self, kV, takeOffAngle_cur, takeOffAngle_new, Z, gammas, spectrum_cur):
        r"""Change the anode take-off angle of a given x-ray source spectra model
        
        This function uses the Philibert absorption correction factor to modify the given spectra for a different take-off angle.
        See the paper referenced in the simulateSpectra function for a detailed description of this factor.
        
        Args:
            kv (scalar): voltage of a bremsstrahlung spectrum
            takeOffAngle_cur (scalar): the take-off angle (degrees) of the input spectra
            takeOffAngle_new (scalar): the take-off angle (degress) of the output spectra
            Z (int): the atomic number of the anode material, must be 29, 42, 74, or 79 (Cu, Mo, W, Au)
            gammas (numpy array): the energy samples at which the input spectra is defined
            spectrum_cur (numpy array): the current x-ray spectra model (from any source; does not have to be created by this package)
        
        Returns:
            The modified x-ray spectra model
        
        """
        #bool changeTakeOffAngle(float kV, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s)
        gammas = self.fixArray(gammas)
        spectrum_cur = self.fixArray(spectrum_cur)
        self.libxrayphysics.changeTakeOffAngle.restype = ctypes.c_bool
        self.libxrayphysics.changeTakeOffAngle.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]
        spectrum_new = spectrum_cur.copy()
        self.libxrayphysics.changeTakeOffAngle(kV, takeOffAngle_cur, takeOffAngle_new, Z, gammas, gammas.size, spectrum_new)
        return spectrum_new
        
    def detectorResponse(self, chemicalFormula, density, thickness, gammas):
        r"""Analytic detector response model
        
        The model is given by
        
        .. math::
           \begin{eqnarray*}
           \gamma \left[ 1 - e^{-\mu(\gamma)L}\right]
           \end{eqnarray*}
           
        where :math:`\mu` is the LAC of the scintillator, :math:`\gamma` is the x-ray energy (keV), and :math:`L` is the scintillator thickness.
        This provides a descent approximation of the detector response.  A detailed model of the photodiode response of the detector will
        provide a better estimate.
        
        Args:
            chemicalFormula (string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            density (scalar): the mass density (g/cm^3) of the scintillator
            thickness (scalar): the thickness (cm) of the scintillator
            gammas (numpy array): the energy samples (keV) at which to calculate the detector response
        
        Returns:
            A model of the detector response
        """
        if density is None:
            density = self.massDensity(chemicalFormula)
        gammas = self.fixArray(gammas)
        detResp = np.ascontiguousarray(np.zeros(gammas.size), dtype=np.float32)
        for n in range(gammas.size):
            detResp[n] = gammas[n]*(1.0-np.exp(-self.mu(chemicalFormula, gammas[n], density)*thickness))
        return detResp
        
    def filterResponse(self, chemicalFormula, density, thickness, gammas):
        r"""X-ray Filter response model
        
        The model is given by
        
        .. math::
           \begin{eqnarray*}
           e^{-\mu(\gamma)L}
           \end{eqnarray*}
           
        where :math:`\mu` is the LAC of the filter material, :math:`\gamma` is the x-ray energy (keV), and :math:`L` is the material thickness.
        
        Args:
            chemicalFormula (string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            density (scalar): the mass density (g/cm^3) of the filter material
            thickness (scalar): the thickness (cm) of the filter material
            gammas (numpy array): the energy samples (keV) at which to calculate the filter response
        
        Returns:
            A model of the response to the given x-ray filter
        
        """
        if density is None:
            density = self.massDensity(chemicalFormula)
        gammas = self.fixArray(gammas)
        filtResp = np.ascontiguousarray(np.zeros(gammas.size), dtype=np.float32)
        for n in range(gammas.size):
            filtResp[n] = np.exp(-self.mu(chemicalFormula, gammas[n], density)*thickness)
        return filtResp
        
    def meanEnergy(self, spectralResponse, gammas):
        """Calculates the mean energy of a given spectra
        
        Args:
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
            
        Returns:
            The mean energy (keV) of the given spectra
        
        """
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        self.libxrayphysics.meanEnergy.restype = ctypes.c_float
        self.libxrayphysics.meanEnergy.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
        return self.libxrayphysics.meanEnergy(spectralResponse, gammas, gammas.size)
        
    def normalizeSpectrum(self, spectralResponse, gammas):
        """Normalizes the given spectra
        
        Args:
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
            
        Returns:
            The normalized spectra
        
        """
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        self.libxrayphysics.normalizeSpectrum.restype = ctypes.c_bool
        self.libxrayphysics.normalizeSpectrum.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
        return self.libxrayphysics.normalizeSpectrum(spectralResponse, gammas, gammas.size)

    def load_spectra(self, fileName):
        """Load spectra from file (txt or npy)

        This function loads two columns of data.  In the context of this software, the first column are the x-ray energies
        and the second column is either the spectra, detector response, LAC values, etc.  This data can either be (N,2) data
        stored in an npy file or an acsii file of two columns of numbers which is a common format used to store spectra in a file.
        
        Args:
            fileName (str): name of file where spectra information is stored
            
        Returns:
            x-ray energy samples, x-ray source spectra model
        
        """
        if fileName.endswith('.npy'):
            data = np.load(fileName)
            if len(data.shape) != 2:
                print('Error: invalid file format!')
                return None, None
            Es = np.array(data[:,0],dtype=np.float32)
            s = np.array(data[:,1],dtype=np.float32)
            return Es, s
        else:
            data = []
            try:
                with open(fileName, 'r') as file:
                    for line in file:
                        if line[0] != '#':
                            two_numbers = [float(num) for num in line.split()]
                            if len(two_numbers) != 2:
                                print('Error: invalid spectra file format, must be two columns of numbers!')
                                return None, None
                            data.append(two_numbers)
            except FileNotFoundError:
                print('Error: failed to open the file ', fileName)
                return None, None
            except Exception as e:
                print('Error occured while loading the data', str(e))
                return None, None
            
            data = np.array(data)
            Es = np.array(data[:,0],dtype=np.float32)
            s = np.array(data[:,1],dtype=np.float32)
            return Es, s
        
    def save_spectra(self, fileName, spectralResponse, gammas):
        """Save spectra to file (txt or npy)
        
        This function saves two columns of data.  In the context of this software, the first column are the x-ray energies
        and the second column is either the spectra, detector response, LAC values, etc.  This data can either be (N,2) data
        stored in an npy file or an acsii file of two columns of numbers which is a common format used to store spectra in a file.
        
        Args:
            fileName (str): name of file where spectra information is to be saved
            spectralResponse (float 32 numpy array): the x-ray spectra, detector response, LAC values, etc.
            gammas (float 32 numpy array): the x-ray energy samples (keV)
        """
        if fileName.endswith('.npy'):
            data = np.zeros((spectralResponse.size,2),dtype=np.float32)
            data[:,0] = gammas
            data[:,1] = spectralResponse
            np.save(fileName, data)
        else:
            try:
                with open(fileName, 'w') as file:
                    for n in range(spectralResponse.size):
                        line = str(gammas[n]) + ' ' + str("{:e}".format(spectralResponse[n])) + '\n'
                        file.write(line)
            except Exception as e:
                print('Error occured while loading the data', str(e))
        
    def effectiveZ(self, chemicalFormula, min_energy=10.0, max_energy=100.0, arealDensity=0.0):
        r"""Calculate the effective atomic number of a material
        
        Args:
            chemicalFormula (string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            min_energy (scalar): the minimum of the energy in which to perform the calculation
            max_energy (scalar): the maximum of the energy in which to perform the calculation
            arealDensity (scalar): the areal density (g/cm^2) to use for the calculation
        
        Returns:
            The effective atomic number
        
        """
        if isinstance(chemicalFormula, str):
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        return self.effectiveLAC(Z, density, thickness, spectralResponse, gammas)
    
    def effectiveLAC(self, Z, density, thickness, spectralResponse, gammas):
        r"""Calculates the effective Linear Attenuation Coefficient (LAC) of a polychromatic x-ray beam through of a material of specified thickness
        
        The model is given by
        
        .. math::
           \begin{eqnarray*}
           \mu_{eff}(L) &:=& \frac{-\log\left(\int \widehat{s}(\gamma) e^{-\mu(\gamma)L} \, d\gamma\right)}{L} \\
           \widehat{s}(\gamma) &:=& \frac{s(\gamma)}{\int s(\gamma') \, d\gamma'}
           \end{eqnarray*}
           
        where :math:`\mu` is the LAC of the material, :math:`\gamma` is the x-ray energy (keV), :math:`L` is the material thickness, and :math:`s` is the x-ray spectra.
        The limit of :math:`L \rightarrow 0` is given by
        
        .. math::
           \begin{eqnarray*}
           \mu_{eff}(0) &:=& \int \widehat{s}(\gamma) \mu(\gamma) \, d\gamma
           \end{eqnarray*}
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            density (scalar): mass density (g/cm^3) of the material
            thickness (scalar): the thickness (cm) of the material
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
        
        Returns:
            The effective Linear Attenuation Coefficient (LAC, (length units)^-1)
        
        """
        #float effectiveAttenuation(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        
        if density is None:
            density = self.massDensity(Z)
        
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
            if sys.version_info[0] == 3:
                chemicalFormula = bytes(str(chemicalFormula), 'ascii')
                
            self.libxrayphysics.effectiveAttenuation_compound.restype = ctypes.c_float
            self.libxrayphysics.effectiveAttenuation_compound.argtypes = [ctypes.c_char_p, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveAttenuation_compound(chemicalFormula, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)*self.LAC_scalar()
        else:
            self.libxrayphysics.effectiveAttenuation.restype = ctypes.c_float
            self.libxrayphysics.effectiveAttenuation.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int]
            return self.libxrayphysics.effectiveAttenuation(Z, density/self.density_scalar(), thickness/self.thickness_scalar(), spectralResponse, gammas, gammas.size)*self.LAC_scalar()
    
    def effectiveEnergy(self, Z, density, thickness, spectralResponse, gammas):
        r"""Calculate the effective energy of a polychromatic beam passing through a material of a given thickness
        
        The model is given by 
        
        .. math::
           \begin{eqnarray*}
           \mu(\gamma_{eff}) &:=& \mu_{eff}(L)
           \end{eqnarray*}
           
        i.e., the effective energy is the energy at which the LAC of the material at this energy is equal to the effective LAC of the material.
        See the description of the effectiveLAC function
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            density (scalar): mass density (g/cm^3) of the material
            thickness (scalar): the thickness (cm) of the material
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
        
        Returns:
            The monochromatic energy (keV) that would provide the same attenuation as the given polychromatic spectra
        
        """
        #float effectiveEnergy(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        r"""Calculate the transmission through a given material
        
        This model is given by
        
        .. math::
           \begin{eqnarray*}
           \frac{\int s(\gamma) e^{-\mu(\gamma)L}  \, d\gamma}{\int s(\gamma') \, d\gamma'}
           \end{eqnarray*}
        
        where :math:`\mu` is the LAC of the material, :math:`\gamma` is the x-ray energy (keV), :math:`L` is the material thickness, and :math:`s` is the x-ray spectra.
        
        Args:
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            density (scalar): mass density (g/cm^3) of the material
            thickness (scalar): the thickness (cm) of the material
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
        
        Returns:
            The transmission (unitless) through the given filter
        
        """
        #float transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        r"""Calculate a beam hardening transfer function of a given material
        
        This function can be used to model beam hardening.  The model is given by
        
        .. math::
           \begin{eqnarray*}
           a_{poly} &=& -\log\left( \int \widehat{s}(\gamma) e^{-a_{mono}\widehat{\sigma}(\gamma)} \, d\gamma \right) \\
           \widehat{\sigma}(\gamma) &:=& \frac{\sigma(\gamma)}{\sigma(\gamma_{ref})} \\
           \widehat{s}(\gamma) &:=& \frac{s(\gamma)}{\int s(\gamma') \, d\gamma'}
           \end{eqnarray*}
        
        where :math:`a_{poly}` is the polychromatic attenuation, :math:`a_{mono}` is the monochromatic attenuation at the reference energy, :math:`\gamma_{ref}`,
        :math:`\gamma` is the x-ray energy (keV), :math:`\sigma` is the mass cross section of the material, and :math:`s` is the x-ray spectra.
        Note that this function returns :math:`a_{poly}` from a uniform sampling of :math:`a_{mono}`.
        
        Args:
            spectralResponse (numpy array): spectra model (if 2D, assumes multiple spectra are given and will generate a table to each spectra)
            gammas (numpy array): energies at which the spectra model is defined
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            referenceEnergy (scalar): the energy (keV) of the monochromatic attenuation (if not specified uses the mean energy of the spectra)
            T_atten (scalar): the sampling rate of the monochromatic attenuation
            N_atten (int): the number of samples of the monochromatic attenuation
        
        Returns:
            A 1D numpy array that maps monochromatic attenuation to polychromatic attenuation,
            the sampling rate of the monochromatic attenuation
        
        """
        
        if len(spectralResponse.shape) == 2:
            # user specified several spectra, so create a lookup table to each spectra
            if spectralResponse.shape[0] == gammas.size:
                spectralResponse = np.transpose(spectralResponse)
            numSpectrum = spectralResponse.shape[0]
            spec = gammas.copy()
            for i in range(numSpectrum):
                spec[:] = spectralResponse[i,:]
                if referenceEnergy <= 0.0:
                    #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
                    referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spec, gammas)
                    print('referenceEnergy = ' + str(referenceEnergy))
                LUT, T_atten = self.setBHlookupTable(spec, gammas, Z, referenceEnergy, T_atten, N_atten)
                if i == 0:
                    LUTs = np.zeros((numSpectrum, LUT.size), dtype=np.float32)
                LUTs[i,:] = LUT[:]
            return LUTs, T_atten
        
        if N_atten <= 0 or T_atten <= 0.0:
            max_lac = 48.0
            T_atten = 1.0e-3 # about 66 counts max
            #T_atten = 2.0e-4 # about 0.06MB; 14 counts on a 16-bit detector
            #T_atten = 1.0e-4 # about 0.12MB
            #T_atten = 1.0e-5 # about 1.2MB
            #T_atten = 1.0e-6 # about 12MB
            N_atten = int(np.ceil(max_lac / T_atten)) + 1
            max_lac = float(N_atten-1)*T_atten
        
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        N_gamma = gammas.size
        if referenceEnergy <= 0.0:
            #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
            referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spectralResponse, gammas)
            print('referenceEnergy = ' + str(referenceEnergy))
        LUT = np.zeros(N_atten, dtype=np.float32)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        r"""Calculate a beam hardening correction transfer function of a given material
        
        This function can be used to correct for single-material beam hardening.  It calculates the inverse of the setBHlookupTable function, i.e.,
        it returns :math:`a_{mono}` from a uniform sampling of :math:`a_{poly}`.
        
        Args:
            spectralResponse (numpy array): spectra model (if 2D, assumes multiple spectra are given and will generate a table to each spectra)
            gammas (numpy array): energies at which the spectra model is defined
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            referenceEnergy (scalar): the energy (keV) of the monochromatic attenuation (if not specified uses the mean energy of the spectra)
            T_atten (scalar): the sampling rate of the polychromatic attenuation
            N_atten (int): the number of samples of the polychromatic attenuation
        
        Returns:
            A 1D numpy array that maps polychromatic attenuation to monochromatic attenuation,
            the sampling rate of the polychromatic attenuation
        
        """
        #bool setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_atten, int N_atten, float referenceEnergy);
        
        if len(spectralResponse.shape) == 2:
            # user specified several spectra, so create a lookup table to each spectra
            if spectralResponse.shape[0] == gammas.size:
                spectralResponse = np.transpose(spectralResponse)
            numSpectrum = spectralResponse.shape[0]
            spec = gammas.copy()
            for i in range(numSpectrum):
                spec[:] = spectralResponse[i,:]
                if referenceEnergy <= 0.0:
                    #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
                    referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spec, gammas)
                    print('referenceEnergy = ' + str(referenceEnergy))
                LUT, T_atten = self.setBHClookupTable(spec, gammas, Z, referenceEnergy, T_atten, N_atten)
                if i == 0:
                    LUTs = np.zeros((numSpectrum, LUT.size), dtype=np.float32)
                LUTs[i,:] = LUT[:]
            return LUTs, T_atten
        
        if N_atten <= 0 or T_atten <= 0.0:
            max_atten = 12.0
            T_atten = 1.0e-3 # about 66 counts max
            #T_atten = 2.0e-4 # about 0.06MB; 14 counts on a 16-bit detector
            #T_atten = 1.0e-4 # about 0.12MB
            #T_atten = 1.0e-5 # about 1.2MB
            #T_atten = 1.0e-6 # about 12MB
            N_atten = int(np.ceil(max_atten / T_atten)) + 1
            max_atten = float(N_atten-1)*T_atten
        
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        N_gamma = gammas.size
        if referenceEnergy <= 0.0:
            #referenceEnergy = self.meanEnergy(spectralResponse, gammas)
            referenceEnergy = self.effectiveEnergy(Z, 1.0, 0.0, spectralResponse, gammas)
            print('referenceEnergy = ' + str(referenceEnergy))
        LUT = np.zeros(N_atten, dtype=np.float32)
        if isinstance(Z, str):
            chemicalFormula = Z
            if chemicalFormula in materialFormulas.keys():
                chemicalFormula = materialFormulas[chemicalFormula]
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
        r"""Calculate the coefficients of a polynomial-based beam hardening correction transfer function of a given material
        
        This function can be used to approximately correct for single-material beam hardening;
        it serves the same purpose as the setBHClookupTable except the correction is approximated by a polynomial which is convenient for some purposes.
        Applying this correction is done as follows:
        
        .. math::
           \begin{eqnarray*}
           a_{mono} &=& \sum_n c_n a_{poly}^n
           \end{eqnarray*}
           
        where :math:`c_n` are the polynomial coefficients calculated by this function, :math:`a_{poly}` is the polychromatic attenuation, :math:`a_{mono}` is the monochromatic attenuation at the reference energy.
        
        Args:
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
            Z (scalar or string): atomic number, chemical formula, mixture of compounds with mass fractions, or member of the material library
            density (scalar): the mass density (g/cm^3) of the matial
            referenceEnergy (scalar): the energy (keV) of the monochromatic attenuation (if not specified uses the mean energy of the spectra)
            maxThickness (scalar): the maximum thickness (cm) that one expects to penetrate of this material (defines the range of values for the polynomial fit)
            order (int): the order of the polynomial
        
        Returns:
            A numpy array of the coefficients of the BHC polynomial
        
        """
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
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
        r"""Calculate a two-material beam hardening correction transfer function of the two given materials
        
        Args:
            spectralResponse (numpy array): spectra model
            gammas (numpy array): energies at which the spectra model is defined
            sigma_1 (numpy array): the mass cross section of the lower attenuating material
            sigma_2 (numpy array): the mass cross section of the higher attenuating material
            referenceEnergy (scalar): the energy (keV) of the monochromatic attenuation (if not specified uses the mean energy of the spectra)
            T_atten (scalar): the sampling rate of the polychromatic attenuation (if unspecified, uses default values)
            N_atten (int): the number of samples of the polychromatic attenuation (if unspecified, uses default values)
        
        Returns:
            A 2D numpy array that maps polychromatic attenuation to monochromatic attenuation
        
        """
        gammas = self.fixArray(gammas)
        spectralResponse = self.fixArray(spectralResponse)
        sigma_1 = self.fixArray(sigma_1)
        sigma_2 = self.fixArray(sigma_2)
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
        r"""Calculate a dual energy decomposition transfer function
        
        This function calculates a look up table to perform dual energy decomposition.  The decomposed data are virtual monochromatic attenuation values
        where the monochromatic energies are given by the referenceEnergies input argument.  We describe how to get the basis coefficients instead of
        monochromatic attenuation values below.
        
        This function solves the following
        
        .. math::
           \begin{eqnarray*}
           \widehat{a}_L, \widehat{a}_H &:=& \text{argmin} \left[ p_L + \log\left( \int \widehat{s}_L(\gamma) e^{-b_L(\gamma)a_L - b_H(\gamma)a_H} \, d\gamma \right) \right]^2 \\
           &+& \left[ p_H + \log\left( \int \widehat{s}_H(\gamma) e^{-b_L(\gamma)a_L - b_H(\gamma)a_H} \, d\gamma \right) \right]^2 
           \end{eqnarray*}
           
        where :math:`p_L` and :math:`p_H` are the measured low and high energy attenuation data, :math:`\widehat{a}_L` and :math:`\widehat{a}_H`
        are the low and high virtual monochromatic attenuation data, :math:`s_L` and :math:`s_H` are the low and high energy spectra, and
        :math:`\gamma_L` and :math:`\gamma_H` are the low and high virtual energies (keV).
        The basis functions are given by
        
        .. math::
           \begin{eqnarray*}
           \begin{bmatrix} b_L(\gamma) \\ b_H(\gamma) \end{bmatrix} &:=& \begin{bmatrix} b_1(\gamma_L) & b_1(\gamma_H) \\ b_2(\gamma_L) & b_2(\gamma_H) \end{bmatrix}^{-1} \begin{bmatrix} b_1(\gamma) \\ b_2(\gamma) \end{bmatrix}
           \end{eqnarray*}
           
        where :math:`b_1` and :math:`b_2` are the given basis functions (e.g., compton/ photoelectric, LAC of two materials, PCA basis, etc.)
        
        Note that the normalized spectra are defined by
        
        .. math::
           \begin{eqnarray*}
           \widehat{s}_L(\gamma) &:=& \frac{s_L(\gamma)}{\int s_L(\gamma') \, d\gamma'} \\
           \widehat{s}_H(\gamma) &:=& \frac{s_H(\gamma)}{\int s_H(\gamma') \, d\gamma'}
           \end{eqnarray*}
           
        If one wishes to get the basis coefficients then just apply the following transformation
        
        .. math::
           \begin{eqnarray*}
           \begin{bmatrix} a_1 \\ a_2 \end{bmatrix} &:=& \begin{bmatrix} b_1(\gamma_L) & b_2(\gamma_L) \\ b_1(\gamma_H) & b_2(\gamma_H) \end{bmatrix}^{-1} \begin{bmatrix} a_L \\ a_H \end{bmatrix}
           \end{eqnarray*}
        
           
        A full description of this algorithm is provided in the following reference:
        
        Champley, Kyle M., Stephen G. Azevedo, Isaac M. Seetho, Steven M. Glenn, Larry D. McMichael, Jerel A. Smith, Jeffrey S. Kallman, William D. Brown, and Harry E. Martz.
        "Method to extract system-independent material properties from dual-energy X-ray CT."
        IEEE Transactions on Nuclear Science 66, no. 3 (2019): 674-686.
        
        Args:
            spectralResponse_L (numpy array): low energy spectra model
            spectralResponse_H (numpy array): high energy spectra model
            gammas (numpy array): energies at which both the spectra models are defined
            basisFunction_1 (numpy array): first basis function
            basisFunction_2 (numpy array): second basis function
            referenceEnergies (2-element numpy array): the low and high energy (keV) energies of the synthesized monochromatic attenuation
            T_atten (scalar): the sampling rate of the polychromatic attenuation (if unspecified, uses default values)
            N_atten (int): the number of samples of the polychromatic attenuation (if unspecified, uses default values)
            
        Returns:
            A 2D numpy array that maps the low and high polychromatic attenuation to the low and high (at referenceEnergies) monochromatic attenuation
        
        """
        #bool generateDEDlookUpTables(float* spectralResponses, float* gammas, int N_gamma, float* referenceEnergies, float* basisFunctions, float* LUT, float T_lac, int N_lac);
        
        gammas = self.fixArray(gammas)
        spectralResponse_L = self.fixArray(spectralResponse_L)
        spectralResponse_H = self.fixArray(spectralResponse_H)
        basisFunction_1 = self.fixArray(basisFunction_1)
        basisFunction_2 = self.fixArray(basisFunction_2)
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
        """Returns the photoelectric basis function"""
        gammas = self.fixArray(gammas)
        basisFcn = gammas.copy()
        basisFcn[:] = gammas[:]**-3.0
        basisFcn *= self.cross_section_scalar()
        return basisFcn
        
    def ComptonBasis(self, gammas):
        """Returns the Compton basis function"""
        gammas = self.fixArray(gammas)
        ELECTRON_REST_MASS_ENERGY = 510.975
        CLASSICAL_ELECTRON_RADIUS = 2.8179403267e-13
        AVOGANDROS_NUMBER = 6.0221414107e23
        
        KNconstant = CLASSICAL_ELECTRON_RADIUS*CLASSICAL_ELECTRON_RADIUS*AVOGANDROS_NUMBER
        two_PI_KNconstant = 2.0*np.pi*KNconstant
    
        if type(gammas) is np.ndarray:
            alpha = gammas.copy() / ELECTRON_REST_MASS_ENERGY
        else:
            alpha = gammas / ELECTRON_REST_MASS_ENERGY
        one_plus_two_alpha = 1.0+2.0*alpha
        log_term = np.log(one_plus_two_alpha)
        basisFcn = (1.0+one_plus_two_alpha)/(2.0*one_plus_two_alpha*one_plus_two_alpha) + ((alpha*alpha-1.0-one_plus_two_alpha)*log_term + 4.0*alpha)/(2.0*alpha*alpha*alpha) * two_PI_KNconstant
        basisFcn *= self.cross_section_scalar()
        return basisFcn
    
    def PCAbases(self, listOfMaterials, gammas):
        """Returns the two basis functions derived from the PCA of several material cross sections"""
        gammas = self.fixArray(gammas)
        N = len(listOfMaterials)
        if N < 2:
            print('Error: must give at least two materials')
            return None, None
        X = np.zeros((N,gammas.size))
        for n in range(N):
            X[n,:] = self.sigma(listOfMaterials[n], gammas)
        U, S, V = np.linalg.svd(X)
        b_1 = V[0,:]
        b_2 = V[1,:]
        
        # Normalize the basis functions (not really necessary)
        b_1 *= b_1[0]
        b_2 *= b_2[0]
        b_1 = b_1/np.sqrt(np.sum(b_1**2))
        b_2 = b_2/np.sqrt(np.sum(b_2**2))
        
        return b_1, b_2

    def print_material_library(self):
        """Prints all materials in the material library"""
        for key in materialFormulas:
            #print(key + ': ' + str(materialFormulas[key]) + ', ' + str(materialDensities[key]) + ' g/cm^3')
            if len(key) < 4:
                print(key + ':\t\t' + str(materialDensities[key]) + ' g/cm^3' + ',\t' + str(materialFormulas[key]))
            elif len(key) < 7:
                print(key + ':\t\t' + str(materialDensities[key]) + ' g/cm^3' + ',\t' + str(materialFormulas[key]))
            elif len(key) < 11:
                print(key + ':\t' + str(materialDensities[key]) + ' g/cm^3' + ',\t' + str(materialFormulas[key]))
            else:
                print(key + ': ' + str(materialDensities[key]) + ' g/cm^3' + ',\t' + str(materialFormulas[key]))
            
            
    def get_material_library(self):
        return materialFormulas, materialDensities

# Define the material library
materialFormulas = dict([('vacuum', '78.084*N2 + 20.946*O2 + 0.9340*Ar + 0.0397*CO2'),
('air', '78.084*N2 + 20.946*O2 + 0.9340*Ar + 0.0397*CO2'),
('foam', 'C8H8'),
('water', 'H2O'),
('bone_cortical', '0.047234*H + 0.144330*C + 0.041990*N + 0.446096*O + 0.002200*Mg + 0.104970*P + 0.003150*S + 0.209930*Ca + 0.000100*Zn'),
('bone_compact','0.063984*H+0.278000*C+0.027000*N+0.410016*O+0.002000*Mg+0.070000*P+0.002000*S+0.147000*Ca'),
('adipose_tissue', '0.119477*H + 0.637240*C + 0.007970*N + 0.232333*O + 0.000500*Na + 0.000020*Mg + 0.000160*P + 0.000730*S + 0.001190*Cl + 0.000320*K + 0.000020*Ca + 0.000020*Fe + 0.000020*Zn'),
('lung', '0.10132562304283013*H + 0.10235810831090612*C + 0.028663471831760927*N + 0.75742799115584314*O + 0.0018408652066471242*Na + 0.00073034326133282642*Mg + 0.00080037617680309744*P + 0.0022510579972587113*S + 0.0026612507878702989*Cl + 0.0019409122287475113*K'),
('ethanol', 'C2H6O'),
('HP', 'H2O2'),
('saltWater', '96.5*H2O + 3.5*NaCl'),
('graphite', 'C'),
('diamond', 'C'),
('paper', 'C6H10O5'),
('rubber', 'C5H8'),
('wood', 'C6H12O6'),
('delrin', 'CH2O'),
('teflon', 'C2F4'),
('PMMA', 'C5O2H8'),
('acrylic', 'C5O2H8'),
('Polyethylene', 'C2H4'),
('PE', 'C2H4'),
('LDPE', 'C2H4'),
('HDPE', 'C2H4'),
('Polypropylene', 'C3H6'),
('PP', 'C3H6'),
('PET', 'C10H8O4'),
('Polycarbonate', 'C15H16O2'),
('PVC', 'C2H3Cl'),
('PVDF', 'C2H2F2'),
('melamine', 'C3H6N6'),
('PETG', 'C14H20O5S'),
('PPSU', 'C24H16O4S'),
('PEEK', 'C19H12O3'),
('silicone', '17.57539*C+3.666654*H+37.47571*O+41.28225*Si'),
('glass', '72*SiO2 + 14.2*Na2O + 10*CaO + 2.5*MgO + 0.6*Al2O3'),
('eglass', '54*SiO2 + 14*Al2O3 + 22*CaO + 10*B2O3'),
('aluminumOxide', 'Al2O3'),
('titaniumDioxide', 'TiO2'),
('brass', 'Cu3Zn2'),
('Inconel', '52.5*Ni + 19*Cr + 18*Fe + 5*Nb + 3*Mo + 1*Ti + 0.5*Al + 1*Co'),
('stainlessSteel', '68.5*Fe+17*Cr+12*Ni+2.5*Mo'),
('Ti64', '6*Al+4*V+90*Ti'),
('Ti5553', '5*Al+5*V+5*Mo+3*Cr+82*Ti'),
('304L', '70*Fe+20*Cr+10*Ni'),
('316L', '68.5*Fe+17*Cr+12*Ni+2.5*Mo'),
('Al6061', '97.35*Al+1.0*Mg+0.6*Si+0.35*Fe+0.3*Cu+0.1*Cr+0.1*Zn+0.1*Ti+0.1*Mn'),
('aermet100', '0.23*C+13.4*Co+11.1*Ni+3.1*Cr+1.2*Mo+70.97*Fe'),
('NM', 'CH3NO2'),
('PTM', 'C5H8N4O12'),
('RDX', '94*C3H6N6O6 + 6*C10F6Cl6H9'),
('PBXN', 'C50H80N57O57'),
('PBXN301', '0.8*C5H8N4O12+0.2*C2H6O1Si1'),
('TATB', 'C6H6N6O6'),
('HMX', 'C4H8N8O8'),
('LX17', '92.5*C6H6N6O6+7.5*C2ClF3'),
('TNT', 'C7H5N3O6'),
('GOS', 'Gd2O2S'),
('CsI', 'CsI'),
('LSO', 'Lu2SiO5'),
('LYSO', 'Lu1.8Y0.2SiO5'),
('NaI', 'NaI'),
('BGO', 'Bi4Ge3O12'),
('GSO', 'Gd2SiO5'),
('CZT', 'Cd0.9Zn0.1Te')])

# g/cm^3
materialDensities = dict([('vacuum', 1.0e-16),
('air', 1.2e-3),
('foam', 0.035),
('water', 1.0),
('bone_cortical', 1.85),
('bone_compact',1.85),
('adipose_tissue', 0.92),
('lung', 0.25),
('ethanol', 0.789),
('HP', 1.450),
('saltWater', 1.025),
('graphite', 1.804),
('diamond', 3.51),
('paper', 0.8),
('rubber', 1.1),
('wood', 0.53),
('delrin', 1.3962),
('teflon', 2.1609),
('PMMA', 1.18),
('acrylic', 1.18),
('Polyethylene', 0.94),
('PE', 0.94),
('LDPE', 0.92),
('HDPE', 0.95),
('Polypropylene', 0.946),
('PP', 0.946),
('PET', 1.38),
('Polycarbonate', 1.21),
('PVC', 1.35),
('PVDF', 1.78),
('melamine', 1.57),
('PETG', 1.4),
('PPSU', 1.28),
('PEEK', 1.32),
('silicone', 1.159604139),
('glass', 2.53),
('eglass', 2.6),
('aluminumOxide', 4.0),
('titaniumDioxide', 4.23),
('brass', 8.73),
('Inconel', 8.1933),
('stainlessSteel', 8.0),
('Ti64', 4.429),
('Ti5553', 4.65),
('304L', 8.0),
('316L', 8.0),
('Al6061', 2.7),
('aermet100', 7.888),
('NM', 1.1371),
('PTM', 0.9),
('RDX', 1.63),
('PBXN', 1.65),
('PBXN301', 1.5615),
('TATB', 1.93),
('HMX', 1.91),
('LX17', 1.9),
('TNT', 1.65),
('GOS', 7.32),
('CsI', 4.51),
('LSO', 7.4),
('LYSO', 7.1),
('NaI', 3.67),
('BGO', 7.13),
('GSO', 6.71),
('CZT', 5.78)])
