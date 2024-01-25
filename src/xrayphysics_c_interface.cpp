
#include "xrayphysics_c_interface.h"
//#include "xsec.h"
//#include "xsource.h"
#include "xrayphysics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

XrayPhysics physics;
//xsecTables.init();

bool simulateSpectra(float kVp, float takeOffAngle, int Z, float* gammas, int N, float* output)
{
    return physics.simulateSpectra(kVp, takeOffAngle, Z, gammas, N, output);
}

bool changeTakeOffAngle(float kVp, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s)
{
    return physics.changeTakeOffAngle(kVp, takeOffAngle_cur, takeOffAngle_new, Z, gammas, N, s);
}

float sigma(float Z, float gamma)
{
    return physics.xsecTables.sigma(Z, gamma);
}

float sigmaCompound(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigma(chemForm, gamma);
}

float sigmae(float Z, float gamma)
{
    return physics.xsecTables.sigma_e(Z, gamma);
}

float sigmaeCompound(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigma_e(chemForm, gamma);
}

float sigmaPE(float Z, float gamma)
{
    return physics.xsecTables.sigmaPE(Z, gamma);
}

float sigmaCompoundPE(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigmaPE(chemForm, gamma);
}

float sigmaCS(float Z, float gamma)
{
    return physics.xsecTables.sigmaCS(Z, gamma);
}

float sigmaCompoundCS(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigmaCS(chemForm, gamma);
}

float sigmaRS(float Z, float gamma)
{
    return physics.xsecTables.sigmaRS(Z, gamma);
}

float sigmaCompoundRS(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigmaRS(chemForm, gamma);
}

float sigmaPP(float Z, float gamma)
{
    return physics.xsecTables.sigmaPP(Z, gamma);
}

float sigmaCompoundPP(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigmaPP(chemForm, gamma);
}

float sigmaTP(float Z, float gamma)
{
    return physics.xsecTables.sigmaTP(Z, gamma);
}

float sigmaCompoundTP(const char* chemForm, float gamma)
{
    return physics.xsecTables.sigmaTP(chemForm, gamma);
}

float meanEnergy(float* spectralResponse, float* gammas, int N)
{
    return physics.meanEnergy(spectralResponse, gammas, N);
}

bool normalizeSpectrum(float* spectralResponse, float* gammas, int N)
{
    return physics.normalizeSpectrum(spectralResponse, gammas, N);
}

float effectiveZ(const char* chemForm, float min_energy, float max_energy, float arealDensity)
{
    return physics.effectiveZ(chemForm, min_energy, max_energy, arealDensity);
}

float effectiveAttenuation(float Ze, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    return physics.effectiveAttenuation(Ze, density, thickness, spectralResponse, gammas, N);
}

float effectiveEnergy(float Ze, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    return physics.effectiveEnergy(Ze, density, thickness, spectralResponse, gammas, N);
}

float effectiveAttenuation_compound(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    return physics.effectiveAttenuation(chemForm, density, thickness, spectralResponse, gammas, N);
}

float effectiveEnergy_compound(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    return physics.effectiveEnergy(chemForm, density, thickness, spectralResponse, gammas, N);
}

float transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    return physics.transmission(Z, density, thickness, spectralResponse, gammas, N);
}

float transmission_compound(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    return physics.transmission(chemForm, density, thickness, spectralResponse, gammas, N);
}

bool setBHlookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    return physics.setBHlookupTable(Ze, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);
}

bool setBHlookupTable_compound(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    return physics.setBHlookupTable(chemForm, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);
}

bool setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    return physics.setBHClookupTable(Ze, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);
}

bool setBHClookupTable_compound(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    return physics.setBHClookupTable(chemForm, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);
}

bool generateDEDlookUpTables(float* spectralResponses, float* gammas, int N_gamma, float* referenceEnergies, float* basisFunctions, float* LUT, float T_lac, int N_lac)
{
    return physics.generateDEDlookUpTables(spectralResponses, gammas, N_gamma, referenceEnergies, basisFunctions, LUT, T_lac, N_lac);
}
