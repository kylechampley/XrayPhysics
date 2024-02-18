#ifndef XRAYPHYSICS_H
#define XRAYPHYSICS_H

#ifdef WIN32
#pragma once
#endif

#define XRAY_PHYSICS_VERSION "1.1"

#include "xsec.h"
#include "xscatter.h"
#include "xsource.h"

class XrayPhysics
{
public:
    XrayPhysics();
    ~XrayPhysics();

    const char* about();

    float atomicMass(int Z);

    bool simulateSpectra(float kVp, float takeOffAngle, int Z, float* gammas, int N, float* output);
    bool changeTakeOffAngle(float kVp, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s);

    float meanEnergy(float* spectralResponse, float* gammas, int N);
    bool normalizeSpectrum(float* spectralResponse, float* gammas, int N);

    float effectiveAttenuation(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
    float effectiveEnergy(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);

    float effectiveAttenuation(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N);
    float effectiveEnergy(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N);

    float transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
    float transmission(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N);

    float effectiveZ(const char* chemForm, float min_energy, float max_energy, float arealDensity = 0.0);

    float incoherentScatterDistribution(float Z, float gamma, float theta);
    float coherentScatterDistribution(float Z, float gamma, float theta);

    float incoherentScatterDistribution(const char* chemForm, float gamma, float theta);
    float coherentScatterDistribution(const char* chemForm, float gamma, float theta);

    float incoherentScatterDistribution_normalizationFactor(float Z, float gamma);
    float coherentScatterDistribution_normalizationFactor(float Z, float gamma);

    float incoherentScatterDistribution_normalizationFactor(const char* chemForm, float gamma);
    float coherentScatterDistribution_normalizationFactor(const char* chemForm, float gamma);

    bool setBHlookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
    bool setBHlookupTable(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);

    bool setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
    bool setBHClookupTable(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);

    bool generateDEDlookUpTables(float* spectralResponses, float* gammas, int N_gamma, float* referenceEnergies, float* basisFunctions, float* LUT, float T_lac, int N_lac);

    bool setTwoMaterialBHClookupTable(float* spectralResponse, float* gammas, int N_gamma, float referenceEnergy, float* sigmas, float* LUT, float T_atten, int N_atten);

    xsec xsecTables;
    xscatter xscatterTables;
    xraySource XraySourceModel;
private:

    bool setBHlookupTable_helper(double* sigma_hat, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
    bool setBHClookupTable_helper(double* sigma_hat, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);

    bool BHCkernel(double& monoAtten, double polyAtten, double* normalizedCrossSection, double* d, int N_gamma);
    float gamma_inv(float val, float* gammas, int N_gamma);
};

//BHC and BH
//pBHC

#endif
