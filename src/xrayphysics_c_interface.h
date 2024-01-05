#ifdef WIN32
    #pragma once

    #ifdef PROJECTOR_EXPORTS
        #define XRAYPHYSICS_API __declspec(dllexport)
    #else
        #define XRAYPHYSICS_API __declspec(dllimport)
    #endif
#else
    #define XRAYPHYSICS_API
#endif

extern "C" XRAYPHYSICS_API bool simulateSpectra(float kVp, float takeOffAngle, int Z, float* gammas, int N, float* output);
extern "C" XRAYPHYSICS_API bool changeTakeOffAngle(float kVp, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s);

extern "C" XRAYPHYSICS_API float sigma(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaCompound(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float sigmae(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaeCompound(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float sigmaPE(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaCompoundPE(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float sigmaCS(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaCompoundCS(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float sigmaRS(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaCompoundRS(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float sigmaPP(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaCompoundPP(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float sigmaTP(float Z, float gamma);
extern "C" XRAYPHYSICS_API float sigmaCompoundTP(const char* chemForm, float gamma);

extern "C" XRAYPHYSICS_API float meanEnergy(float* spectralResponse, float* gammas, int N);

extern "C" XRAYPHYSICS_API float effectiveZ(const char* chemForm, float min_energy, float max_energy, float arealDensity);

extern "C" XRAYPHYSICS_API float effectiveAttenuation(float Ze, float density, float thickness, float* spectralResponse, float* gammas, int N);
extern "C" XRAYPHYSICS_API float effectiveEnergy(float Ze, float density, float thickness, float* spectralResponse, float* gammas, int N);

extern "C" XRAYPHYSICS_API float effectiveAttenuation_compound(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N);
extern "C" XRAYPHYSICS_API float effectiveEnergy_compound(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N);

extern "C" XRAYPHYSICS_API float transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N);
extern "C" XRAYPHYSICS_API float transmission_compound(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N);

extern "C" XRAYPHYSICS_API bool setBHlookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
extern "C" XRAYPHYSICS_API bool setBHlookupTable_compound(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);

extern "C" XRAYPHYSICS_API bool setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
extern "C" XRAYPHYSICS_API bool setBHClookupTable_compound(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy);
