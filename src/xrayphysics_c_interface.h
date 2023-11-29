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
