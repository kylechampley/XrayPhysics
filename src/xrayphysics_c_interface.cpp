
#include "xrayphysics_c_interface.h"
#include "xsec.h"
#include "xsource.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

xsec xsecTables;
//xsecTables.init();

bool simulateSpectra(float kVp, float takeOffAngle, int Z, float* gammas, int N, float* output)
{
    xraySource XraySourceModel;
    XraySourceModel.init(&xsecTables);
    float* s = XraySourceModel.simulateSpectra(kVp, takeOffAngle, Z, gammas, N);
    if (s != NULL)
    {
        for (int i = 0; i < N; i++)
            output[i] = s[i];
        free(s);
        return true;
    }
    else
        return false;    
}

bool changeTakeOffAngle(float kVp, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s)
{
    xraySource XraySourceModel;
    XraySourceModel.init(&xsecTables);
    float* factors = XraySourceModel.takeOffAngleConversionFactor(kVp, takeOffAngle_cur, takeOffAngle_new, Z, gammas, N);
    if (factors != NULL)
    {
        for (int i = 0; i < N; i++)
            s[i] *= factors[i];
        free(factors);
        return true;
    }
    else
        return false;

}

float sigma(float Z, float gamma)
{
    return xsecTables.sigma(Z, gamma);
}

float sigmaCompound(const char* chemForm, float gamma)
{
    return xsecTables.sigma(chemForm, gamma);
}

float sigmaPE(float Z, float gamma)
{
    return xsecTables.sigmaPE(Z, gamma);
}

float sigmaCompoundPE(const char* chemForm, float gamma)
{
    return xsecTables.sigmaPE(chemForm, gamma);
}

float sigmaCS(float Z, float gamma)
{
    return xsecTables.sigmaCS(Z, gamma);
}

float sigmaCompoundCS(const char* chemForm, float gamma)
{
    return xsecTables.sigmaCS(chemForm, gamma);
}

float sigmaRS(float Z, float gamma)
{
    return xsecTables.sigmaRS(Z, gamma);
}

float sigmaCompoundRS(const char* chemForm, float gamma)
{
    return xsecTables.sigmaRS(chemForm, gamma);
}

float sigmaPP(float Z, float gamma)
{
    return xsecTables.sigmaPP(Z, gamma);
}

float sigmaCompoundPP(const char* chemForm, float gamma)
{
    return xsecTables.sigmaPP(chemForm, gamma);
}

float sigmaTP(float Z, float gamma)
{
    return xsecTables.sigmaTP(Z, gamma);
}

float sigmaCompoundTP(const char* chemForm, float gamma)
{
    return xsecTables.sigmaTP(chemForm, gamma);
}
