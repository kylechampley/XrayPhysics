#include "xrayphysics.h"
#include <string>

XrayPhysics::XrayPhysics()
{
    XraySourceModel.init(&xsecTables);
}

XrayPhysics::~XrayPhysics()
{
}

bool XrayPhysics::simulateSpectra(float kVp, float takeOffAngle, int Z, float* gammas, int N, float* output)
{
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

bool XrayPhysics::changeTakeOffAngle(float kVp, float takeOffAngle_cur, float takeOffAngle_new, int Z, float* gammas, int N, float* s)
{
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

float XrayPhysics::meanEnergy(float* spectralResponse, float* gammas, int N)
{
    if (spectralResponse == NULL || gammas == NULL || N <= 0)
        return 0.0;
    if (N == 1)
        return gammas[0];
    double retVal = 0.0;
    double accum = 0.0;
    for (int i = 0; i < N; i++)
    {
        double T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        retVal += gammas[i] * T_phi * spectralResponse[i];
        accum += T_phi * spectralResponse[i];
    }
    return float(retVal / accum);
}

float XrayPhysics::effectiveAttenuation(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    if (spectralResponse == NULL || gammas == NULL || N <= 0)
        return 0.0;
    if (N == 1)
        return density*xsecTables.sigma(Z, gammas[0]);

    double retVal = 0.0;
    double accum = 0.0;
    for (int i = 0; i < N; i++)
    {
        float T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        if (thickness == 0.0)
        {
            retVal += T_phi * spectralResponse[i] * density * xsecTables.sigma(Z, gammas[i]);
        }
        else
        {
            retVal += T_phi * spectralResponse[i] * exp(-density * thickness * xsecTables.sigma(Z, gammas[i]));
        }
        accum += T_phi * spectralResponse[i];
    }
    retVal = retVal / accum;

    if (thickness != 0.0)
        retVal = -log(retVal) / thickness;
    return float(retVal);
}

float XrayPhysics::effectiveEnergy(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    float mu_eff = effectiveAttenuation(Z, density, thickness, spectralResponse, gammas, N);
    return xsecTables.sigma_inv(Z, mu_eff /density);
}

float XrayPhysics::effectiveAttenuation(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    if (spectralResponse == NULL || gammas == NULL || N <= 0)
        return 0.0;
    if (N == 1)
        return density * xsecTables.sigma(chemForm, gammas[0]);

    double retVal = 0.0;
    double accum = 0.0;
    for (int i = 0; i < N; i++)
    {
        float T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        if (thickness == 0.0)
        {
            retVal += T_phi * spectralResponse[i] * density * xsecTables.sigma(chemForm, gammas[i]);
        }
        else
        {
            retVal += T_phi * spectralResponse[i] * exp(-density * thickness * xsecTables.sigma(chemForm, gammas[i]));
        }
        accum += T_phi * spectralResponse[i];
    }
    retVal = retVal / accum;

    if (thickness != 0.0)
        retVal = -log(retVal) / thickness;
    return float(retVal);
}

float XrayPhysics::effectiveEnergy(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    float mu_eff = effectiveAttenuation(chemForm, density, thickness, spectralResponse, gammas, N);
    return xsecTables.sigma_inv(chemForm, mu_eff / density);
}

float XrayPhysics::transmission(float Z, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    if (spectralResponse == NULL || gammas == NULL || N <= 0)
        return 1.0;
    if (N == 1)
        return exp(-thickness * density * xsecTables.sigma(Z, gammas[0]));

    double retVal = 0.0;
    double accum = 0.0;
    for (int i = 0; i < N; i++)
    {
        float T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        if (thickness == 0.0)
            retVal += T_phi * spectralResponse[i] * density * xsecTables.sigma(Z, gammas[i]);
        else
            retVal += T_phi * spectralResponse[i] * exp(-density * thickness * xsecTables.sigma(Z, gammas[i]));
        accum += T_phi * spectralResponse[i];
    }
    retVal = retVal / accum;

    return float(retVal);
}

float XrayPhysics::transmission(const char* chemForm, float density, float thickness, float* spectralResponse, float* gammas, int N)
{
    if (spectralResponse == NULL || gammas == NULL || N <= 0)
        return 1.0;
    if (N == 1)
        return exp(-thickness * density * xsecTables.sigma(chemForm, gammas[0]));

    double retVal = 0.0;
    double accum = 0.0;
    for (int i = 0; i < N; i++)
    {
        float T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        if (thickness == 0.0)
            retVal += T_phi * spectralResponse[i] * density * xsecTables.sigma(chemForm, gammas[i]);
        else
            retVal += T_phi * spectralResponse[i] * exp(-density * thickness * xsecTables.sigma(chemForm, gammas[i]));
        accum += T_phi * spectralResponse[i];
    }
    retVal = retVal / accum;

    return float(retVal);
}

bool XrayPhysics::setBHClookupTable_helper(double* sigma_hat, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    if (spectralResponse == NULL || gammas == NULL || N_gamma <= 0)
        return false;
    if (LUT == NULL || T_lac <= 0.0 || N_lac <= 0)
        return false;

    double* d = new double[N_gamma];
    double accum = 0.0;
    for (int i = 0; i < N_gamma; i++)
    {
        double T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N_gamma - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        d[i] = T_phi * spectralResponse[i];
        accum += T_phi * spectralResponse[i];
    }
    for (int i = 0; i < N_gamma; i++)
        d[i] = d[i] / accum;

    bool retVal = true;

    int N_iter = 10;
    double tol = 1.0e-7;

    LUT[0] = 0.0;
    for (int i = 1; i < N_lac; i++)
    {
        double curLAC = double(i) * T_lac;
        double p_cur = LUT[i - 1];
        if (i == 1)
            p_cur = curLAC;
        double modelError, denom, theExponent, expFactor;
        for (int n = 0; n < N_iter; n++)
        {
            modelError = 0.0;
            denom = 0.0;
            for (int l = 0; l < N_gamma; l++)
            {
                theExponent = sigma_hat[l];
                expFactor = exp(-theExponent * p_cur);
                modelError += d[l] * expFactor;
                denom += d[l] * theExponent * expFactor;
            }
            //modelError *= T_gamma;
            //denom *= T_gamma;

            denom = denom / modelError; // new method
            modelError = (curLAC + log(modelError)); // new method
            if (fabs(modelError) < tol || fabs(denom) < 1.0e-15)
                break;
            p_cur = p_cur + modelError / denom;
        }
        if (isnan(p_cur))
        {
            printf("XrayPhysics::setBHClookupTable: nan value encountered, quitting!\n");
            retVal = false;
            break;
        }
        if (p_cur >= 0.0)
            LUT[i] = p_cur;
        else
            LUT[i] = LUT[i - 1];
        //printf("%f => %f\n", curLAC, LUT[i]);
    }
    delete[] d;
    return retVal;
}

bool XrayPhysics::setBHClookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    if (spectralResponse == NULL || gammas == NULL || N_gamma <= 0)
        return false;
    if (LUT == NULL || T_lac <= 0.0 || N_lac <= 0)
        return false;

    double* sigma_hat = new double[N_gamma];
    double sigma_refE = xsecTables.sigma(Ze, referenceEnergy);
    for (int i = 0; i < N_gamma; i++)
        sigma_hat[i] = xsecTables.sigma(Ze, gammas[i]) / sigma_refE;

    bool retVal = setBHClookupTable_helper(sigma_hat, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);

    delete[] sigma_hat;
    return retVal;
}

bool XrayPhysics::setBHClookupTable(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    if (spectralResponse == NULL || gammas == NULL || N_gamma <= 0)
        return false;
    if (LUT == NULL || T_lac <= 0.0 || N_lac <= 0)
        return false;

    double* sigma_hat = new double[N_gamma];
    double sigma_refE = xsecTables.sigma(chemForm, referenceEnergy);
    for (int i = 0; i < N_gamma; i++)
        sigma_hat[i] = xsecTables.sigma(chemForm, gammas[i]) / sigma_refE;

    //printf("referenceEnergy = %f, sigma_refE = %f\n", referenceEnergy, sigma_refE);

    bool retVal = setBHClookupTable_helper(sigma_hat, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);

    delete[] sigma_hat;
    return retVal;
}

bool XrayPhysics::setBHlookupTable(float Ze, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    if (spectralResponse == NULL || gammas == NULL || N_gamma <= 0)
        return false;
    if (LUT == NULL || T_lac <= 0.0 || N_lac <= 0)
        return false;

    double* sigma_hat = new double[N_gamma];
    double sigma_refE = xsecTables.sigma(Ze, referenceEnergy);
    for (int i = 0; i < N_gamma; i++)
        sigma_hat[i] = xsecTables.sigma(Ze, gammas[i]) / sigma_refE;

    bool retVal = setBHlookupTable_helper(sigma_hat, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);

    delete[] sigma_hat;
    return retVal;
}

bool XrayPhysics::setBHlookupTable(const char* chemForm, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    if (spectralResponse == NULL || gammas == NULL || N_gamma <= 0)
        return false;
    if (LUT == NULL || T_lac <= 0.0 || N_lac <= 0)
        return false;

    double* sigma_hat = new double[N_gamma];
    double sigma_refE = xsecTables.sigma(chemForm, referenceEnergy);
    for (int i = 0; i < N_gamma; i++)
        sigma_hat[i] = xsecTables.sigma(chemForm, gammas[i]) / sigma_refE;

    //printf("referenceEnergy = %f, sigma_refE = %f\n", referenceEnergy, sigma_refE);

    bool retVal = setBHlookupTable_helper(sigma_hat, spectralResponse, gammas, N_gamma, LUT, T_lac, N_lac, referenceEnergy);

    delete[] sigma_hat;
    return retVal;
}

bool XrayPhysics::setBHlookupTable_helper(double* sigma_hat, float* spectralResponse, float* gammas, int N_gamma, float* LUT, float T_lac, int N_lac, float referenceEnergy)
{
    if (spectralResponse == NULL || gammas == NULL || N_gamma <= 0)
        return false;
    if (LUT == NULL || T_lac <= 0.0 || N_lac <= 0)
        return false;

    double* d = new double[N_gamma];
    double accum = 0.0;
    for (int i = 0; i < N_gamma; i++)
    {
        double T_phi;
        if (i == 0)
            T_phi = gammas[i + 1] - gammas[i];
        else if (i == N_gamma - 1)
            T_phi = gammas[i] - gammas[i - 1];
        else
            T_phi = 0.5 * (gammas[i + 1] - gammas[i - 1]);

        d[i] = T_phi * spectralResponse[i];
        accum += T_phi * spectralResponse[i];
    }
    for (int i = 0; i < N_gamma; i++)
        d[i] = d[i] / accum;

    bool retVal = true;

    int N_iter = 10;
    double tol = 1.0e-7;

    LUT[0] = 0.0;
    for (int i = 1; i < N_lac; i++)
    {
        double curLAC = double(i) * T_lac;
        double accum = 0.0;
        for (int l = 0; l < N_gamma; l++)
        {
            accum += d[l]*exp(-sigma_hat[l] * curLAC);
        }
        LUT[i] = -log(accum);
        //printf("%f => %f\n", curLAC, LUT[i]);
    }
    delete[] d;
    return retVal;
}
