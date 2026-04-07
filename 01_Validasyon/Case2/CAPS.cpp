#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

// =============================================================================
// LOOP DIŞI SKALER DEĞİŞKENLER
// =============================================================================

// blowing_factor
double bf_phi, bf_Omega;

// piroliz
double pyr_I, pyr_ai, pyr_anew, pyr_bkt, pyr_exp;

// basinc denklemi
double p_dxl, p_dxr, p_dxi;
double p_phi, p_atime, p_K, p_ae, p_aw;
double p_aeff_old, p_dalpha, p_Sptemp, p_Spporo, p_Slin;

// sicaklik denklemi
double t_dxl, t_dxr, t_dxi;
double t_cp, t_phi, t_ke, t_kw;
double t_rhogas, t_rhoc, t_atime, t_ae, t_aw;
double t_mdotgas, t_hgrad;
double t_hsolid, t_hpyro, t_Spyro, t_Stotal;

// =============================================================================
// küçük yardımcılar
// =============================================================================

double inverseAverage(double v1, double d1, double v2, double d2)
{
    d1 = (fabs(d1));
    d2 = (fabs(d2));
    return (v1 * d2 + v2 * d1) / (d1 + d2);
}

// =============================================================================
// BLOWING FACTOR  (Case 5: enthalpy-based)
// phi = 2*lambda * m_dot / (rho_ue_CH0)
// Omega = phi / (exp(phi) - 1)
// h_eff_out = Omega * rho_ue_CH   [kg/m2/s * J/kg = W/m2, enthalpy-based]
// =============================================================================
double blowing_factor(double m_dot, double rho_ue_CH0, double& h_eff_out) //BU KISIM İNCELENMELİ DİKKATLİCE
{
    bf_phi = 2.0 * 0.5 * m_dot / rho_ue_CH0;

    if (fabs(bf_phi) < 1e-8)
        bf_Omega = 1.0 - 0.5 * bf_phi;
    else
        bf_Omega = bf_phi / (exp(bf_phi) - 1.0);

    h_eff_out = bf_Omega * rho_ue_CH0;   // blowing-corrected rho_ue_CH [kg/m2/s]
    return bf_Omega;
}

// =============================================================================
// PATANKAR TDMA
// =============================================================================
vector<double> thomas_patankar(const vector<double>& a,
                              const vector<double>& b,
                              const vector<double>& c,
                              const vector<double>& d)
{
    int N = (int)a.size();
    vector<double> P(N), Q(N), T(N);

    P[0] = b[0] / a[0];
    Q[0] = d[0] / a[0];
    for (int i = 1; i < N; i++)
    {
        double den = a[i] - c[i] * P[i-1];
        P[i] = b[i] / den;
        Q[i] = (d[i] + c[i] * Q[i-1]) / den;
    }
    T[N-1] = Q[N-1];
    for (int i = N-2; i >= 0; i--)
        T[i] = P[i] * T[i+1] + Q[i];
    return T;
}

// =============================================================================
// TERMOFIZIKSEL TABLOLAR
// =============================================================================

static inline double lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}

double interp1_linear(const vector<double>& x, const vector<double>& y, double xq)
{
    const int N = (int)x.size();
    if (N == 0) return 0.0;
    if (N == 1) return y[0];

    if (xq <= x[0])     return y[0];
    if (xq >= x[N - 1]) return y[N - 1];

    int i = 0;
    while (i < N - 1 && xq > x[i + 1]) i++;

    const double t = (xq - x[i]) / (x[i + 1] - x[i]);
    return lerp(y[i], y[i + 1], t);
}

struct ThermalTable
{
    vector<double> T;
    vector<double> cp;
    vector<double> k;
    vector<double> eps;
};

ThermalTable tbl_v, tbl_c;

static void push_row_thermal(ThermalTable& tbl, double T, double cp, double k, double eps)
{
    tbl.T.push_back(T);
    tbl.cp.push_back(cp);
    tbl.k.push_back(k);
    tbl.eps.push_back(eps);
}

void init_virgin_table()
{
    tbl_v = ThermalTable{};
    push_row_thermal(tbl_v, 256,  879.2, 0.3975, 0.8);
    push_row_thermal(tbl_v, 298,  983.9, 0.4025, 0.8);
    push_row_thermal(tbl_v, 444, 1298, 0.4162, 0.8);
    push_row_thermal(tbl_v, 556, 1465, 0.453, 0.8);
    push_row_thermal(tbl_v, 644, 1570, 0.4698, 0.8);
    push_row_thermal(tbl_v, 833, 1717, 0.486, 0.8);
    push_row_thermal(tbl_v, 1111,1863, 0.5234, 0.8);
    push_row_thermal(tbl_v, 1389,1934, 0.5601, 0.8);
    push_row_thermal(tbl_v, 1667,1980, 0.6978, 0.8);
    push_row_thermal(tbl_v, 1944,1989, 0.8723, 0.8);
    push_row_thermal(tbl_v, 2222,2001, 1.109, 0.8);
    push_row_thermal(tbl_v, 2778,2010, 1.751, 0.8);
    push_row_thermal(tbl_v, 3333,2010, 2.779, 0.8);
}

void init_char_table()
{
    tbl_c = ThermalTable{};
    push_row_thermal(tbl_c, 256,  7327, 0.3975, 0.9);
    push_row_thermal(tbl_c, 298,  7829, 0.4025, 0.9);
    push_row_thermal(tbl_c, 444, 1093, 0.4162, 0.9);
    push_row_thermal(tbl_c, 556, 1319, 0.453, 0.9);
    push_row_thermal(tbl_c, 644, 1432, 0.4698, 0.9);
    push_row_thermal(tbl_c, 833, 1675, 0.486, 0.9);
    push_row_thermal(tbl_c, 1111,1842, 0.5234, 0.9);
    push_row_thermal(tbl_c, 1389,1968, 0.5601, 0.9);
    push_row_thermal(tbl_c, 1667,2052, 0.605, 0.9);
    push_row_thermal(tbl_c, 1944,2093, 0.729, 0.9);
    push_row_thermal(tbl_c, 2222,2110, 0.9221, 0.9);
    push_row_thermal(tbl_c, 2778,2135, 1.458, 0.9);
    push_row_thermal(tbl_c, 3333,2152, 2.318, 0.9);
}

double cp_v_T(double T) { return interp1_linear(tbl_v.T, tbl_v.cp, T); }
double k_v_T (double T) { return interp1_linear(tbl_v.T, tbl_v.k,  T); }
double cp_c_T(double T) { return interp1_linear(tbl_c.T, tbl_c.cp, T); }
double k_c_T (double T) { return interp1_linear(tbl_c.T, tbl_c.k,  T); }

double cp_mix(double T, double a) { return cp_v_T(T)*(1-a) + cp_c_T(T)*a; }
double k_mix (double T, double a) { return k_v_T(T) *(1-a) + k_c_T(T) *(a); }

double eps_surf(double T, double a)
{
    const double eps_v = interp1_linear(tbl_v.T, tbl_v.eps, T);
    const double eps_c = interp1_linear(tbl_c.T, tbl_c.eps, T);
    return eps_v*(1.0-a) + eps_c*a;
}

// =============================================================================
// PIROLIZ GAZI TABLOLARI
// =============================================================================

struct PyroGasTable
{
    vector<double> T;
    vector<double> cp;   // J/kg/K
    vector<double> mu;   // Pa.s
    vector<double> h;    // J/kg
    vector<double> MW;   // kg/kmol
};

PyroGasTable tbl_pg;

static void push_row_pyro_gas(PyroGasTable& tbl, double T, double cp, double mu)
{
    tbl.T.push_back(T);
    tbl.cp.push_back(cp);
    tbl.mu.push_back(mu);
}
void init_pyrogas_table()
{
    tbl_pg = PyroGasTable{};
    // pyrolysisModel.csv — P = 1.0 atm kesiti
    // kolonlar: T(K), h(kJ/kg), mu(millipoise→Pa·s), cp(kJ/kg·K)
    // push_row_pyro_gas(tbl_pg, T, cp [J/kgK], mu [Pa·s], h [J/kg])
    auto push = [&](double T, double cp_kJ, double mu_mp, double h_kJ, double MW_kgkmol){
        tbl_pg.T.push_back(T);
        tbl_pg.cp.push_back(cp_kJ * 1000.0);   // kJ/kgK -> J/kgK
        tbl_pg.mu.push_back(mu_mp * 1e-4);     // millipoise -> Pa.s
        tbl_pg.h.push_back(h_kJ * 1000.0);     // kJ/kg -> J/kg
        tbl_pg.MW.push_back(MW_kgkmol);        // kg/kmol
    };
        // push(200, 1.512, 8.688e-02, -7247, 22.00);
        // push(225, 1.534, 9.666e-02, -7208, 22.00);
        // push(250, 1.560, 1.065e-01, -7170, 22.00);
        // push(275, 1.592, 1.162e-01, -7130, 22.00);
        // push(300, 1.631, 1.257e-01, -7090, 22.00);
        // push(325, 1.676, 1.351e-01, -7049, 22.00);
        // push(350, 1.726, 1.444e-01, -7006, 22.00);
        // push(375, 1.783, 1.534e-01, -6963, 21.99);
        // push(400, 1.847, 1.623e-01, -6917, 21.99);
        // push(425, 1.921, 1.710e-01, -6870, 21.99);
        // push(450, 2.008, 1.796e-01, -6821, 21.98);
        // push(475, 2.112, 1.879e-01, -6770, 21.97);
        // push(500, 2.241, 1.962e-01, -6715, 21.95);
        // push(525, 2.402, 2.042e-01, -6657, 21.92);
        // push(550, 2.605, 2.122e-01, -6595, 21.87);
        // push(575, 2.867, 2.201e-01, -6527, 21.80);
        // push(600, 3.210, 2.279e-01, -6451, 21.71);
        // push(625, 3.671, 2.356e-01, -6365, 21.59);
        // push(650, 4.301, 2.432e-01, -6266, 21.42);
        // push(675, 5.171, 2.509e-01, -6148, 21.19);
        // push(700, 6.351, 2.586e-01, -6005, 20.89);
        // push(725, 7.882, 2.664e-01, -5827, 20.50);
        // push(750, 9.748, 2.744e-01, -5608, 19.99);
        // push(775, 11.85, 2.826e-01, -5338, 19.37);
        // push(800, 14.03, 2.909e-01, -5014, 18.64);
        // push(825, 16.01, 2.993e-01, -4638, 17.84);
        // push(850, 17.44, 3.076e-01, -4219, 17.00);
        // push(875, 17.89, 3.157e-01, -3775, 16.19);
        // push(900, 17.01, 3.235e-01, -3335, 15.46);
        // push(925, 14.77, 3.309e-01, -2936, 14.86);
        // push(950, 11.65, 3.378e-01, -2605, 14.41);
        // push(975, 8.558, 3.444e-01, -2353, 14.12);
        // push(1000, 6.222, 3.506e-01, -2170, 13.95);
        // push(1025, 4.784, 3.566e-01, -2034, 13.85);
        // push(1050, 4.014, 3.625e-01, -1925, 13.80);
        // push(1075, 3.644, 3.683e-01, -1830, 13.78);
        // push(1100, 3.509, 3.740e-01, -1741, 13.76);
        // push(1125, 3.561, 3.797e-01, -1653, 13.75);
        // push(1150, 3.901, 3.854e-01, -1561, 13.74);
        // push(1175, 4.807, 3.913e-01, -1454, 13.71);
        // push(1200, 6.235, 3.979e-01, -1316, 13.64);
        // push(1225, 7.447, 4.053e-01, -1144, 13.53);
        // push(1250, 8.140, 4.131e-01, -948.0, 13.40);
        // push(1275, 8.479, 4.212e-01, -739.7, 13.26);
        // push(1300, 8.634, 4.292e-01, -525.5, 13.11);
        // push(1325, 8.707, 4.371e-01, -308.7, 12.97);
        // push(1350, 8.767, 4.448e-01, -90.28, 12.84);
        // push(1375, 8.860, 4.521e-01, 129.9, 12.71);
        // push(1400, 9.024, 4.589e-01, 353.3, 12.58);
        // push(1425, 9.288, 4.651e-01, 582.0, 12.46);
        // push(1450, 9.674, 4.708e-01, 818.7, 12.34);
        // push(1475, 10.19, 4.757e-01, 1067, 12.22);
        // push(1500, 10.82, 4.799e-01, 1329, 12.10);
        // push(1525, 11.52, 4.832e-01, 1608, 11.98);
        // push(1550, 12.15, 4.859e-01, 1905, 11.86);
        // push(1575, 12.53, 4.880e-01, 2214, 11.73);
        // push(1600, 12.38, 4.899e-01, 2527, 11.61);
        // push(1625, 11.51, 4.920e-01, 2827, 11.50);
        // push(1650, 10.10, 4.946e-01, 3098, 11.40);
        // push(1675, 8.590, 4.981e-01, 3331, 11.32);
        // push(1700, 7.338, 5.021e-01, 3529, 11.26);
        // push(1725, 6.422, 5.066e-01, 3701, 11.21);
        // push(1750, 5.774, 5.112e-01, 3853, 11.17);
        // push(1775, 5.312, 5.160e-01, 3991, 11.14);
        // push(1800, 4.974, 5.209e-01, 4119, 11.12);
        // push(1825, 4.723, 5.258e-01, 4241, 11.10);
        // push(1850, 4.533, 5.306e-01, 4356, 11.08);
        // push(1875, 4.389, 5.355e-01, 4468, 11.07);
        // push(1900, 4.281, 5.404e-01, 4576, 11.06);
        // push(1925, 4.200, 5.453e-01, 4682, 11.05);
        // push(1950, 4.142, 5.502e-01, 4786, 11.04);
        // push(1975, 4.103, 5.550e-01, 4889, 11.03);
        // push(2000, 4.078, 5.599e-01, 4991, 11.02);
        // push(2025, 4.068, 5.647e-01, 5093, 11.02);
        // push(2050, 4.068, 5.696e-01, 5195, 11.01);
        // push(2075, 4.080, 5.744e-01, 5297, 11.01);
        // push(2100, 4.101, 5.792e-01, 5399, 11.00);
        // push(2125, 4.130, 5.840e-01, 5502, 11.00);
        // push(2150, 4.169, 5.888e-01, 5605, 11.00);
        // push(2175, 4.215, 5.936e-01, 5710, 10.99);
        // push(2200, 4.270, 5.984e-01, 5816, 10.99);
        // push(2225, 4.333, 6.032e-01, 5924, 10.98);
        // push(2250, 4.405, 6.080e-01, 6033, 10.98);
        // push(2275, 4.484, 6.128e-01, 6144, 10.97);
        // push(2300, 4.573, 6.176e-01, 6257, 10.96);
        // push(2325, 4.670, 6.224e-01, 6373, 10.96);
        // push(2350, 4.777, 6.272e-01, 6491, 10.95);
        // push(2375, 4.893, 6.320e-01, 6612, 10.94);
        // push(2400, 5.020, 6.368e-01, 6736, 10.93);
        // push(2425, 5.157, 6.416e-01, 6863, 10.92);
        // push(2450, 5.305, 6.464e-01, 6994, 10.91);
        // push(2475, 5.465, 6.512e-01, 7128, 10.90);
        // push(2500, 5.636, 6.560e-01, 7267, 10.89);
        // push(2525, 5.821, 6.609e-01, 7410, 10.88);
        // push(2550, 6.018, 6.657e-01, 7558, 10.86);
        // push(2575, 6.229, 6.705e-01, 7711, 10.85);
        // push(2600, 6.455, 6.754e-01, 7870, 10.83);
        // push(2625, 6.695, 6.803e-01, 8034, 10.81);
        // push(2650, 6.951, 6.852e-01, 8204, 10.80);
        // push(2675, 7.224, 6.901e-01, 8382, 10.77);
        // push(2700, 7.513, 6.951e-01, 8566, 10.75);
        // push(2725, 7.820, 7.000e-01, 8757, 10.73);
        // push(2750, 8.145, 7.050e-01, 8957, 10.70);
        // push(2775, 8.489, 7.100e-01, 9165, 10.67);
        // push(2800, 8.854, 7.151e-01, 9382, 10.64);
        // push(2825, 9.238, 7.201e-01, 9608, 10.61);
        // push(2850, 9.644, 7.252e-01, 9844, 10.58);
        // push(2875, 10.07, 7.304e-01, 10090, 10.54);
        // push(2900, 10.52, 7.356e-01, 10350, 10.50);
        // push(2925, 11.00, 7.408e-01, 10620, 10.46);
        // push(2950, 11.50, 7.461e-01, 10900, 10.42);
        // push(2975, 12.02, 7.514e-01, 11190, 10.37);
        // push(3000, 12.57, 7.567e-01, 11500, 10.33);
        // push(3025, 13.15, 7.621e-01, 11820, 10.27);
        // push(3050, 13.76, 7.676e-01, 12160, 10.22);
        // push(3075, 14.39, 7.731e-01, 12510, 10.16);
        // push(3100, 15.05, 7.787e-01, 12880, 10.11);
        // push(3125, 15.73, 7.843e-01, 13260, 10.04);
        // push(3150, 16.45, 7.900e-01, 13660, 9.978);
        // push(3175, 17.19, 7.957e-01, 14080, 9.910);
        // push(3200, 17.96, 8.015e-01, 14520, 9.839);
        // push(3225, 18.75, 8.074e-01, 14980, 9.766);
        // push(3250, 19.56, 8.133e-01, 15460, 9.689);
        // push(3275, 20.39, 8.193e-01, 15960, 9.610);
        // push(3300, 21.23, 8.253e-01, 16480, 9.529);
        // push(3325, 22.09, 8.313e-01, 17020, 9.444);
        // push(3350, 22.94, 8.374e-01, 17580, 9.357);
        // push(3375, 23.80, 8.436e-01, 18170, 9.268);
        // push(3400, 24.65, 8.497e-01, 18770, 9.177);
        // push(3425, 25.49, 8.558e-01, 19400, 9.085);
        // push(3450, 26.31, 8.620e-01, 20050, 8.990);
        // push(3475, 27.09, 8.680e-01, 20720, 8.894);
        // push(3500, 27.85, 8.741e-01, 21400, 8.797);
        // push(3525, 28.56, 8.801e-01, 22110, 8.699);
        // push(3550, 29.22, 8.859e-01, 22830, 8.601);
        // push(3575, 29.83, 8.917e-01, 23570, 8.503);
        // push(3600, 30.39, 8.974e-01, 24320, 8.405);
        // push(3625, 30.88, 9.029e-01, 25090, 8.307);
        // push(3650, 31.31, 9.082e-01, 25870, 8.210);
        // push(3675, 31.67, 9.133e-01, 26650, 8.113);
        // push(3700, 31.96, 9.183e-01, 27450, 8.018);
        // push(3725, 32.19, 9.230e-01, 28250, 7.925);
        // push(3750, 32.35, 9.276e-01, 29060, 7.833);
        // push(3775, 32.44, 9.319e-01, 29870, 7.743);
        // push(3800, 32.47, 9.360e-01, 30680, 7.655);
        // push(3825, 32.44, 9.398e-01, 31490, 7.569);
        // push(3850, 32.35, 9.435e-01, 32300, 7.485);
        // push(3875, 32.20, 9.469e-01, 33110, 7.404);
        // push(3900, 32.00, 9.502e-01, 33910, 7.325);
        // push(3925, 31.75, 9.532e-01, 34710, 7.248);
        // push(3950, 31.46, 9.561e-01, 35500, 7.174);
        // push(3975, 31.13, 9.587e-01, 36280, 7.103);
        //////////////////////////////////////////
        // CHEMICAL EQUILIBRIUM PYROLYSIS GAS

        // push(200, 1.512, 8.688e-02, -7247, 22.00);
        // push(300, 1.631, 1.257e-01, -7090, 22.00);
        // push(400, 1.847, 1.623e-01, -6917, 21.99);
        // push(500, 2.241, 1.962e-01, -6715, 21.95);
        // push(600, 3.210, 2.279e-01, -6451, 21.71);
        // push(700, 6.351, 2.586e-01, -6005, 20.89);
        // push(800, 14.03, 2.909e-01, -5014, 18.64);
        // push(900, 17.01, 3.235e-01, -3335, 15.46);
        // push(1000, 6.222, 3.506e-01, -2170, 13.95);
        // push(1100, 3.509, 3.740e-01, -1741, 13.76);
        // push(1200, 6.235, 3.979e-01, -1316, 13.64);
        // push(1300, 8.634, 4.292e-01, -525.5, 13.11);
        // push(1400, 9.024, 4.589e-01, 353.3, 12.58);
        // push(1500, 10.82, 4.799e-01, 1329, 12.10);
        // push(1600, 12.38, 4.899e-01, 2527, 11.61);
        // push(1700, 7.338, 5.021e-01, 3529, 11.26);
        // push(1800, 4.974, 5.209e-01, 4119, 11.12);
        // push(1900, 4.281, 5.404e-01, 4576, 11.06);
        // push(2000, 4.078, 5.599e-01, 4991, 11.02);
        // FROZEN PYROLYSIS GAS

        push(200, 1.873, 7.392e-02, -4930, 17.89);
        push(300, 2.111, 1.086e-01, -4730, 17.89);
        push(400, 2.327, 1.423e-01, -4508, 17.89);
        push(500, 2.516, 1.745e-01, -4266, 17.89);
        push(600, 2.679, 2.052e-01, -4006, 17.89);
        push(700, 2.821, 2.345e-01, -3731, 17.89);
        push(800, 2.947, 2.624e-01, -3442, 17.89);
        push(900, 3.060, 2.892e-01, -3142, 17.89);
        push(1000, 3.164, 3.150e-01, -2830, 17.89);
        push(1100, 3.259, 3.398e-01, -2509, 17.89);
        push(1200, 3.346, 3.638e-01, -2179, 17.89);
        push(1300, 3.426, 3.870e-01, -1840, 17.89);
        push(1400, 3.499, 4.096e-01, -1494, 17.89);
        push(1500, 3.566, 4.314e-01, -1141, 17.89);
        push(1600, 3.626, 4.527e-01, -781.2, 17.89);
        push(1700, 3.682, 4.734e-01, -415.7, 17.89);
        push(1800, 3.733, 4.936e-01, -44.98, 17.89);
        push(1900, 3.779, 5.134e-01, 142.2, 17.89);
        push(2000, 3.822, 5.327e-01, 710.7, 17.89);

        
}

// Char solid enthalpy tablosu (thermalProperties.csv, SI)
static const vector<double> hc_T_arr   = {256, 298, 444, 556, 644, 833,
                                           1111, 1389, 1667, 1944, 2222, 2778, 3333};
static const vector<double> hc_val_arr = {-32160, 0, 137300, 271300, 393600, 687000,
                                           1175000, 1705000, 2263000, 28439000, 3422000,
                                           4602000, 5793000};

double h_c_T(double T) { return interp1_linear(hc_T_arr, hc_val_arr, T); }  // J/kg
double h_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.h, T); }   // J/kg

// Virgin solid enthalpy tablosu (thermalProperties.csv, SI)
static const vector<double> hv_T_arr   = {256,   298,   444,   556,   644,   833,
                                          1111,  1389,  1667,  1944,  2222,  2778,  3333};
static const vector<double> hv_val_arr = {-896700, -857100, -690100, -536500, -401600, -91240,
                                           405900,  933400, 1477000, 2028000, 2583000,
                                          3697000, 4813000};

double h_v_T(double T) { return interp1_linear(hv_T_arr, hv_val_arr, T); }  // J/kg

double cp_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.cp, T); }
double mu_g_T(double T) { return 1.845e-5; }  // sabit hava viskozitesi
double MW_g_T(double T)
{
    return interp1_linear(tbl_pg.T, tbl_pg.MW, T);
}

// Qp=0: TACOT tablosunda heat of pyrolysis ayri tanimlanmiyor (Ewing Eq.26)
// double Q_p_T(double /*T*/) { return 0.0; }


// Qp: heat of pyrolysis [J/kg] — TACOT tablosundan, endotermik → negatif



// =============================================================================
// BPRIME TABLOSU
// =============================================================================

struct BprimeTable
{
    vector<double> bg;
    vector<vector<double>> Bc;
    vector<vector<double>> Hw;    // J/kg
};

BprimeTable bpt;

void load_bprime_table(const string& fname)
{
    ifstream fin(fname);
    if (!fin.is_open()) {
        printf("ERROR: '%s' acilamiyor!\n", fname.c_str());
        exit(1);
    }

    double prev_bg = -1e99;
    int ibg = -1;

    string line;
    while (getline(fin, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        double bg, tw, bc, hw;
        if (sscanf(line.c_str(), "%lf %lf %lf %lf", &bg, &tw, &bc, &hw) != 4) continue;

        if (bg != prev_bg)
        {
            bpt.bg.push_back(bg);
            bpt.Bc.push_back(vector<double>());
            bpt.Hw.push_back(vector<double>());
            ibg++;
            prev_bg = bg;
        }

        bpt.Bc[ibg].push_back(bc);
        bpt.Hw[ibg].push_back(hw);
    }

    const int nBg  = (int)bpt.bg.size();
    const int nRow = (nBg > 0) ? (int)bpt.Bc[0].size() : 0;
    printf("[TABLE] %s: %d Bg grubu, %d satir/grup\n", fname.c_str(), nBg, nRow);
}

static void bg_bracket(double Bg, int& i0, int& i1, double& t)
{
    const int N = (int)bpt.bg.size();
    if (N == 0) { i0 = i1 = 0; t = 0.0; return; }

    if (Bg <= bpt.bg[0])     { i0 = 0;   i1 = 0;   t = 0.0; return; }
    if (Bg >= bpt.bg[N - 1]) { i0 = N-1; i1 = N-1; t = 0.0; return; }

    int i = 0;
    while (i < N - 1 && Bg > bpt.bg[i + 1]) i++;

    i0 = i;
    i1 = i + 1;
    t  = (Bg - bpt.bg[i0]) / (bpt.bg[i1] - bpt.bg[i0]);
}

// Tw grid: 250, 275, 300, ..., 4000  (adim 25 K, 151 satir/grup)
static const double TW_MIN  = 250.0;
static const double TW_STEP = 25.0;

static double tw_row_interp(const vector<double>& row, double Tw)
{
    const int N = (int)row.size();
    if (N == 0) return 0.0;

    double idx = (Tw - TW_MIN) / TW_STEP;
    if (idx < 0.0)           idx = 0.0;
    if (idx > (double)(N-1)) idx = (double)(N-1);

    const int j0 = (int)idx;
    const int j1 = (j0 < N-1) ? j0+1 : j0;
    const double tj = idx - j0;
    return lerp(row[j0], row[j1], tj);
}

double lookup_Bc(double Bg, double Tw)
{
    int i0, i1; double t;
    bg_bracket(Bg, i0, i1, t);
    return lerp(tw_row_interp(bpt.Bc[i0], Tw),
                tw_row_interp(bpt.Bc[i1], Tw), t);
}

double lookup_Hw(double Bg, double Tw)   // J/kg doner
{
    int i0, i1; double t;
    bg_bracket(Bg, i0, i1, t);
    return lerp(tw_row_interp(bpt.Hw[i0], Tw),
                tw_row_interp(bpt.Hw[i1], Tw), t);
}

// =============================================================================
// NEWTON-RAPHSON on Tw  — Ewing Eq. 86 + 87 (unity Le)
//
//   T_chem = (-1 - B') * h_w(Bg,Tw)  +  B'_c(Bg,Tw) * h_c(Tw)  +  B'_g * h_g(Tw)
//   B' = B'_g + B'_c
//
//   f(Tw) = rho_ue_CH * (H_r + T_chem)
//          - eps * sigma * (Tw^4 - T_surr^4)
//          - k_surf/dx * (Tw - T1)
//          = 0
// =============================================================================
double solve_Tw_NR(double Bg,
                   double rho_ue_CH,
                   double H_r,
                   double k_surf,
                   double dx_surf,
                   double T1,
                   double emissivity,
                   double sigma_SB,
                   double T_surr,
                   double Tw_guess)
{
    const double dTw   = 1.0;      // finite-diff adimi [K]
    const double Tw_min = 300.0;
    const double Tw_max = 4000.0;
    double Tw = Tw_guess;

    auto eval_f = [&](double Twx) -> double
    {
        double Bc    = lookup_Bc(Bg, Twx);
        double hw    = lookup_Hw(Bg, Twx);
        double Bp    = Bg + Bc;
        double hc    = h_c_T(Twx);
        double hg    = h_g_T(Twx);
        double Tchem = (-1.0 - Bp)*hw + Bc*hc + Bg*hg;

        return rho_ue_CH * (H_r + Tchem)
             - emissivity * sigma_SB * (pow(Twx, 4.0) - pow(T_surr, 4.0))
             - (k_surf / dx_surf) * (Twx - T1);
    };

    for (int iter = 0; iter < 300; iter++)
    {
        double f    = eval_f(Tw);
        double Tw2  = min(Tw + dTw, Tw_max);
        double dfdT = (eval_f(Tw2) - f) / dTw;

        if (fabs(dfdT) < 1e-30) break;

        double step = -f / dfdT;
        step = max(-200.0, min(200.0, step));
        Tw  += step;
        Tw   = max(Tw_min, min(Tw_max, Tw));

        if (fabs(step) < 0.01) break;
    }
    return Tw;
}

// =============================================================================
// YARDIMCI: phi-agirlikli bulk yogunluk hesabi
// =============================================================================
static double calc_rho_bulk(
    double alpha0, double alpha1,
    double alpha2,
    double Gamma,
    const vector<double>& rho_v_comp,
    const vector<double>& rho_c_comp,
    double phi_v, double phi_c,
    double rho_intr_v, double rho_intr_c)
{
    double rho_intr =
        Gamma * ( (rho_v_comp[0] - (rho_v_comp[0]-rho_c_comp[0])*alpha0)
                 +(rho_v_comp[1] - (rho_v_comp[1]-rho_c_comp[1])*alpha1) )
      + (1.0-Gamma)*( rho_v_comp[2] - (rho_v_comp[2]-rho_c_comp[2])*alpha2 );

    double alpha_intr = (rho_intr_v - rho_intr) / (rho_intr_v - rho_intr_c);

    double phi_i = phi_v*(1.0 - alpha_intr) + phi_c*alpha_intr;
    return (1.0 - phi_i) * rho_intr;
}

// =============================================================================
// MAIN
// =============================================================================
int main()
{
    cout << string(80, '=') << "\n";
    cout << "ABLASYON ISI TRANSFERI v16 — Case 5 (Pyrolysis+Pressure+T+Erosion)\n";
    cout << string(80, '=') << "\n";

    // =========================================================================
    // PARAMETRELER
    // =========================================================================
    const double t_end = 1;
    const double dt    = 0.0001;
    const int    nstep = (int)round(t_end / dt);

    const int    N_nodes  = 300;
    const double L_domain = 0.1;

    vector<double> x(N_nodes);
    for (int m = 0; m < N_nodes; m++)
        x[m] = m * L_domain / (N_nodes - 1);

    const int N_comp = 3;
    vector<double> B_arr   = {1.200e4,   4.480e9,   0.0};
    vector<double> Psi_arr = {3.0,       3.0,       0.0};
    vector<double> E_arr   = {71.14e6,   169.98e6,  0.0};

    const double Gamma  = 0.5;
    const double R_univ = 8314.0; //J/kmol/K
    // const double MW_air = 12; //kg/kmol 
    // const double R_air  = R_univ / MW_air; // J/kg/K

    const double T_reac[3] = {333.3, 555.6, 5556.0};

    vector<double> E_over_R = {8.556e3, 2.044e4, 0.0};

    vector<double> rho_v_comp = {300.0,  900.0,  1600.0};
    vector<double> rho_c_comp = {  0.0,  600.0,  1600.0};

    vector<double> B_eff_arr = B_arr;
    for (int c = 0; c < N_comp; c++)
    {
        if (B_arr[c] == 0.0) { B_eff_arr[c] = 0.0; continue; }
        const double rv = rho_v_comp[c];
        const double rc = rho_c_comp[c];
        const double dr = rv - rc;
        if (fabs(rv) < 1e-14 || fabs(dr) < 1e-14) { B_eff_arr[c] = 0.0; continue; }
        B_eff_arr[c] = B_arr[c] * pow(dr / rv, Psi_arr[c] - 1.0);
    }

    const double phi_v      = 0.1,   phi_c      = 0.1;
    const double K_v        = 1e-13, K_c = 1e-13;  // permeability [m^2]

    const double rho_intr_v = Gamma*(rho_v_comp[0]+rho_v_comp[1])
                            + (1.0-Gamma)*rho_v_comp[2];
    const double rho_intr_c = Gamma*(rho_c_comp[0]+rho_c_comp[1])
                            + (1.0-Gamma)*rho_c_comp[2];

    const double rho_virgin     = (1.0-phi_v) * rho_intr_v;
    const double rho_char_total = (1.0-phi_c) * rho_intr_c;

    const double T_back  = 300.0, P_surf = 200000, P_back = 100000;

    // =========================================================================
    // Case 5 SINIR KOSULU PARAMETRELERİ (Ewing 2013, Sec. IV.E)
    // =========================================================================
    const double rho_ue_CH0 = 0.3;        // enthalpy-based HTC [kg/m2/s], blowing yokken
    // const double H_recovery = 2.5e7;      // recovery enthalpy [J/kg] — Case 2.3
    const double H_recovery=1.5e6;
    const double t_ramp_end = 0.1;        // HTC ramp süresi [s]: 0 -> rho_ue_CH0 lineer

    const double sigma_SB =5.670374419e-8;
    const double T_surr   = 300.0;

    const int    max_iter = 400;
    const double eps_T    = 0.1;
    const double eps_sdot = 1e-9;

    printf("\n[INIT] dt=%.4f s, t_end=%.1f s, adim=%d\n", dt, t_end, nstep);
    printf("  rho_intr_v=%.1f  rho_intr_c=%.1f kg/m3\n", rho_intr_v, rho_intr_c);
    printf("  rho_virgin=%.1f  rho_char=%.1f kg/m3  (phi dahil)\n",
           rho_virgin, rho_char_total);
    printf("  dx=%.3f mm, N=%d\n", (x[1]-x[0])*1e3, N_nodes);
    printf("  [Case5] rho_ue_CH0=%.3f kg/m2/s  H_r=%.3e J/kg  ramp=%.2fs\n",
           rho_ue_CH0, H_recovery, t_ramp_end);

    init_virgin_table();
    init_char_table();
    init_pyrogas_table();
    load_bprime_table("bprime_table.txt");

    // =========================================================================
    // BASLANGIC KOSULLARI
    // =========================================================================
    vector<double>         T_old(N_nodes, T_back);
    vector<double>         T_prev(N_nodes, T_back);
    vector<double>         P_old(N_nodes, P_back);
    vector<double>         rho_solid_old(N_nodes, rho_virgin);
    vector<vector<double>> alpha_old(N_comp, vector<double>(N_nodes, 0.0));

    double recession_total = 0.0;

    double T_wall        = T_old[0];
    double m_dot_g       = 0.0;
    double k_surf;
    double h_eff         = rho_ue_CH0;   // blowing-corrected rho_ue_CH [kg/m2/s]
    double m_dot_surface = 0.0;

    // =========================================================================
    // CIKTI DOSYALARI
    // =========================================================================
    ofstream fout("ablation_history.csv");
    fout << "time,Twall,T1,P0,mdot,rho_ue_CH,sdot_mm_s,recession_mm,thickness_mm,"
            "Bg,Bc,Tw_nr,eps_surf,k_surf,iter,converged,resid_pct,mdot_g,mdot_c,T_24mm\n";

    const int save_every = 2;

    printf("\n[LOOP] Basliyor...\n\n");

    // =========================================================================
    // LOOP DISI VEKTORLER
    // =========================================================================
    vector<vector<double>> alpha_new(N_comp, vector<double>(N_nodes, 0.0));
    vector<double> rho_solid_new(N_nodes);
    vector<double> alpha_eff(N_nodes);
    vector<double> alpha_eff_old(N_nodes, 0.0);
    vector<double> drho_dt(N_nodes);
    vector<double> k_node(N_nodes);
    vector<double> dT_dt(N_nodes);
    vector<double> a_P(N_nodes), b_P(N_nodes), c_P(N_nodes), d_P(N_nodes);
    vector<double> P_new(N_nodes);
    vector<double> a_T(N_nodes), b_T(N_nodes), c_T(N_nodes), d_T(N_nodes);
    vector<double> T_new(N_nodes);

    double sdot   = 0.0;
    double Bg_now = 0.0, Bc_now = 0.0;

    int n_converged     = 0;
    int n_not_converged = 0;
    int max_iters_used  = 0;

    // =========================================================================
    // TIME LOOP
    // =========================================================================
    for (int n = 1; n <= nstep; n++)
    {

        
        double time    = n * dt;
        double dx_surf = x[1] - x[0];

        // Case 5: lineer ramp — t <= 0.1s arası 0'dan rho_ue_CH0'a çıkar
        double ramp = (time < t_ramp_end) ? (time / t_ramp_end) : 1.0;
        double rho_ue_CH_now = ramp * rho_ue_CH0;
        // =====================================================================
        // m_dot_g (pyroliz gaz uretimi entegrali) ve m_dot_c
        // mdot_c (bir adim onceki Bc ile, timestep basinda)
        double mdot_c_out = Bc_now * h_eff;   // [kg/m2/s]
        alpha_eff_old = alpha_eff;

        // =====================================================================
        // 1. PIROLIZ
        // =====================================================================
        for (int c = 0; c < N_comp; c++)
        {
            if (B_eff_arr[c] == 0.0) { alpha_new[c] = alpha_old[c]; continue; }
            for (int i = 0; i < N_nodes; i++)
            {
                if (T_old[i] <= T_reac[c]) {
                    alpha_new[c][i] = alpha_old[c][i];
                    continue;
                }
                pyr_I  = B_eff_arr[c] * exp(-E_over_R[c] / T_old[i]) * dt;
                pyr_ai = alpha_old[c][i];
                if (Psi_arr[c] == 1.0)
                {
                    pyr_anew = 1.0 - (1.0 - pyr_ai) * exp(-pyr_I);
                }
                else
                {
                    pyr_exp  = 1.0 / (1.0 - Psi_arr[c]);
                    pyr_bkt  = pyr_I*(Psi_arr[c]-1.0) + pow(1.0-pyr_ai, 1.0-Psi_arr[c]);
                    pyr_anew = 1.0 - pow(pyr_bkt, pyr_exp);
                }
                alpha_new[c][i] = pyr_anew;
            }
        }

        // =====================================================================
        // 2. BULK YOGUNLUK, alpha_eff, k_node
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
        {
            rho_solid_new[i] = calc_rho_bulk(
                alpha_new[0][i], alpha_new[1][i], alpha_new[2][i],
                Gamma, rho_v_comp, rho_c_comp,
                phi_v, phi_c, rho_intr_v, rho_intr_c);

            alpha_eff[i] = (rho_virgin - rho_solid_new[i])
                         / (rho_virgin - rho_char_total);

            double phi_i = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i];
            drho_dt[i] = 0.0;
            double weight[3] = {Gamma, Gamma, 1.0-Gamma};
            for (int c = 0; c < N_comp; c++)
                drho_dt[i] += weight[c] * (rho_v_comp[c]-rho_c_comp[c])
                            * (alpha_new[c][i]-alpha_old[c][i]) / dt;
            drho_dt[i] *= (1.0 - phi_i);
            k_node[i] = k_mix(T_old[i], alpha_eff[i]);
        }

        // mdot_g_total: bu adimin drho_dt'sinden integral
        // rho azaliyor → drho_dt > 0 (alpha artiyor, yogunluk azaliyor)
        // gaz uretimi = -d(rho_solid)/dt * dx → pozitif olmali
        double mdot_g_total = 0.0;
        for (int i = 0; i < N_nodes; i++)
        {
            double dx_i;
            if (i == 0)
                dx_i = 0.5*(x[1] - x[0]);
            else if (i == N_nodes-1)
                dx_i = 0.5*(x[N_nodes-1] - x[N_nodes-2]);
            else
                dx_i = 0.5*(x[i+1] - x[i-1]);
            mdot_g_total += drho_dt[i] * dx_i;   // drho_dt pozitif → rho azaliyor
        }

        // =====================================================================
        // 3. dT/dt
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
            dT_dt[i] = (T_old[i] - T_prev[i]) / dt;

        // =====================================================================
        // 4. BASINC
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
        {
            a_P[i] = 0.0; b_P[i] = 0.0; c_P[i] = 0.0; d_P[i] = 0.0;
        }
        a_P[0] = 1.0;         d_P[0]         = P_surf;
        a_P[N_nodes-1] = 1.0; d_P[N_nodes-1] = P_back;

        for (int i = 1; i < N_nodes-1; i++)
        {
            p_dxl = x[i]   - x[i-1];
            p_dxr = x[i+1] - x[i];
            p_dxi = 0.5*(p_dxl + p_dxr);

            p_phi   = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i];
            double R_air=R_univ /MW_g_T(T_old[i]);
            p_atime = p_phi / (R_air*T_old[i]) * p_dxi / dt;

            double K_i  = K_v*(1.0-alpha_eff[i])   + K_c*alpha_eff[i];
            double K_im = K_v*(1.0-alpha_eff[i-1]) + K_c*alpha_eff[i-1];
            double K_ip = K_v*(1.0-alpha_eff[i+1]) + K_c*alpha_eff[i+1];

            double P_e = 0.5*(P_old[i] + P_old[i+1]);
            double P_w = 0.5*(P_old[i] + P_old[i-1]);

            double K_e = 2.0*K_i*K_ip / (K_i + K_ip);
            double K_w = 2.0*K_i*K_im / (K_i + K_im);
            p_ae = P_e / (R_air*T_old[i]) * K_e / mu_g_T(T_old[i]) / p_dxr;
            p_aw = P_w / (R_air*T_old[i]) * K_w / mu_g_T(T_old[i]) / p_dxl;

            p_aeff_old = (rho_virgin - rho_solid_old[i]) / (rho_virgin - rho_char_total);
            p_dalpha   = (alpha_eff[i] - p_aeff_old) / dt;
            p_Sptemp   = p_phi / (R_air*T_old[i]*T_old[i]) * dT_dt[i];
            p_Spporo   = (phi_v-phi_c) / (R_air*T_old[i]) * p_dalpha;
            p_Slin     = p_Sptemp + p_Spporo;
            double d_alpha_eff_dt = (alpha_eff[i] - alpha_eff_old[i]) / dt;
            a_P[i] = p_atime + p_ae + p_aw - p_Slin*p_dxi;
            b_P[i] = p_ae;
            c_P[i] = p_aw;
            d_P[i] = (rho_virgin-rho_char_total) * d_alpha_eff_dt * p_dxi + p_atime*P_old[i];
        }

        P_new = thomas_patankar(a_P, b_P, c_P, d_P);

        // =====================================================================
        // 5. BLOWING (ilk tahmin)
        // =====================================================================
        {
            double K_0    = K_v*(1.0-alpha_eff[0]) + K_c*alpha_eff[0];
            double K_1    = K_v*(1.0-alpha_eff[1]) + K_c*alpha_eff[1];
            double K_surf_tmp = 2.0*K_0*K_1 / (K_0 + K_1);
            double R_air=R_univ /MW_g_T(T_old[0]);
            m_dot_surface = (P_new[0]/(R_air*T_old[0]))
                          * (K_surf_tmp/mu_g_T(T_old[0]))
                          * (P_new[1]-P_new[0]) / dx_surf;
        }
        // m_dot_surface*=10.0;  // Case 5'te ilk iterasyonda HTC çok düşük geliyor, ramp etkisiyle uyumlu olması için 10x artırarak başlatıyorum  
        // if (rho_ue_CH_now > 0.0)
        //     blowing_factor(m_dot_surface, rho_ue_CH_now, h_eff);
        // else
        //     h_eff = 0.0;
        m_dot_g = m_dot_surface;

        // =====================================================================
        // 6. ITERATIF YUZEY + SICAKLIK
        // =====================================================================
        k_surf = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);

        double T1        = T_old[1];
        double sdot_iter = sdot;
        int    iter_count = 0;

        for (int it = 0; it < max_iter; it++)
        {
            iter_count = it + 1;

            // (i) mdot / h_eff guncelle
            {
                double K_0    = K_v*(1.0-alpha_eff[0]) + K_c*alpha_eff[0];
                double K_1    = K_v*(1.0-alpha_eff[1]) + K_c*alpha_eff[1];
                double K_surf2 = 2.0*K_0*K_1 / (K_0 + K_1);
                double R_air=R_univ /MW_g_T(T_wall);
                m_dot_surface = (P_new[0]/(R_air*T_wall))
                              * (K_surf2/mu_g_T(T_wall))
                              * (P_new[1]-P_new[0]) / dx_surf;
            }
            // if (rho_ue_CH_now > 0.0)
            //     blowing_factor(m_dot_surface, rho_ue_CH_now, h_eff);
            // else
            //     h_eff = 0.0;
            m_dot_g = m_dot_surface;

            double T_wall_new = T_wall;

            // (a) NR on Tw — Bg = mdot_g_total / rho_ue_CH_now
            Bg_now = (rho_ue_CH_now > 0.0) ? (mdot_g_total / h_eff) : 0.0;
            if (Bg_now < 0.0) Bg_now = 0.0;

            // blowing correction: toplam yuzey kitlesi = mdot_g + mdot_c
            double mdot_total_blow = mdot_g_total + Bc_now * rho_ue_CH_now;
            if (mdot_total_blow < 0.0) mdot_total_blow = 0.0;
            if (rho_ue_CH_now > 0.0)
                blowing_factor(mdot_total_blow, rho_ue_CH_now, h_eff);
            else
                h_eff = 0.0;

            double emissivity_now = eps_surf(T_wall, alpha_eff[0]);
            T_wall_new = solve_Tw_NR(Bg_now, h_eff, H_recovery,
                                     k_surf, dx_surf, T1,
                                     emissivity_now, sigma_SB, T_surr,
                                     T_wall);

            Bc_now = lookup_Bc(Bg_now, T_wall_new);

            // sdot = B'c * rho_ue_CH / rho_c  (enthalpy-based, Ewing Eq.52)
            sdot = (rho_ue_CH_now > 0.0) ? (Bc_now * h_eff / rho_char_total) : 0.0;
            if (sdot < 0.0) sdot = 0.0;

            double sdot_old_iter = sdot_iter;
            sdot_iter = sdot;

            // (ii) m_dot_surface güncelle (T_wall_new ile)
            {
                double K_0    = K_v*(1.0-alpha_eff[0]) + K_c*alpha_eff[0];
                double K_1    = K_v*(1.0-alpha_eff[1]) + K_c*alpha_eff[1];
                double K_surf2 = 2.0*K_0*K_1 / (K_0 + K_1);
                double R_air=R_univ /MW_g_T(T_wall_new);
                m_dot_surface = (P_new[0]/(R_air*T_wall_new))
                              * (K_surf2/mu_g_T(T_wall_new))
                              * (P_new[1]-P_new[0]) / dx_surf;
            }

            // (b) SICAKLIK TDMA
            for (int i = 0; i < N_nodes; i++)
                a_T[i] = b_T[i] = c_T[i] = d_T[i] = 0.0;

            a_T[0] = 1.0; d_T[0] = T_wall_new;

            for (int i = 1; i < N_nodes-1; i++)
            {
                t_dxl = x[i]   - x[i-1];
                t_dxr = x[i+1] - x[i];
                t_dxi = 0.5*(t_dxl + t_dxr);

                t_cp  = cp_mix(T_old[i], alpha_eff[i]);
                t_phi = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i];

                t_ke = 2.0*k_node[i]*k_node[i+1] / (k_node[i]+k_node[i+1]);
                t_kw = 2.0*k_node[i]*k_node[i-1] / (k_node[i]+k_node[i-1]);
                double R_air=R_univ /MW_g_T(T_old[i]);
                t_rhogas = P_new[i] / (R_air*T_old[i]);
                t_rhoc   = rho_solid_new[i]*t_cp
                         + t_rhogas*t_phi*cp_g_T(T_old[i]);
                t_atime  = t_rhoc * t_dxi / dt;
                t_ae     = t_ke / t_dxr;
                t_aw     = t_kw / t_dxl;

                double K_i2 = K_v*(1.0-alpha_eff[i]) + K_c*alpha_eff[i];
                t_mdotgas = t_rhogas * (K_i2/mu_g_T(T_old[i]))
                          * (P_new[i+1]-P_new[i-1]) / (t_dxl+t_dxr);
                t_hgrad   = cp_g_T(T_old[i])
                          * (T_old[i+1]-T_old[i-1]) / (t_dxl+t_dxr);

                double d_alpha_eff_dt = (alpha_eff[i] - alpha_eff_old[i]) / dt;

                // Ewing Eq.30: hs = hv*(1-alpha) + hc*alpha
                // Ewing Eq.34: h* = hs + rho_s*(hv-hc)/(rho_v-rho_c)
                // Ewing Eq.35 (Qp=0): S_pyro = (h* - hg)*(rho_v-rho_c)*dalpha/dt
                double hv_i    = h_v_T(T_old[i]);
                double hc_i    = h_c_T(T_old[i]);
                double hg_i    = h_g_T(T_old[i]);
                double hs_i    = hv_i*(1.0 - alpha_eff[i]) + hc_i*alpha_eff[i];
                // double hstar_i = hs_i + rho_solid_new[i]*(hv_i - hc_i)
                //                        / (rho_virgin - rho_char_total);
                // // Qp=0 (TACOT) → only h* and hg terms remain
                // t_Spyro  = (hstar_i - hg_i)
                //            * (rho_virgin - rho_char_total) * d_alpha_eff_dt;

                double hbar_i = (rho_virgin*hv_i - rho_char_total*hc_i)
                            / (rho_virgin - rho_char_total);
                // Volkan Eq.2.5-2.6: S_pyr = -drho_dt * (hg - hbar)
                t_Spyro = -drho_dt[i] * (hg_i - hbar_i);

                t_Stotal = t_Spyro + t_mdotgas*t_hgrad;

                a_T[i] = t_atime + t_ae + t_aw;
                b_T[i] = t_ae;
                c_T[i] = t_aw;
                d_T[i] = t_Stotal*t_dxi + t_atime*T_old[i];
            }

            // Arka yuzey: adiabatik
            {
                int jj   = N_nodes - 1;
                t_dxl    = x[jj] - x[jj-1];
                t_cp     = cp_mix(T_old[jj], alpha_eff[jj]);
                t_phi    = phi_v*(1.0-alpha_eff[jj]) + phi_c*alpha_eff[jj];
                double R_air=R_univ /MW_g_T(T_old[jj]);
                t_rhogas = P_new[jj] / (R_air*T_old[jj]);
                t_rhoc   = rho_solid_new[jj]*t_cp
                         + t_rhogas*t_phi*cp_g_T(T_old[jj]);
                t_atime  = t_rhoc * (0.5*t_dxl) / dt;
                t_kw     = 2.0*k_node[jj]*k_node[jj-1] / (k_node[jj]+k_node[jj-1]);
                a_T[jj]  = t_atime + t_kw/t_dxl;
                b_T[jj]  = 0.0;
                c_T[jj]  = t_kw/t_dxl;
                d_T[jj]  = t_atime * T_old[jj];
            }

            T_new = thomas_patankar(a_T, b_T, c_T, d_T);

            double T1_new = T_new[1];
            double dT1    = fabs(T1_new   - T1);
            double dsdot  = fabs(sdot_iter - sdot_old_iter);

            T_wall = T_wall_new;
            T1     = T1_new;

            if (dT1 < eps_T && dsdot < eps_sdot) break;
        }

        bool iter_converged = (iter_count < max_iter);
        if (iter_converged) n_converged++;
        else                n_not_converged++;
        if (iter_count > max_iters_used) max_iters_used = iter_count;

        // Yuzey enerji dengesi residual (enthalpy-based)
        // Yüzey enerji dengesi residual (Eq. 86+87)
        double Tw_f   = T_wall;
        double Bc_f   = lookup_Bc(Bg_now, Tw_f);
        double hw_f   = lookup_Hw(Bg_now, Tw_f);
        double Bp_f   = Bg_now + Bc_f;
        double hc_f   = h_c_T(Tw_f);
        double hg_f   = h_g_T(Tw_f);
        double Tchem_f = (-1.0 - Bp_f)*hw_f + Bc_f*hc_f + Bg_now*hg_f;

        double q_conv_s  = h_eff * (H_recovery + Tchem_f);
        double q_rad_s   = 1 * sigma_SB * pow(Tw_f, 4.0);
        double q_cond_s  = k_surf * (Tw_f - T_new[1]) / dx_surf;
        double resid_s   = q_conv_s - q_rad_s - q_cond_s;
        double resid_pct_s = (fabs(q_conv_s) > 1e-10) ? (resid_s / q_conv_s * 100.0) : 0.0;

        // =====================================================================
        // 7. MOVING BOUNDARY + REMAP
        // =====================================================================
        double recession_step = 0;
        recession_total += recession_step;

        vector<double> x_old = x;

        if (recession_step > 0.0)
        {
            if ((L_domain - recession_total) < 2.0*dx_surf)
            {
                printf("\nWARNING: Malzeme tamamen eridi!\n");
                break;
            }

            for (int jj = 0; jj < N_nodes; jj++)
                x[jj] += recession_step * (N_nodes-1-jj) / (N_nodes-1);

            for (int jj = 0; jj < N_nodes-1; jj++)
            {
                double d1 = x[jj]       - x_old[jj];
                double d2 = x_old[jj+1] - x[jj];

                T_new[jj]         = inverseAverage(T_new[jj],         d1, T_new[jj+1],         d2);
                P_new[jj]         = inverseAverage(P_new[jj],         d1, P_new[jj+1],         d2);
                rho_solid_new[jj] = inverseAverage(rho_solid_new[jj], d1, rho_solid_new[jj+1], d2);

                for (int c = 0; c < N_comp; c++)
                    alpha_new[c][jj] = inverseAverage(alpha_new[c][jj], d1, alpha_new[c][jj+1], d2);
            }

            T_new[0] = T_wall;
            P_new[0] = P_surf;

            for (int jj = 0; jj < N_nodes; jj++)
            {
                rho_solid_new[jj] = calc_rho_bulk(
                    alpha_new[0][jj], alpha_new[1][jj], alpha_new[2][jj],
                    Gamma, rho_v_comp, rho_c_comp,
                    phi_v, phi_c, rho_intr_v, rho_intr_c);

                alpha_eff[jj] = (rho_virgin - rho_solid_new[jj])
                              / (rho_virgin - rho_char_total);

                k_node[jj] = k_mix(T_old[jj], alpha_eff[jj]);
            }
            k_surf = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);
        }

        if (n*dt==0.1)
        {
            for (int i = 0; i < N_nodes; i++)
            {
                cout<<x[i] << "," << P_old[i]  << "\n";    
            }
            
        }

        // 24mm'deki sicaklik — malzeme koordinatinda sabit nokta
        // x[0] = yuzey (kayiyor), 24mm orijinal derinlik
        // En yakin nodu bul
        double T_24mm = T_new[0];
        {
            double target = 0.024;
            double best = 1e99;
            for (int i = 0; i < N_nodes; i++)
            {
                double dist = fabs(x[i] - target);
                if (dist < best) { best = dist; T_24mm = T_new[i]; }
            }
        }

        // =====================================================================
        // 8. UPDATE
        // =====================================================================
        T_prev        = T_old;
        T_old         = T_new;
        P_old         = P_new;
        alpha_old     = alpha_new;
        rho_solid_old = rho_solid_new;

        // =====================================================================
        // OUTPUT
        // =====================================================================
        double thickness_now = L_domain - recession_total;

        if ((n % (save_every*10) == 0 || n == 1 || n == nstep))
        {
            fout << time << "," << T_wall << "," << T_old[1] << ","
                 << P_new[0] << "," << m_dot_surface << "," << h_eff << ","
                 << sdot*1e3 << "," << recession_total*1e3 << ","
                 << thickness_now*1e3 << ","
                 << Bg_now << "," << Bc_now << "," << T_wall << ","
                 << eps_surf(T_wall, alpha_eff[0]) << "," << k_surf << ","
                 << iter_count << ","
                 << (iter_converged ? 1 : 0) << ","
                 << resid_pct_s << ","
                 << mdot_g_total << "," << mdot_c_out << "," << T_24mm << "\n";
        }
        // if ((n % (save_every*10)) == 0 || n == 1 || n == nstep)
        // {
        //     printf("t=%7.1fs | Tw=%7.1fK | T1=%7.1fK | rho_ue_CH=%6.4f | "
        //            "Bg=%.3f Bc=%.4f L=%5.2f eps=%.3f k=%.3f | "
        //            "sdot=%7.4fmm/s | erim=%7.4fmm | "
        //            "it=%d %s | resid=%+.3f%%\n",
        //            time, T_wall, T_old[1], h_eff,
        //            Bg_now, Bc_now, L_prev,
        //            eps_surf(T_wall, alpha_eff[0]), k_surf,
        //            sdot*1e3, recession_total*1e3,
        //            iter_count,
        //            iter_converged ? "[CONV]" : "[WARN:MAX_ITER]",
        //            resid_pct_s);
        // }
    }

    fout.close();

    printf("\n%s\n", string(80, '=').c_str());
    printf("TAMAMLANDI\n");
    printf("Final Twall         : %.2f K\n", T_wall);
    printf("Final T[1]          : %.2f K\n", T_old[1]);
    printf("Final P[0]          : %.3f kPa\n", P_old[0]/1000.0);
    printf("Final mdot          : %.3e kg/m2s\n", m_dot_surface);
    printf("Final rho_ue_CH_eff : %.4f kg/m2s\n", "dikkat!");
    printf("Final s_dot         : %.4f mm/s\n", sdot*1e3);
    printf("Final Bg            : %.4f\n", Bg_now);
    printf("Final Bc            : %.4f\n", Bc_now);
    printf("Toplam erime        : %.4f mm\n", recession_total*1e3);
    printf("Kalan kalinlik      : %.4f mm\n", (L_domain-recession_total)*1e3);

    double h_w_f  = cp_g_T(T_wall) * T_wall;
    double q_conv = h_eff * (H_recovery - h_w_f);
    double q_rad  = -eps_surf(T_wall, alpha_eff[0]) * sigma_SB * (pow(T_surr,4) - pow(T_wall,4));
    double q_chem = h_eff * 1;
    double q_cond = k_surf * (T_wall - T_old[1]) / (x[1]-x[0]);
    double resid  = q_conv + q_rad + q_chem - q_cond;

    printf("\nYuzey enerji dengesi (kW/m2):\n");
    printf("  q_conv = %+12.3f kW/m2\n", q_conv/1000.0);
    printf("  q_rad  = %+12.3f kW/m2\n", q_rad/1000.0);
    printf("  q_chem = %+12.3f kW/m2\n", q_chem/1000.0);
    printf("  q_cond = %+12.3f kW/m2\n", q_cond/1000.0);
    printf("  resid  = %+12.3f kW/m2\n", resid/1000.0);
    printf("  Convergence Rate = %.6f %%\n", (fabs(q_conv)>1e-10)?(resid/q_conv)*100.0:0.0);

    int n_total_steps = n_converged + n_not_converged;
    printf("\nYAKINSAMA OZETI:\n");
    printf("  Toplam zaman adimi      : %d\n", n_total_steps);
    printf("  Yakinsayan adim         : %d  (%.1f%%)\n",
           n_converged,
           n_total_steps > 0 ? 100.0 * n_converged / n_total_steps : 0.0);
    printf("  Yakinsamayan adim [WARN]: %d  (%.1f%%)\n",
           n_not_converged,
           n_total_steps > 0 ? 100.0 * n_not_converged / n_total_steps : 0.0);
    printf("  Kullanilan maks iter    : %d / %d\n", max_iters_used, max_iter);
    if (n_not_converged == 0)
        printf("  Sonuc: TUM ADIMLAR YAKINSADI.\n");
    else
        printf("  Sonuc: DIKKAT — %d adimda yakinsama saglanamadi!\n", n_not_converged);
    printf("%s\n", string(80, '=').c_str());
    double dx_surf = x[1] - x[0];

    return 0;
}