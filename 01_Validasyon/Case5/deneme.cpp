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
double blowing_factor(double m_dot, double rho_ue_CH0, double& h_eff_out)
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
    push_row_thermal(tbl_v, 256,  879, 0.398, 0.8);
    push_row_thermal(tbl_v, 298,  984, 0.403, 0.8);
    push_row_thermal(tbl_v, 444, 1300, 0.416, 0.8);
    push_row_thermal(tbl_v, 556, 1470, 0.453, 0.8);
    push_row_thermal(tbl_v, 644, 1570, 0.470, 0.8);
    push_row_thermal(tbl_v, 833, 1720, 0.486, 0.8);
    push_row_thermal(tbl_v, 1111,1860, 0.523, 0.8);
    push_row_thermal(tbl_v, 1389,1930, 0.560, 0.8);
    push_row_thermal(tbl_v, 1667,1980, 0.698, 0.8);
    push_row_thermal(tbl_v, 1944,1990, 0.872, 0.8);
    push_row_thermal(tbl_v, 2222,2000, 1.110, 0.8);
    push_row_thermal(tbl_v, 2778,2010, 1.750, 0.8);
    push_row_thermal(tbl_v, 3333,2010, 2.780, 0.8);
}

void init_char_table()
{
    tbl_c = ThermalTable{};
    push_row_thermal(tbl_c, 256,  733, 0.398, 0.9);
    push_row_thermal(tbl_c, 298,  783, 0.403, 0.9);
    push_row_thermal(tbl_c, 444, 1090, 0.416, 0.9);
    push_row_thermal(tbl_c, 556, 1320, 0.453, 0.9);
    push_row_thermal(tbl_c, 644, 1430, 0.470, 0.9);
    push_row_thermal(tbl_c, 833, 1680, 0.486, 0.9);
    push_row_thermal(tbl_c, 1111,1840, 0.523, 0.9);
    push_row_thermal(tbl_c, 1389,1970, 0.560, 0.9);
    push_row_thermal(tbl_c, 1667,2050, 0.605, 0.9);
    push_row_thermal(tbl_c, 1944,2090, 0.729, 0.9);
    push_row_thermal(tbl_c, 2222,2110, 0.922, 0.9);
    push_row_thermal(tbl_c, 2778,2140, 1.460, 0.9);
    push_row_thermal(tbl_c, 3333,2150, 2.320, 0.9);
}

double cp_v_T(double T) { return interp1_linear(tbl_v.T, tbl_v.cp, T); }
double k_v_T (double T) { return interp1_linear(tbl_v.T, tbl_v.k,  T); }
double cp_c_T(double T) { return interp1_linear(tbl_c.T, tbl_c.cp, T); }
double k_c_T (double T) { return interp1_linear(tbl_c.T, tbl_c.k,  T); }

double cp_mix(double T, double a) { return cp_v_T(T)*(1.0-a) + cp_c_T(T)*a; }
double k_mix (double T, double a) { return k_v_T(T) *(1.0-a) + k_c_T(T) *a; }

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
    vector<double> cp;
    vector<double> mu;
    vector<double> h;    // ← YENİ: kJ/kg → J/kg ile kullanılacak
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
    auto push = [&](double T, double cp_kJ, double mu_mp, double h_kJ){
        tbl_pg.T.push_back(T);
        tbl_pg.cp.push_back(cp_kJ * 1000.0);
        tbl_pg.mu.push_back(mu_mp * 1e-4);
        tbl_pg.h.push_back(h_kJ * 1000.0);   // J/kg
    };

    push(200,  1.512, 8.688e-02, -7247.0);
    push(225,  1.534, 9.666e-02, -7208.0);
    push(250,  1.560, 1.065e-01, -7170.0);
    push(275,  1.592, 1.162e-01, -7130.0);
    push(300,  1.631, 1.257e-01, -7090.0);
    push(325,  1.676, 1.351e-01, -7049.0);
    push(350,  1.726, 1.444e-01, -7006.0);
    push(375,  1.783, 1.534e-01, -6963.0);
    push(400,  1.847, 1.623e-01, -6917.0);
    push(425,  1.921, 1.710e-01, -6870.0);
    push(450,  2.008, 1.796e-01, -6821.0);
    push(475,  2.112, 1.879e-01, -6770.0);
    push(500,  2.241, 1.962e-01, -6715.0);
    push(525,  2.402, 2.042e-01, -6657.0);
    push(550,  2.605, 2.122e-01, -6595.0);
    push(575,  2.867, 2.201e-01, -6527.0);
    push(600,  3.210, 2.279e-01, -6451.0);
    push(625,  3.671, 2.356e-01, -6365.0);
    push(650,  4.301, 2.432e-01, -6266.0);
    push(675,  5.171, 2.509e-01, -6148.0);
    push(700,  6.351, 2.586e-01, -6005.0);
    push(725,  7.882, 2.664e-01, -5827.0);
    push(750,  9.748, 2.744e-01, -5608.0);
    push(775, 11.850, 2.826e-01, -5338.0);
    push(800, 14.030, 2.909e-01, -5014.0);
    push(825, 16.010, 2.993e-01, -4638.0);
    push(850, 17.120, 3.077e-01, -4222.0);
    push(875, 16.650, 3.161e-01, -3786.0);
    push(900, 14.380, 3.243e-01, -3360.0);
    push(925, 11.260, 3.322e-01, -2982.0);
    push(950,  8.607, 3.397e-01, -2678.0);
    push(975,  6.823, 3.467e-01, -2444.0);
    push(1000,  5.748, 3.533e-01, -2261.0);
    push(1025,  5.064, 3.596e-01, -2107.0);
    push(1050,  4.620, 3.655e-01, -1972.0);
    push(1075,  4.304, 3.711e-01, -1849.0);
    push(1100,  4.054, 3.763e-01, -1735.0);
    push(1125,  3.849, 3.812e-01, -1628.0);
    push(1150,  3.677, 3.858e-01, -1527.0);
    push(1175,  3.530, 3.900e-01, -1430.0);
    push(1200,  3.401, 3.940e-01, -1337.0);
    push(1225,  3.287, 3.977e-01, -1248.0);
    push(1250,  3.186, 4.011e-01, -1162.0);
    push(1275,  3.096, 4.043e-01, -1079.0);
    push(1300,  3.017, 4.072e-01,  -998.5);
    push(1325,  2.949, 4.099e-01,  -921.0);
    push(1350,  2.893, 4.124e-01,  -845.5);
    push(1375,  2.851, 4.147e-01,  -771.6);
    push(1400,  2.824, 4.167e-01,  -697.7);
    push(1425,  2.813, 4.185e-01,  -623.3);
    push(1450,  2.818, 4.202e-01,  -547.9);
    push(1475,  2.839, 4.218e-01,  -471.0);
    push(1500,  2.875, 4.232e-01,  -392.3);
    push(1525,  2.921, 4.246e-01,  -312.3);
    push(1550,  2.970, 4.260e-01,  -230.8);
    push(1575,  3.016, 4.274e-01,  -148.5);
    push(1600,  3.050, 4.287e-01,   -66.0);
    push(1625,  3.068, 4.300e-01,    17.0);
    push(1650,  3.065, 4.313e-01,    99.7);
    push(1675,  3.042, 4.325e-01,   181.6);
    push(1700,  3.002, 4.337e-01,   262.4);
    push(1725,  2.950, 4.348e-01,   341.8);
    push(1750,  2.892, 4.358e-01,   419.4);
    push(1775,  2.833, 4.368e-01,   495.4);
    push(1800,  2.777, 4.377e-01,   569.8);
    push(1825,  2.727, 4.386e-01,   642.7);
    push(1850,  2.684, 4.394e-01,   714.4);
    push(1875,  2.650, 4.401e-01,   784.8);
    push(1900,  2.626, 4.408e-01,   854.3);
    push(1925,  2.612, 4.415e-01,   923.0);
    push(1950,  2.609, 4.421e-01,   991.3);
    push(1975,  2.617, 4.427e-01,  1059.7);
    push(2000,  2.636, 4.432e-01,  1128.7);
    push(2025,  2.665, 4.438e-01,  1198.5);
    push(2050,  2.703, 4.443e-01,  1269.4);
    push(2075,  2.749, 4.448e-01,  1341.4);
    push(2100,  2.802, 4.453e-01,  1414.7);
    push(2125,  2.861, 4.458e-01,  1489.3);
    push(2150,  2.925, 4.462e-01,  1565.2);
    push(2175,  2.993, 4.466e-01,  1642.5);
    push(2200,  3.065, 4.471e-01,  1721.2);
    push(2225,  3.142, 4.475e-01,  1801.5);
    push(2250,  3.224, 4.479e-01,  1883.4);
    push(2275,  3.312, 4.483e-01,  1967.1);
    push(2300,  3.407, 4.487e-01,  2052.7);
    push(2325,  3.511, 4.492e-01,  2140.4);
    push(2350,  3.624, 4.496e-01,  2230.3);
    push(2375,  3.748, 4.500e-01,  2322.6);
    push(2400,  3.882, 4.504e-01,  2417.5);
    push(2425,  4.029, 4.508e-01,  2514.9);
    push(2450,  4.190, 4.512e-01,  2615.3);
    push(2475,  4.364, 4.516e-01,  2718.7);
    push(2500,  4.552, 4.520e-01,  2825.3);
    push(2550,  4.969, 4.528e-01,  3046.6);
    push(2600,  5.436, 4.537e-01,  3282.4);
    push(2650,  5.948, 4.545e-01,  3533.9);
    push(2700,  6.503, 4.554e-01,  3802.4);
    push(2750,  7.096, 4.562e-01,  4088.7);
    push(2800,  7.720, 4.571e-01,  4393.7);
    push(2850,  8.367, 4.580e-01,  4718.4);
    push(2900,  9.026, 4.589e-01,  5063.5);
    push(2950,  9.685, 4.598e-01,  5429.0);
    push(3000, 10.330, 4.607e-01,  5815.0);
    push(3050, 10.940, 4.617e-01,  6221.0);
    push(3100, 11.500, 4.626e-01,  6647.0);
    push(3200, 12.430, 4.646e-01,  7559.0);
    push(3300, 13.190, 4.666e-01,  8541.0);
    push(3400, 13.760, 4.687e-01,  9580.0);
    push(3500, 14.090, 4.708e-01, 10660.0);
    push(3600, 14.170, 4.730e-01, 11770.0);
    push(3700, 13.990, 4.752e-01, 12890.0);
    push(3800, 13.550, 4.775e-01, 14010.0);
    push(3900, 12.870, 4.798e-01, 15120.0);
    push(3975, 12.230, 4.816e-01, 16010.0);
}

// Char solid enthalpy tablosu (thermalProperties.csv, SI)
static const vector<double> hc_T_arr   = {256, 298, 444, 556, 644, 833,
                                           1111, 1389, 1667, 1944, 2222, 2778, 3333};
static const vector<double> hc_val_arr = {-32200, 0, 137000, 271000, 394000, 687000,
                                           1180000, 1710000, 2260000, 2840000, 3420000,
                                           4600000, 5790000};

double h_c_T(double T) { return interp1_linear(hc_T_arr, hc_val_arr, T); }  // J/kg
double h_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.h, T); }   // J/kg

double cp_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.cp, T); }
double mu_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.mu, T); }

static const vector<double> Qp_T   = {256, 298, 444, 556, 644, 833,
                                       1111, 1389, 1667, 1944, 2222, 2778, 3333};
static const vector<double> Qp_val = {-864540, -857000, -827400, -807800, -795200, -778240,
                                       -769100, -772600, -786000, -808000, -837000, -903000, -977000};

double Q_p_T(double T) { return interp1_linear(Qp_T, Qp_val, T); }

// =============================================================================
// BPRIME TABLOSU
// =============================================================================

struct BprimeTable
{
    vector<double> bg;
    vector<vector<double>> Tw;
    vector<vector<double>> Bc;
    vector<vector<double>> Hw;    // ← YENİ: J/kg (bprime_table.txt 4. kolon)
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
            bpt.Tw.push_back(vector<double>());
            bpt.Bc.push_back(vector<double>());
            bpt.Hw.push_back(vector<double>());
            ibg++;
            prev_bg = bg;
        }

        bpt.Tw[ibg].push_back(tw);
        bpt.Bc[ibg].push_back(bc);
        bpt.Hw[ibg].push_back(hw);   // J/kg olarak geliyor
    }

    const int nBg  = (int)bpt.bg.size();
    const int nRow = (nBg > 0) ? (int)bpt.Tw[0].size() : 0;
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

static double row_interp(const vector<double>& row, double L)
{
    const int N = (int)row.size();
    if (N == 0) return 0.0;
    if (N == 1) return row[0];

    double idx = L - 1.0;
    if (idx < 0.0) idx = 0.0;
    if (idx > (double)(N - 1)) idx = (double)(N - 1);

    const int j0 = (int)idx;
    const int j1 = (j0 < N - 1) ? (j0 + 1) : j0;
    const double tj = idx - j0;

    return lerp(row[j0], row[j1], tj);
}

double lookup_Tw(double Bg, double L)
{
    int i0, i1; double t;
    bg_bracket(Bg, i0, i1, t);
    return lerp(row_interp(bpt.Tw[i0], L), row_interp(bpt.Tw[i1], L), t);
}

double lookup_Bc(double Bg, double L)
{
    int i0, i1; double t;
    bg_bracket(Bg, i0, i1, t);
    return lerp(row_interp(bpt.Bc[i0], L), row_interp(bpt.Bc[i1], L), t);
}

double lookup_Hw(double Bg, double L)    // J/kg döner
{
    int i0, i1; double t;
    bg_bracket(Bg, i0, i1, t);
    return lerp(row_interp(bpt.Hw[i0], L), row_interp(bpt.Hw[i1], L), t);
}

// =============================================================================
// NEWTON-RAPHSON on L  (Ewing Eq.86-90)
// Case 5: enthalpy-based enerji dengesi
//   q_cond = rho_ue_CH*(H_r - h_w) + eps*sigma*(T_surr^4 - Tw^4) + rho_ue_CH*Tchem/cp_g
//   h_w = cp_g * Tw  (simplified, unity Le)
//   f(L) = k_surf/dx*(Tw-T1) - rho_ue_CH*(H_r - cp_g*Tw)
//          - eps*sigma*(T_surr^4-Tw^4) - rho_ue_CH*Tchem
// =============================================================================
// =============================================================================
// NEWTON-RAPHSON on L  — Ewing Eq. 86 + 87 (unity Le)
//
//   T_chem = (-1 - B') * h_w  +  B'_c * h_c(Tw)  +  B'_g * h_g(Tw)
//   B' = B'_g + B'_c
//
//   q''_cond = rho_ue_CH * (H_r + T_chem)
//            + alpha_w * q_rad_inc
//            - eps_w * sigma * Tw^4
//
//   Denge: q''_cond = k_surf/dx * (Tw - T1)
//
//   f(L) = rho_ue_CH * (H_r + T_chem(Bg,L))
//        + alpha_w * q_rad_inc
//        - eps_w * sigma * Tw(Bg,L)^4
//        - k_surf / dx * (Tw(Bg,L) - T1)
//        = 0
//
//   Tablolardan:
//     Tw(Bg,L)  → lookup_Tw
//     Bc(Bg,L)  → lookup_Bc   (= B'_c)
//     Hw(Bg,L)  → lookup_Hw   (= h_w, J/kg)
//     B'_g      → Bg (giriş parametresi, tablo ekseni)
//     h_c(Tw)   → h_c_T
//     h_g(Tw)   → h_g_T
// =============================================================================
double solve_L_NR(double Bg,
                  double rho_ue_CH,
                  double H_r,
                  double k_surf,
                  double dx_surf,
                  double T1,
                  double emissivity,
                  double alpha_w,
                  double q_rad_inc,
                  double sigma_SB,
                  double L_guess)
{
    const double dL    = 0.01;
    const double L_min = 1.0;
    const double L_max = (double)bpt.Tw[0].size();
    double L = L_guess;

    auto eval_f = [&](double Lx) -> double
    {
        double Tw   = lookup_Tw(Bg, Lx);
        double Bc   = lookup_Bc(Bg, Lx);        // B'_c
        double hw   = lookup_Hw(Bg, Lx);        // h_w  [J/kg]
        double Bp   = Bg + Bc;                  // B' = B'_g + B'_c

        double hc   = h_c_T(Tw);               // char solid enthalpi [J/kg]
        double hg   = h_g_T(Tw);               // pyro gas enthalpi   [J/kg]
        double T_surr     = 300.0;
        // Eq. 87
        double Tchem = (-1.0 - Bp) * hw  +  Bc * hc  +  Bg * hg;

        // Eq. 86 residual
        // Eq. 86 residual
        return rho_ue_CH * (H_r + Tchem)
            - emissivity * sigma_SB * (pow(Tw, 4.0) - pow(T_surr, 4.0))
            - (k_surf / dx_surf) * (Tw - T1);
    };

    for (int iter = 0; iter < 100; iter++)
    {
        double f    = eval_f(L);
        double L2   = min(L + dL, L_max);
        double dfdL = (eval_f(L2) - f) / dL;

        if (fabs(dfdL) < 1e-30) break;

        double step = -f / dfdL;
        step = max(-10.0, min(10.0, step));
        L   += step;
        L    = max(L_min, min(L_max, L));

        if (fabs(step) < 1e-6) break;
    }
    return L;
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
    const double t_end = 60.0;
    const double dt    = 0.001;
    const int    nstep = (int)round(t_end / dt);

    const int    N_nodes  = 1001;
    const double L_domain = 0.05;

    vector<double> x(N_nodes);
    for (int m = 0; m < N_nodes; m++)
        x[m] = m * L_domain / (N_nodes - 1);

    const int N_comp = 3;
    vector<double> B_arr   = {1.200e4,  4.480e9,  0.0};
    vector<double> Psi_arr = {3.0,      3.0,      0.0};
    vector<double> E_arr   = {71.14e6,  169.98e6, 0.0};

    const double Gamma  = 0.5;
    const double R_univ = 8314.0; //J/kmol/K
    const double MW_air = 28.851; //kg/kmol 
    const double R_air  = R_univ / MW_air; // J/kg/K

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

    const double phi_v = 0.8,  phi_c = 0.8;
    const double K_v   = 1.6e-11, K_c = 2.0e-11;

    const double rho_intr_v = Gamma*(rho_v_comp[0]+rho_v_comp[1])
                            + (1.0-Gamma)*rho_v_comp[2];
    const double rho_intr_c = Gamma*(rho_c_comp[0]+rho_c_comp[1])
                            + (1.0-Gamma)*rho_c_comp[2];

    const double rho_virgin     = (1.0-phi_v) * rho_intr_v;
    const double rho_char_total = (1.0-phi_c) * rho_intr_c;

    const double T_back  = 300.0, P_surf = 101325.0, P_back = 101325.0;

    // =========================================================================
    // Case 5 SINIR KOSULU PARAMETRELERİ (Ewing 2013, Sec. IV.E)
    // =========================================================================
    const double rho_ue_CH0 = 0.3;        // enthalpy-based HTC [kg/m2/s], blowing yokken
    const double H_recovery = 2.5e7;      // recovery enthalpy [J/kg]
    const double t_ramp_end = 0.1;        // HTC ramp süresi [s]: 0 -> rho_ue_CH0 lineer

    const double sigma_SB = 5.67e-8;
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
    double h_eff         = 0.0;   // blowing-corrected rho_ue_CH [kg/m2/s]
    double m_dot_surface = 0.0;

    // =========================================================================
    // CIKTI DOSYALARI
    // =========================================================================
    ofstream fout("ablation_history.csv");
    fout << "time,Twall,T1,P0,mdot,rho_ue_CH,sdot_mm_s,recession_mm,thickness_mm,"
            "Bg,Bc,L_nr,eps_surf,k_surf,iter,converged,resid_pct,mdot_g,mdot_c,T_24mm\n";

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
    double L_prev = 1.0;
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

        // =====================================================================
        // 3. dT/dt
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
            dT_dt[i] = (T_old[i] - T_prev[i]) / dt;

        // =====================================================================
        // 4. BASINC
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
            a_P[i] = b_P[i] = c_P[i] = d_P[i] = 0.0;

        a_P[0] = 1.0; d_P[0] = P_surf;
        {
            const int j = N_nodes - 1;
            a_P[j] = 1.0; b_P[j] = 0.0; c_P[j] = 1.0; d_P[j] = 0.0;
        }

        for (int i = 1; i < N_nodes-1; i++)
        {
            p_dxl = x[i]   - x[i-1];
            p_dxr = x[i+1] - x[i];
            p_dxi = 0.5*(p_dxl + p_dxr);

            p_phi   = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i];
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
            m_dot_surface = (P_new[0]/(R_air*T_old[0]))
                          * (K_surf_tmp/mu_g_T(T_old[0]))
                          * (P_new[1]-P_new[0]) / dx_surf;
        }
        // m_dot_surface*=10.0;  // Case 5'te ilk iterasyonda HTC çok düşük geliyor, ramp etkisiyle uyumlu olması için 10x artırarak başlatıyorum  
        if (rho_ue_CH_now > 0.0)
            blowing_factor(m_dot_surface, rho_ue_CH_now, h_eff);
        else
            h_eff = 0.0;
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
                m_dot_surface = (P_new[0]/(R_air*T_wall))
                              * (K_surf2/mu_g_T(T_wall))
                              * (P_new[1]-P_new[0]) / dx_surf;
            }
            if (rho_ue_CH_now > 0.0)
                blowing_factor(m_dot_surface, rho_ue_CH_now, h_eff);
            else
                h_eff = 0.0;
            m_dot_g = m_dot_surface;
            double T_wall_new = T_wall;

            // (a) NR on L — Bg = m_dot_g / rho_ue_CH_now  (enthalpy-based B'g)
            double cp_g_wall = cp_g_T(T_wall);
            Bg_now = m_dot_g / rho_ue_CH_now;
            if (Bg_now < 0.0) Bg_now = 0.0;

            double emissivity_now = eps_surf(T_wall, alpha_eff[0]);
            double L_cur = solve_L_NR(Bg_now, h_eff, H_recovery,
                                    k_surf, dx_surf, T1,
                                    emissivity_now,
                                    1.0,    // alpha_w — gray body
                                    0.0,    // q_rad_inc — Case 5'te yok
                                    sigma_SB,
                                    L_prev);
            L_prev = L_cur;
            T_wall_new = lookup_Tw(Bg_now, L_cur);

            Bc_now = lookup_Bc(Bg_now, L_cur);

            // sdot = B'c * rho_ue_CH / rho_c  (enthalpy-based, Ewing Eq.52)
            sdot = (rho_ue_CH_now > 0.0) ? (Bc_now * h_eff / rho_char_total) : 0.0;
            if (sdot < 0.0) sdot = 0.0;

            double sdot_old_iter = sdot_iter;
            sdot_iter = sdot;

            // (ii) mdot / h_eff tekrar guncelle
            {
                double K_0    = K_v*(1.0-alpha_eff[0]) + K_c*alpha_eff[0];
                double K_1    = K_v*(1.0-alpha_eff[1]) + K_c*alpha_eff[1];
                double K_surf2 = 2.0*K_0*K_1 / (K_0 + K_1);
                m_dot_surface = (P_new[0]/(R_air*T_wall_new))
                              * (K_surf2/mu_g_T(T_wall_new))
                              * (P_new[1]-P_new[0]) / dx_surf;
            }
            if (rho_ue_CH_now > 0.0)
                blowing_factor(m_dot_surface, rho_ue_CH_now, h_eff);
            else
                h_eff = 0.0;
            m_dot_g = m_dot_surface;

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
                t_Spyro  = -(Q_p_T(T_old[i])) * (rho_virgin-rho_char_total) * d_alpha_eff_dt;
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
        double Bc_f   = lookup_Bc(Bg_now, L_prev);
        double hw_f   = lookup_Hw(Bg_now, L_prev);
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
        double recession_step = sdot * dt;
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

        // =====================================================================
        // m_dot_g (pyroliz gaz uretimi entegrali) ve m_dot_c
        // =====================================================================
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
            mdot_g_total += (-drho_dt[i] * dx_i);
        }
        double mdot_c_out = Bc_now * h_eff;   // [kg/m2/s]

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
                 << Bg_now << "," << Bc_now << "," << L_prev << ","
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
    printf("Final rho_ue_CH_eff : %.4f kg/m2s\n", h_eff);
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