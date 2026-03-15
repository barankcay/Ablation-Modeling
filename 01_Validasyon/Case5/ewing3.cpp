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
    bf_phi = 2.0 * 1.0 * m_dot / rho_ue_CH0;

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
    push_row_pyro_gas(tbl_pg, 200,  1512.0, 8.688e-06);
    push_row_pyro_gas(tbl_pg, 325,  1676.0, 1.351e-05);
    push_row_pyro_gas(tbl_pg, 450,  2008.0, 1.796e-05);
    push_row_pyro_gas(tbl_pg, 575,  2867.0, 2.201e-05);
    push_row_pyro_gas(tbl_pg, 700,  6351.0, 2.586e-05);
    push_row_pyro_gas(tbl_pg, 825, 16010.0, 2.993e-05);
    push_row_pyro_gas(tbl_pg, 950, 11650.0, 3.378e-05);
    push_row_pyro_gas(tbl_pg, 1075, 3644.0, 3.683e-05);
    push_row_pyro_gas(tbl_pg, 1200, 6235.0, 3.979e-05);
    push_row_pyro_gas(tbl_pg, 1325, 8707.0, 4.371e-05);
    push_row_pyro_gas(tbl_pg, 1450, 9674.0, 4.708e-05);
    push_row_pyro_gas(tbl_pg, 1575,12530.0, 4.880e-05);
    push_row_pyro_gas(tbl_pg, 1700, 7338.0, 5.021e-05);
    push_row_pyro_gas(tbl_pg, 1825, 4723.0, 5.258e-05);
    push_row_pyro_gas(tbl_pg, 1950, 4142.0, 5.502e-05);
    push_row_pyro_gas(tbl_pg, 2075, 4080.0, 5.744e-05);
    push_row_pyro_gas(tbl_pg, 2200, 4270.0, 5.984e-05);
    push_row_pyro_gas(tbl_pg, 2325, 4670.0, 6.224e-05);
    push_row_pyro_gas(tbl_pg, 2450, 5305.0, 6.464e-05);
    push_row_pyro_gas(tbl_pg, 2575, 6229.0, 6.705e-05);
    push_row_pyro_gas(tbl_pg, 2700, 7513.0, 6.951e-05);
    push_row_pyro_gas(tbl_pg, 2825, 9238.0, 7.201e-05);
    push_row_pyro_gas(tbl_pg, 2950,11500.0, 7.461e-05);
    push_row_pyro_gas(tbl_pg, 3075,14390.0, 7.731e-05);
    push_row_pyro_gas(tbl_pg, 3200,17960.0, 8.015e-05);
    push_row_pyro_gas(tbl_pg, 3325,22090.0, 8.313e-05);
    push_row_pyro_gas(tbl_pg, 3450,26310.0, 8.620e-05);
    push_row_pyro_gas(tbl_pg, 3575,29830.0, 8.917e-05);
    push_row_pyro_gas(tbl_pg, 3700,31960.0, 9.183e-05);
    push_row_pyro_gas(tbl_pg, 3825,32440.0, 9.398e-05);
    push_row_pyro_gas(tbl_pg, 3950,31460.0, 9.561e-05);
}

double cp_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.cp, T); }
double mu_g_T(double T) { return interp1_linear(tbl_pg.T, tbl_pg.mu, T); }

// Piroliz gazi entalpisi: h(T) = integral_T0^T cp(T') dT'  (trapez kuralı)
// h_g_T referans noktasi: tbl_pg.T[0] (T_ref=200K → h=0)
static vector<double> h_g_table;
void init_h_g_table()
{
    const int N = (int)tbl_pg.T.size();
    h_g_table.resize(N, 0.0);
    for (int i = 1; i < N; i++)
    {
        double dT     = tbl_pg.T[i] - tbl_pg.T[i-1];
        double cp_avg = 0.5*(tbl_pg.cp[i] + tbl_pg.cp[i-1]);
        h_g_table[i]  = h_g_table[i-1] + cp_avg * dT;
    }
}
double h_g_T(double T) { return interp1_linear(tbl_pg.T, h_g_table, T); }

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
    vector<vector<double>> Tchem;
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

        double bg, tw, bc, tchem;
        if (sscanf(line.c_str(), "%lf %lf %lf %lf", &bg, &tw, &bc, &tchem) != 4)
            continue;

        if (bg != prev_bg)
        {
            bpt.bg.push_back(bg);
            bpt.Tw.push_back(vector<double>());
            bpt.Bc.push_back(vector<double>());
            bpt.Tchem.push_back(vector<double>());
            ibg++;
            prev_bg = bg;
        }

        bpt.Tw[ibg].push_back(tw);
        bpt.Bc[ibg].push_back(bc);
        bpt.Tchem[ibg].push_back(tchem);
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

double lookup_Tchem(double Bg, double L)
{
    int i0, i1; double t;
    bg_bracket(Bg, i0, i1, t);
    return lerp(row_interp(bpt.Tchem[i0], L), row_interp(bpt.Tchem[i1], L), t);
}

// =============================================================================
// NEWTON-RAPHSON on L  (Ewing Eq.86-90)
// Case 5: enthalpy-based enerji dengesi
//   q_cond = rho_ue_CH*(H_r - h_w) + eps*sigma*(T_surr^4 - Tw^4) + rho_ue_CH*Tchem/cp_g
//   h_w = cp_g * Tw  (simplified, unity Le)
//   f(L) = k_surf/dx*(Tw-T1) - rho_ue_CH*(H_r - cp_g*Tw)
//          - eps*sigma*(T_surr^4-Tw^4) - rho_ue_CH*Tchem
// =============================================================================
double solve_L_NR(double Bg, double rho_ue_CH, double k_surf,
                  double T1, double H_recovery,
                  double emissivity, double sigma_SB, double T_surr,
                  double dx_surf, double cp_g, double L_guess)
{
    double L = L_guess;

    for (int iter = 0; iter < 100; iter++)
    {
        double Tw    = lookup_Tw(Bg, L);
        double Tchem = lookup_Tchem(Bg, L);
        double h_w   = h_g_T(Tw);   // entegre entalpi: h(T) = integral cp dT

        double f = (k_surf/dx_surf)*(Tw - T1)
                 - rho_ue_CH*(H_recovery - h_w)
                 - emissivity*sigma_SB*(pow(T_surr,4) - pow(Tw,4));

        double dL   = 0.01;
        double Tw2  = lookup_Tw(Bg, L + dL);
        double Tc2  = lookup_Tchem(Bg, L + dL);
        double h_w2 = h_g_T(Tw2);

        double f2 = (k_surf/dx_surf)*(Tw2 - T1)
                  - rho_ue_CH*(H_recovery - h_w2)
                  - emissivity*sigma_SB*(pow(T_surr,4) - pow(Tw2,4));

        double dfdL = (f2 - f) / dL;
        if (fabs(dfdL) < 1e-30) break;

        double step = -f / dfdL;
        if (step >  10.0) step =  10.0;
        if (step < -10.0) step = -10.0;
        L += step;
        if (L < 1.0)                       L = 1.0;
        if (L > (double)bpt.Tw[0].size())  L = (double)bpt.Tw[0].size();

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
    const double dt    = 0.01;
    const int    nstep = (int)round(t_end / dt);

    const int    N_nodes  = 501;
    const double L_domain = 0.05;

    vector<double> x(N_nodes);
    for (int m = 0; m < N_nodes; m++)
        x[m] = m * L_domain / (N_nodes - 1);

    const int N_comp = 3;
    vector<double> B_arr   = {1.200e4,  4.480e9,  0.0};
    vector<double> Psi_arr = {3.0,      3.0,      0.0};
    vector<double> E_arr   = {71.14e6,  169.98e6, 0.0};

    const double Gamma  = 0.5;
    const double R_univ = 8314.0;
    const double MW_air = 28.97;
    const double R_air  = R_univ / MW_air;

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

    const double phi_v = 0.8,  phi_c = 0.85;
    const double K_v   = 1.6e-11, K_c = 2.0e-11;

    const double rho_intr_v = Gamma*(rho_v_comp[0]+rho_v_comp[1])
                            + (1.0-Gamma)*rho_v_comp[2];
    const double rho_intr_c = Gamma*(rho_c_comp[0]+rho_c_comp[1])
                            + (1.0-Gamma)*rho_c_comp[2];

    const double rho_virgin     = (1.0-phi_v) * rho_intr_v;
    const double rho_char_total = (1.0-phi_c) * rho_intr_c;

    const double T_back  = 298.0, P_surf = 101325.0, P_back = 101325.0;

    // =========================================================================
    // Case 5 SINIR KOSULU PARAMETRELERİ (Ewing 2013, Sec. IV.E)
    // =========================================================================
    const double rho_ue_CH0 = 0.3;        // enthalpy-based HTC [kg/m2/s], blowing yokken
    const double H_recovery = 2.5e7;      // recovery enthalpy [J/kg]
    const double t_ramp_end = 0.1;        // HTC ramp süresi [s]: 0 -> rho_ue_CH0 lineer

    const double sigma_SB = 5.67e-8;
    const double T_surr   = 300.0;

    const int    max_iter = 40;
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
    init_h_g_table();
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
            "Bg,Bc,L_nr,eps_surf,k_surf,iter,converged,resid_pct\n";

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
            for (int c = 0; c < N_comp; c++)
                drho_dt[i] += (rho_v_comp[c]-rho_c_comp[c])
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
        if (rho_ue_CH_now > 0.0)
            blowing_factor(m_dot_surface + Bc_now * rho_ue_CH_now, rho_ue_CH_now, h_eff);
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
                blowing_factor(m_dot_surface + Bc_now * rho_ue_CH_now, rho_ue_CH_now, h_eff);
            else
                h_eff = 0.0;
            m_dot_g = m_dot_surface;
            double T_wall_new = T_wall;

            // (a) NR on L — Bg = m_dot_g / rho_ue_CH_now  (enthalpy-based B'g)
            double cp_g_wall = cp_g_T(T_wall);
            Bg_now = m_dot_g / rho_ue_CH0;
            if (Bg_now < 0.0) Bg_now = 0.0;

            double L_cur = solve_L_NR(Bg_now, h_eff, k_surf, T1,
                                      H_recovery, 0.0, sigma_SB,
                                      T_surr, dx_surf, cp_g_wall, L_prev);
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
                blowing_factor(m_dot_surface + Bc_now * rho_ue_CH_now, rho_ue_CH_now, h_eff);
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
        double h_w_now   = h_g_T(T_wall);
        double q_conv_s  = h_eff * (H_recovery - h_w_now);
        double q_rad_s   = 0.0;   // Case 5: radyasyon yok
        double q_chem_s  = h_eff * lookup_Tchem(Bg_now, L_prev);
        double q_cond_s  = k_surf * (T_wall - T_new[1]) / (x[1] - x[0]);
        double resid_s   = q_conv_s + q_rad_s - q_cond_s;   // Tchem diagnostik, NR'de yok
        double resid_pct_s = (fabs(q_conv_s) > 1e-10) ? (resid_s / q_conv_s) * 100.0 : 0.0;

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
                 << resid_pct_s << "\n";
        }
        if ((n % (save_every*10)) == 0 || n == 1 || n == nstep)
        {
            printf("t=%7.1fs | Tw=%7.1fK | T1=%7.1fK | rho_ue_CH=%6.4f | "
                   "Bg=%.3f Bc=%.4f L=%5.2f eps=%.3f k=%.3f | "
                   "sdot=%7.4fmm/s | erim=%7.4fmm | "
                   "it=%d %s | resid=%+.3f%%\n",
                   time, T_wall, T_old[1], h_eff,
                   Bg_now, Bc_now, L_prev,
                   eps_surf(T_wall, alpha_eff[0]), k_surf,
                   sdot*1e3, recession_total*1e3,
                   iter_count,
                   iter_converged ? "[CONV]" : "[WARN:MAX_ITER]",
                   resid_pct_s);
        }
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

    double h_w_f  = h_g_T(T_wall);
    double q_conv = h_eff * (H_recovery - h_w_f);
    double q_rad  = 0.0;   // Case 5: radyasyon yok
    double q_chem = h_eff * lookup_Tchem(Bg_now, L_prev);
    double q_cond = k_surf * (T_wall - T_old[1]) / (x[1]-x[0]);
    double resid  = q_conv + q_rad - q_cond;   // Tchem diagnostik, NR'de yok

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

    return 0;
}