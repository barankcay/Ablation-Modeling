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

// distance-weighted linear average (distance must be positive)
double inverseAverage(double v1, double d1, double v2, double d2)
{
    d1 = (fabs(d1));
    d2 = (fabs(d2));
    return (v1 * d2 + v2 * d1) / (d1 + d2);
}

// =============================================================================
// BLOWING FACTOR
// Omega = phi / (exp(phi) - 1)
// =============================================================================
double blowing_factor(double m_dot, double rho_e, double u_e,
                      double h0, double cp_g, double& h_eff_out)
{
    double C_h0  = h0 / (rho_e * u_e * cp_g);
    double denom = rho_e * u_e * C_h0;

    bf_phi = 2.0 * 0.5 * m_dot / denom;

    // phi ~ 0 yakınında 0/0 sayısal sorun olmasın
    if (fabs(bf_phi) < 1e-8)
        bf_Omega = 1.0 - 0.5 * bf_phi;  // 1. mertebe Taylor
    else
        bf_Omega = bf_phi / (exp(bf_phi) - 1.0);

    h_eff_out = bf_Omega * h0;
    return bf_Omega;
}

// =============================================================================
// PATANKAR TDMA
// a[i]*T[i] = b[i]*T[i+1] + c[i]*T[i-1] + d[i]
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
// =============================================================================
// =============================================================================
// ============================================================
// TERMOFIZIKSEL TABLOLAR  (thermalProperties.csv - SI)
// Virgin eps=0.8  /  Char eps=0.9
// ============================================================
struct ThermalTable{vector<double>T,cp,k,eps;}tbl_v,tbl_c;
void init_virgin_table(){
    double d[][4]={{256,879,0.398,0.8},{298,984,0.403,0.8},
        {444,1300,0.416,0.8},{556,1470,0.453,0.8},{644,1570,0.470,0.8},
        {833,1720,0.486,0.8},{1111,1860,0.523,0.8},{1389,1930,0.560,0.8},
        {1667,1980,0.698,0.8},{1944,1990,0.872,0.8},{2222,2000,1.110,0.8},
        {2778,2010,1.750,0.8},{3333,2010,2.780,0.8}};
    for(auto&r:d){tbl_v.T.push_back(r[0]);tbl_v.cp.push_back(r[1]);
                  tbl_v.k.push_back(r[2]);tbl_v.eps.push_back(r[3]);}}
void init_char_table(){
    double d[][4]={{256,733,0.398,0.9},{298,783,0.403,0.9},
        {444,1090,0.416,0.9},{556,1320,0.453,0.9},{644,1430,0.470,0.9},
        {833,1680,0.486,0.9},{1111,1840,0.523,0.9},{1389,1970,0.560,0.9},
        {1667,2050,0.605,0.9},{1944,2090,0.729,0.9},{2222,2110,0.922,0.9},
        {2778,2140,1.460,0.9},{3333,2150,2.320,0.9}};
    for(auto&r:d){tbl_c.T.push_back(r[0]);tbl_c.cp.push_back(r[1]);
                  tbl_c.k.push_back(r[2]);tbl_c.eps.push_back(r[3]);}}
double interp1(const vector<double>&xv,const vector<double>&yv,double x){
    int N=(int)xv.size();
    if(x<=xv[0])return yv[0];if(x>=xv[N-1])return yv[N-1];
    for(int i=0;i<N-1;i++)if(x<xv[i+1]){
        double t=(x-xv[i])/(xv[i+1]-xv[i]);return yv[i]+t*(yv[i+1]-yv[i]);}
    return yv[N-1];}
double cp_v_T(double T){return interp1(tbl_v.T,tbl_v.cp,T);}
double k_v_T (double T){return interp1(tbl_v.T,tbl_v.k, T);}
double cp_c_T(double T){return interp1(tbl_c.T,tbl_c.cp,T);}
double k_c_T (double T){return interp1(tbl_c.T,tbl_c.k, T);}
double cp_mix(double T,double a){return cp_v_T(T)*(1.0-a)+cp_c_T(T)*a;}
double k_mix (double T,double a){return k_v_T(T)*(1.0-a) +k_c_T(T)*a;}
double eps_surf(double T,double a){
    return interp1(tbl_v.T,tbl_v.eps,T)*(1.0-a)
          +interp1(tbl_c.T,tbl_c.eps,T)*a;}

// ============================================================
// PIROLIZ GAZI TABLOLARI  (pyrolysisModel.csv - 1 atm, SI)
// cp_g [J/kgK]  ve  mu_g [Pa.s]  T-bagli
// ============================================================
struct PyroGasTable{vector<double>T,cp,mu;}tbl_pg;
void init_pyrogas_table(){
    double d[][3]={
        {200,1512.0,8.688e-06},{325,1676.0,1.351e-05},{450,2008.0,1.796e-05},
        {575,2867.0,2.201e-05},{700,6351.0,2.586e-05},{825,16010.0,2.993e-05},
        {950,11650.0,3.378e-05},{1075,3644.0,3.683e-05},{1200,6235.0,3.979e-05},
        {1325,8707.0,4.371e-05},{1450,9674.0,4.708e-05},{1575,12530.0,4.880e-05},
        {1700,7338.0,5.021e-05},{1825,4723.0,5.258e-05},{1950,4142.0,5.502e-05},
        {2075,4080.0,5.744e-05},{2200,4270.0,5.984e-05},{2325,4670.0,6.224e-05},
        {2450,5305.0,6.464e-05},{2575,6229.0,6.705e-05},{2700,7513.0,6.951e-05},
        {2825,9238.0,7.201e-05},{2950,11500.0,7.461e-05},{3075,14390.0,7.731e-05},
        {3200,17960.0,8.015e-05},{3325,22090.0,8.313e-05},{3450,26310.0,8.620e-05},
        {3575,29830.0,8.917e-05},{3700,31960.0,9.183e-05},{3825,32440.0,9.398e-05},
        {3950,31460.0,9.561e-05}};
    for(auto&r:d){tbl_pg.T.push_back(r[0]);tbl_pg.cp.push_back(r[1]);
                  tbl_pg.mu.push_back(r[2]);}}
double cp_g_T(double T){return interp1(tbl_pg.T,tbl_pg.cp,T);}
double mu_g_T(double T){return interp1(tbl_pg.T,tbl_pg.mu,T);}

// BPRIME TABLOSU — bprime_table.txt dosyasindan yuklenir
// Format: Bg   Tw[K]   Bc[-]   Tchem[J/kg]   (# yorum, bos satirlar atlanir)
// =============================================================================
struct BprimeTable {
    vector<double>         bg_vals;
    vector<vector<double>> Tw, Bc, Tchem;
} bpt;

void load_bprime_table(const string& fname)
{
    ifstream fin(fname);
    if (!fin.is_open()) {
        printf("ERROR: '%s' acilamiyor!\n", fname.c_str()); exit(1);
    }
    double prev_bg=-1e99; int ibg=-1;
    string line;
    while (getline(fin,line)) {
        if (line.empty()||line[0]=='#') continue;
        double bg,tw,bc,tchem;
        if (sscanf(line.c_str(),"%lf %lf %lf %lf",&bg,&tw,&bc,&tchem)!=4) continue;
        if (bg!=prev_bg) {
            bpt.bg_vals.push_back(bg);
            bpt.Tw.push_back({}); bpt.Bc.push_back({}); bpt.Tchem.push_back({});
            ibg++; prev_bg=bg;
        }
        bpt.Tw[ibg].push_back(tw);
        bpt.Bc[ibg].push_back(bc);
        bpt.Tchem[ibg].push_back(tchem);
    }
    printf("[TABLE] %s: %d Bg grubu, %d satir/grup, Tw=[%.0f..%.0f]K\n",
           fname.c_str(),(int)bpt.bg_vals.size(),(int)bpt.Tw[0].size(),
           bpt.Tw[0].front(),bpt.Tw[0].back());
}

void find_bg_bracket(double Bg, int& ilo, int& ihi, double& frac)
{
    int N=(int)bpt.bg_vals.size();
    if (Bg<=bpt.bg_vals[0])   {ilo=0;  ihi=0;  frac=0.0;return;}
    if (Bg>=bpt.bg_vals[N-1]) {ilo=N-1;ihi=N-1;frac=0.0;return;}
    for (int i=0;i<N-1;i++)
        if (bpt.bg_vals[i+1]>=Bg) {
            ilo=i; ihi=i+1;
            frac=(Bg-bpt.bg_vals[i])/(bpt.bg_vals[i+1]-bpt.bg_vals[i]);
            return;
        }
}

double interp_row(const vector<double>& v, double L)
{
    int N=(int)v.size();
    double idx=max(0.0,min((double)(N-1),L-1.0));
    int jlo=(int)idx, jhi=min(jlo+1,N-1);
    return v[jlo]+(idx-jlo)*(v[jhi]-v[jlo]);
}

double lookup_Tw(double Bg, double L) {
    int ilo,ihi;double frac;find_bg_bracket(Bg,ilo,ihi,frac);
    return interp_row(bpt.Tw[ilo],L)+frac*(interp_row(bpt.Tw[ihi],L)-interp_row(bpt.Tw[ilo],L));
}
double lookup_Bc(double Bg, double L) {
    int ilo,ihi;double frac;find_bg_bracket(Bg,ilo,ihi,frac);
    return interp_row(bpt.Bc[ilo],L)+frac*(interp_row(bpt.Bc[ihi],L)-interp_row(bpt.Bc[ilo],L));
}
double lookup_Tchem(double Bg, double L) {
    int ilo,ihi;double frac;find_bg_bracket(Bg,ilo,ihi,frac);
    return interp_row(bpt.Tchem[ilo],L)+frac*(interp_row(bpt.Tchem[ihi],L)-interp_row(bpt.Tchem[ilo],L));
}


// =============================================================================
// NEWTON-RAPHSON on L  (Ewing Eq.86-90)
//
// solve_Twall_NR ile birebir ayni yapi.
// Degisken: Tw yerine L (satir numarasi).
// Tw(L) ve Tchem(L) tablodan okunuyor.
//
//   f(L) = k/dx*(Tw(L) - T1)
//         - h_eff*(T_rec - Tw(L))
//         - eps*sigma*(T_surr^4 - Tw(L)^4)
//         - (h_eff/cp_g)*Tchem(L)
//        = 0
// =============================================================================
double solve_L_NR(double Bg, double h_eff, double k_surf,
                  double T1, double T_recovery,
                  double emissivity, double sigma_SB, double T_surr,
                  double dx_surf, double cp_g, double L_guess)
{
    double L = L_guess;

    for (int iter = 0; iter < 100; iter++)
    {
        double Tw = lookup_Tw(Bg, L);
        double Tc = lookup_Tchem(Bg, L);

        double f = (k_surf/dx_surf)*(Tw - T1)
                 - h_eff*(T_recovery - Tw)
                 - emissivity*sigma_SB*(pow(T_surr,4) - pow(Tw,4))
                 - (h_eff/cp_g)*Tc;

        double dL  = 0.01;
        double Tw2 = lookup_Tw(Bg, L + dL);
        double Tc2 = lookup_Tchem(Bg, L + dL);

        double f2  = (k_surf/dx_surf)*(Tw2 - T1)
                   - h_eff*(T_recovery - Tw2)
                   - emissivity*sigma_SB*(pow(T_surr,4) - pow(Tw2,4))
                   - (h_eff/cp_g)*Tc2;

        double dfdL = (f2 - f) / dL;
        if (fabs(dfdL) < 1e-30) break;

        double step = -f / dfdL;
        if (step >  10.0) step =  10.0;
        if (step < -10.0) step = -10.0;
        L += step;
        if (L < 1.0)            L = 1.0;
        if (L > (double)bpt.Tw[0].size()) L = (double)bpt.Tw[0].size();

        if (fabs(step) < 1e-6) break;
    }
    return L;
}


// =============================================================================
// MAIN
// =============================================================================
int main()
{
    cout << string(80, '=') << "\n";
    cout << "ABLASYON ISI TRANSFERI v15 — Ewing NR-L + T-bagli termofizik + transport CSV\n";
    cout << string(80, '=') << "\n";

    // =========================================================================
    // PARAMETRELER
    // =========================================================================
    const double t_end = 500.0;
    const double dt    = 0.01;
    const int    nstep = (int)round(t_end / dt);

    const int    N_nodes  = 100;
    const double L_domain = 0.05;

    vector<double> x(N_nodes);
    for (int m = 0; m < N_nodes; m++)
        x[m] = m * L_domain / (N_nodes - 1);

    const int N_comp = 3;
    vector<double> B_arr   = {1.40e4,   9.75e8,   0.0};
    vector<double> Psi_arr = {3.0,      3.0,      0.0};
    vector<double> E_arr   = {71.14e6,  169.98e6, 0.0};
    const double R_univ    = 8.314;

    vector<double> E_over_R(N_comp);
    for (int c = 0; c < N_comp; c++)
        E_over_R[c] = (E_arr[c] / 1000.0) / R_univ;

    vector<double> rho_v_comp = {325.015,  973.926,  2066.380};
    vector<double> rho_c_comp = {0.0,      518.998,  2066.380};
    const double Gamma = 0.422;

    double rho_virgin     = Gamma*(rho_v_comp[0]+rho_v_comp[1])
                          + (1-Gamma)*rho_v_comp[2];
    double rho_char_total = Gamma*(rho_c_comp[0]+rho_c_comp[1])
                          + (1-Gamma)*rho_c_comp[2];

    // k ve cp: T-bagli tablolardan (init_virgin/char_table)
    const double Q_p   = 0.5e7;

    // transportProperties.csv
    const double phi_v      = 0.8,   phi_c      = 0.85;
    const double K_v        = 1.6e-11, K_c = 2.0e-11;  // permeability [m^2]
    // pyrolysisModel.csv -> T-bagli (cp_g_T, mu_g_T fonksiyonlari)

    const double T_back = 300.0, P_surf = 101325.0, P_back = 101325.0;

    const double h_0_external = 200.0;
    const double T_recovery   = 8000.0;
    const double rho_e        = 1.2,  u_e = 1500.0;

    const double sigma_SB     = 5.67e-8;
    const double T_surr = 300.0;

    // iterasyon parametreleri
    const int    max_iter  = 40;
    const double eps_T     = 0.1;   // K
    const double eps_sdot  = 1e-9;   // m/s


    // --------------------------------

    printf("\n[INIT] dt=%.4f s, t_end=%.1f s, adim=%d\n", dt, t_end, nstep);
    printf("  rho_virgin=%.1f  rho_char=%.1f kg/m3\n", rho_virgin, rho_char_total);
    printf("  dx=%.3f mm, N=%d\n", (x[1]-x[0])*1e3, N_nodes);

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

    double T_wall        = T_back;
    double m_dot_g       = 0.0;
    double k_surf        ;
    double h_eff         = h_0_external;
    double m_dot_surface = 0.0;

    // =========================================================================
    // CSV
    // =========================================================================
    ofstream fout("ablation_history.csv");
    // debug kolonları eklendi: Tw_stat, sdot, dP01, dx_surf, k_surf
    fout << "time,Twall,T1,P0,mdot,heff,sdot_mm_s,recession_mm,thickness_mm,Bg,Bc,L_nr,eps_surf,k_surf,iter\n";
    const int save_every = 2;

    printf("\n[LOOP] Basliyor...\n\n");

    // =========================================================================
    // LOOP DISI VEKTORLER
    // =========================================================================
    vector<vector<double>> alpha_new(N_comp, vector<double>(N_nodes, 0.0));
    vector<double> rho_solid_new(N_nodes);
    vector<double> alpha_eff(N_nodes);
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

    // =========================================================================
    // TIME LOOP
    // =========================================================================
    for (int n = 1; n <= nstep; n++)
    {
        double time    = n * dt;
        double dx_surf = x[1] - x[0];

        // =====================================================================
        // 1. PIROLIZ
        // =====================================================================
        for (int c = 0; c < N_comp; c++)
        {
            if (B_arr[c] == 0.0) { alpha_new[c] = alpha_old[c]; continue; }
            for (int i = 0; i < N_nodes; i++)
            {
                pyr_I   = B_arr[c] * exp(-E_over_R[c] / T_old[i]) * dt;
                pyr_ai  = alpha_old[c][i];
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

        for (int i = 0; i < N_nodes; i++)
        {
            rho_solid_new[i] =
                Gamma * ( (rho_v_comp[0]-(rho_v_comp[0]-rho_c_comp[0])*alpha_new[0][i])
                         +(rho_v_comp[1]-(rho_v_comp[1]-rho_c_comp[1])*alpha_new[1][i]) )
              + (1-Gamma)*( rho_v_comp[2]-(rho_v_comp[2]-rho_c_comp[2])*alpha_new[2][i] );

            alpha_eff[i] = (rho_virgin - rho_solid_new[i]) / (rho_virgin - rho_char_total);

            drho_dt[i] = 0.0;
            for (int c = 0; c < N_comp; c++)
                drho_dt[i] += (rho_v_comp[c]-rho_c_comp[c])
                            * (alpha_new[c][i]-alpha_old[c][i]) / dt;

            k_node[i] = k_mix(T_old[i],alpha_eff[i]);
        }

        // =====================================================================
        // 2. dT/dt (pressure source için)
        // Not: T_prev her step sonunda T_old yapılacak (daha temiz)
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
            dT_dt[i] = (T_old[i] - T_prev[i]) / dt;

        // =====================================================================
        // 3. BASINC
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

            p_phi   = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];
            p_atime = p_phi / (R_univ*T_old[i]) * p_dxi / dt;
            double K_i = K_v*(1.0-alpha_eff[i])+K_c*alpha_eff[i];
            p_K     = P_old[i] / (R_univ*T_old[i]) * K_i/mu_g_T(T_old[i]);
            p_ae    = p_K / p_dxr;
            p_aw    = p_K / p_dxl;

            p_aeff_old = (rho_virgin-rho_solid_old[i]) / (rho_virgin-rho_char_total);
            p_dalpha   = (alpha_eff[i] - p_aeff_old) / dt;
            p_Sptemp   = p_phi / (R_univ*T_old[i]*T_old[i]) * dT_dt[i];
            p_Spporo   = (phi_v-phi_c) / (R_univ*T_old[i]) * p_dalpha;
            p_Slin     = p_Sptemp + p_Spporo;

            a_P[i] = p_atime + p_ae + p_aw - p_Slin*p_dxi;
            b_P[i] = p_ae;
            c_P[i] = p_aw;
            d_P[i] = drho_dt[i]*p_dxi + p_atime*P_old[i];
        }

        P_new = thomas_patankar(a_P, b_P, c_P, d_P);

        // =====================================================================
        // 4. BLOWING (dış iterasyon dışında ilk tahmin)
        // =====================================================================
        m_dot_surface = (P_new[0]/(R_univ*T_old[0]))
                      * (K_v/mu_g_T(T_old[0])) * (P_new[1]-P_new[0]) / dx_surf;
        blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g_T(T_old[0]), h_eff);
        m_dot_g = m_dot_surface;

        // =====================================================================
        // 5. ITERATIF YUZEY + SICAKLIK COZUMU
        // Patch:
        //  - histerezis (dT_hys)
        //  - sdot under-relax (ur_sdot)
        //  - mdot/heff iterasyon içinde Twall_new ile güncelle (coupling)
        // =====================================================================
        k_surf = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);

        double T1         = T_old[1];
        double sdot_iter  = sdot;
        int    iter_count = 0;

        for (int it = 0; it < max_iter; it++)
        {
            iter_count = it + 1;

            // --- (i) iterasyon içi: mdot/heff güncelle (Twall’a bağlı yüzey yoğunluğu) ---
            // önce Twall tahmini olarak mevcut T_wall kullan
            double Tw_for_rho = T_wall;
            m_dot_surface = (P_new[0]/(R_univ*Tw_for_rho))
                          * (K_v/mu_g_T(T_wall)) * (P_new[1]-P_new[0]) / dx_surf;
            blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g_T(T_wall), h_eff);
            m_dot_g = m_dot_surface;

            // --- a) EWING: NR on L ---
            double T_wall_new = T_wall;

            // B'g = m_dot_g / (rho_e*u_e*C_M),  C_M = h_0/(rho_e*u_e*cp_g)
            double cp_g_wall=cp_g_T(T_wall);
            Bg_now = m_dot_g * cp_g_wall / h_0_external;
            if (Bg_now < 0.0) Bg_now = 0.0;

            // solve_L_NR: orijinal NR gibi, degisken L
            double emissivity_now = eps_surf(T_wall,alpha_eff[0]);
            double L_cur = solve_L_NR(Bg_now, h_eff, k_surf, T1,
                                      T_recovery, emissivity_now, sigma_SB,
                                      T_surr, dx_surf, cp_g_wall, L_prev);
            L_prev = L_cur;

            T_wall_new = lookup_Tw(Bg_now, L_cur);
            Bc_now     = lookup_Bc(Bg_now, L_cur);

            // sdot = B'c * rho_e*u_e*C_M / rho_char  (Ewing Eq.47)
            sdot = Bc_now * (h_0_external / cp_g_wall) / rho_char_total;
            if (sdot < 0.0) sdot = 0.0;
            // --- (ii) sdot under-relax ---
            double sdot_old_iter = sdot_iter;
            sdot_iter = sdot;
            

            // --- (iii) mdot/heff’i yeni Twall ile tekrar bağla (çok pahalı değil) ---
            Tw_for_rho = T_wall_new;
            m_dot_surface = (P_new[0]/(R_univ*Tw_for_rho))
                          * (K_v/mu_g_T(T_wall_new)) * (P_new[1]-P_new[0]) / dx_surf;
            blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g_T(T_wall_new), h_eff);
            m_dot_g = m_dot_surface;

            // --- b) SICAKLIK TDMA ---
            for (int i = 0; i < N_nodes; i++)
            {
                a_T[i] = 0.0; b_T[i] = 0.0; c_T[i] = 0.0; d_T[i] = 0.0;
            }
            a_T[0] = 1.0; d_T[0] = T_wall_new;

            for (int i = 1; i < N_nodes-1; i++)
            {
                t_dxl = x[i]   - x[i-1];
                t_dxr = x[i+1] - x[i];
                t_dxi = 0.5*(t_dxl + t_dxr);

                t_cp  = cp_mix(T_old[i],alpha_eff[i]);
                t_phi = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];

                t_ke = 2.0*k_node[i]*k_node[i+1] / (k_node[i]+k_node[i+1]);
                t_kw = 2.0*k_node[i]*k_node[i-1] / (k_node[i]+k_node[i-1]);

                t_rhogas = P_new[i] / (R_univ*T_old[i]);
                t_rhoc   = rho_solid_new[i]*t_cp + t_rhogas*t_phi*cp_g_T(T_old[i]);
                t_atime  = t_rhoc * t_dxi / dt;
                t_ae     = t_ke / t_dxr;
                t_aw     = t_kw / t_dxl;

                double K_i2 = K_v*(1.0-alpha_eff[i])+K_c*alpha_eff[i];
                t_mdotgas = t_rhogas * (K_i2/mu_g_T(T_old[i]))
                          * (P_new[i+1]-P_new[i-1]) / (t_dxl+t_dxr);
                t_hgrad   = cp_g_T(T_old[i]) * (T_old[i+1]-T_old[i-1]) / (t_dxl+t_dxr);

                t_hsolid = t_cp * T_old[i];
                t_hpyro  = t_hsolid + rho_solid_new[i]*(cp_v_T(T_old[i])-cp_c_T(T_old[i]))*T_old[i]
                                     / (rho_virgin-rho_char_total);
                t_Spyro  = -(Q_p - t_hpyro + cp_g_T(T_old[i])*T_old[i]) * drho_dt[i];
                t_Stotal = t_Spyro - t_mdotgas*t_hgrad;

                a_T[i] = t_atime + t_ae + t_aw;
                b_T[i] = t_ae;
                c_T[i] = t_aw;
                d_T[i] = t_Stotal*t_dxi + t_atime*T_old[i];
            }

            // arka yuzey: adiabatik
            {
                int j    = N_nodes - 1;
                t_dxl    = x[j] - x[j-1];
                t_cp     = cp_mix(T_old[j],alpha_eff[j]);
                t_phi    = phi_v*(1-alpha_eff[j]) + phi_c*alpha_eff[j];
                t_rhogas = P_new[j] / (R_univ*T_old[j]);
                t_rhoc   = rho_solid_new[j]*t_cp + t_rhogas*t_phi*cp_g_T(T_old[j]);
                t_atime  = t_rhoc * (0.5*t_dxl) / dt;
                t_kw     = 2.0*k_node[j]*k_node[j-1] / (k_node[j]+k_node[j-1]);
                a_T[j]   = t_atime + t_kw/t_dxl;
                b_T[j]   = 0.0;
                c_T[j]   = t_kw/t_dxl;
                d_T[j]   = t_atime * T_old[j];
            }

            T_new = thomas_patankar(a_T, b_T, c_T, d_T);

            // --- c) T1 guncelle ---
            double T1_new = T_new[1];

            // --- d) convergence ---
            double dT1   = fabs(T1_new   - T1);
            double dsdot = fabs(sdot_iter - sdot_old_iter);

            T_wall = T_wall_new;
            T1     = T1_new;

            if (dT1 < eps_T && dsdot < eps_sdot) break;
        }


        // =====================================================================
        // 6. MOVING BOUNDARY + REMAP
        // Patch:
        //  - remap distance fabs+clamp (inverseAverage zaten bunu yapıyor)
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

            for (int j = 0; j < N_nodes; j++)
                x[j] += recession_step * (N_nodes-1-j) / (N_nodes-1);

            for (int j = 0; j < N_nodes-1; j++)
            {
                double d1 = x[j] - x_old[j];
                double d2 = x_old[j+1] - x[j];

                T_new[j]         = inverseAverage(T_new[j],         d1, T_new[j+1],         d2);
                P_new[j]         = inverseAverage(P_new[j],         d1, P_new[j+1],         d2);
                rho_solid_new[j] = inverseAverage(rho_solid_new[j], d1, rho_solid_new[j+1], d2);

                for (int c = 0; c < N_comp; c++)
                    alpha_new[c][j] = inverseAverage(alpha_new[c][j], d1, alpha_new[c][j+1], d2);
            }

            T_new[0] = T_wall;
            P_new[0] = P_surf;
            // k_node'u remap sonrası alpha_new ile güncelle
            for (int j = 0; j < N_nodes; j++)
            {
                rho_solid_new[j] =
                    Gamma * ( (rho_v_comp[0]-(rho_v_comp[0]-rho_c_comp[0])*alpha_new[0][j])
                            +(rho_v_comp[1]-(rho_v_comp[1]-rho_c_comp[1])*alpha_new[1][j]) )
                + (1-Gamma)*( rho_v_comp[2]-(rho_v_comp[2]-rho_c_comp[2])*alpha_new[2][j] );
                alpha_eff[j] = (rho_virgin - rho_solid_new[j]) / (rho_virgin - rho_char_total);
                k_node[j]    = k_mix(T_old[j],alpha_eff[j]);
            }
            k_surf = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);
        }

        // =====================================================================
        // 7. UPDATE
        // Patch:
        //  - T_prev her step sonunda T_old yapılır (dT/dt tanımı net)
        // =====================================================================
        T_prev = T_old;      // <<< PATCH (önceki karışıklığı kaldır)

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
                 << eps_surf(T_wall,alpha_eff[0]) << "," << k_surf << ","
                 << iter_count
                 << "\n";
        }

        if ((n % (save_every*10)) == 0 || n == 1 || n == nstep)
        {
            printf("t=%7.1fs | Tw=%7.1fK | T1=%7.1fK | heff=%6.1f | "
                   "Bg=%.3f Bc=%.4f L=%5.2f eps=%.3f k=%.3f | "
                   "sdot=%7.4fmm/s | erim=%7.4fmm | it=%d\n",
                   time, T_wall, T_old[1], h_eff,
                   Bg_now, Bc_now, L_prev,
                   eps_surf(T_wall,alpha_eff[0]), k_surf,
                   sdot*1e3, recession_total*1e3, iter_count);
        }
    }

    fout.close();

    printf("\n%s\n", string(80, '=').c_str());
    printf("TAMAMLANDI\n");
    printf("Final Twall       : %.2f K\n", T_wall);
    printf("Final T[0]        : %.2f K\n", T_old[0]);
    printf("Final T[1]        : %.2f K\n", T_old[1]);
    printf("Final P[0]        : %.3f kPa\n", P_old[0]/1000.0);
    printf("Final mdot        : %.3e kg/m2s\n", m_dot_surface);
    printf("Final h_eff       : %.2f W/m2K\n", h_eff);
    printf("Final s_dot       : %.4f mm/s\n", sdot*1e3);
    printf("Final Bg          : %.4f\n", Bg_now);
    printf("Final Bc          : %.4f\n", Bc_now);
    printf("Toplam erime      : %.4f mm\n", recession_total*1e3);
    printf("Kalan kalinlik    : %.4f mm\n", (L_domain-recession_total)*1e3);

    // =============================================================================
    // Yüzey enerji dengesi (W/m^2)
    // Not: burada /1e3 YOK, doğrudan W/m^2 basıyoruz (kafa karışmasın)
    // =============================================================================
    double eps_f  = eps_surf(T_wall,alpha_eff[0]);
    double q_conv = h_eff*(T_recovery - T_wall);
    double q_rad  = eps_f*sigma_SB*(pow(T_surr,4) - pow(T_wall,4));
    double q_chem = (h_eff/cp_g_T(T_wall))*lookup_Tchem(Bg_now,L_prev);
    double q_cond = k_surf*(T_wall - T_old[1])/(x[1]-x[0]);
    double resid  = q_conv + q_rad + q_chem - q_cond;

    printf("\nYuzey enerji dengesi — Ewing tutarli (kW/m2):\n");
    printf("  q_conv = %+12.3f kW/m2\n", q_conv/1000.0);
    printf("  q_rad  = %+12.3f kW/m2\n", q_rad/1000.0);
    printf("  q_chem = %+12.3f kW/m2\n", q_chem/1000.0);
    printf("  q_cond = %+12.3f kW/m2\n", q_cond/1000.0);
    printf("  resid  = %+12.6f kW/m2  (~0 olmali)\n", resid/1000.0);
    printf("%s\n", string(80, '=').c_str());


    


    return 0;
}