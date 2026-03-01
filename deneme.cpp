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
// YARDIMCI
// =============================================================================
double inverseAverage(double v1, double d1, double v2, double d2)
{
    return ((v1/d1)+(v2/d2)) / ((1.0/d1)+(1.0/d2));
}

// =============================================================================
// BLOWING FACTOR
// =============================================================================
void blowing_factor(double m_dot, double rho_e, double u_e,
                    double h0, double cp_g, double& h_eff_out)
{
    double C_h0  = h0 / (rho_e * u_e * cp_g);
    double denom = rho_e * u_e * C_h0;
    bf_phi   = 2.0 * 0.5 * m_dot / denom;
    bf_Omega = bf_phi / (exp(bf_phi) - 1.0);
    h_eff_out = bf_Omega * h0;
}

// =============================================================================
// PATANKAR TDMA
// =============================================================================
vector<double> thomas_patankar(const vector<double>& a,
                               const vector<double>& b,
                               const vector<double>& c,
                               const vector<double>& d, const int N)
{
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
// EWING TCHEM TABLOSU
//
// Gerçek kullanımda: bPrimeTable.csv'den okuyun (1 atm, hava, unity Le)
// Burada placeholder değerler var — fiziksel büyüklük sırası doğru,
// ama doğrulama için gerçek tablo değerlerini kullanın.
//
// Tablo yapısı (Ewing Table 3):
//   Her Bg değeri için N_ROWS satır: (Tw, Bc, Tchem)
//   Tchem [J/kg], unity Le, Eq.(87): Tchem = -(1-B')*hw + Bc*hc + Bg*hg
// =============================================================================

const int N_BG   = 5;   // kaç farklı Bg değeri var
const int N_ROWS = 10;  // her Bg seksiyonunda kaç T satırı var

// B'g değerleri (sıralı)
double tbl_Bg[N_BG] = {0.00, 0.05, 0.10, 0.20, 0.40};

// Her (Bg, satır) için duvar sıcaklığı [K]
double tbl_Tw[N_BG][N_ROWS] = {
    {300,  600,  900, 1200, 1500, 1800, 2100, 2400, 2700, 3000},
    {300,  600,  900, 1200, 1500, 1800, 2100, 2400, 2700, 3000},
    {300,  600,  900, 1200, 1500, 1800, 2100, 2400, 2700, 3000},
    {300,  600,  900, 1200, 1500, 1800, 2100, 2400, 2700, 3000},
    {300,  600,  900, 1200, 1500, 1800, 2100, 2400, 2700, 3000},
};

// Her (Bg, satır) için B'c (char consumption B-prime)
double tbl_Bc[N_BG][N_ROWS] = {
    {0.087, 0.088, 0.090, 0.110, 0.160, 0.175, 0.176, 0.180, 0.190, 0.210},
    {0.060, 0.065, 0.070, 0.095, 0.145, 0.158, 0.162, 0.168, 0.178, 0.198},
    {0.020, 0.030, 0.045, 0.075, 0.130, 0.145, 0.150, 0.158, 0.168, 0.188},
    {0.010, 0.015, 0.025, 0.055, 0.110, 0.128, 0.135, 0.145, 0.156, 0.176},
    {0.005, 0.008, 0.015, 0.035, 0.085, 0.105, 0.115, 0.128, 0.142, 0.163},
};

// Her (Bg, satır) için Tchem [J/kg]
double tbl_Tchem[N_BG][N_ROWS] = {
    {-2.7e5, -5.4e5, -7.2e5, -8.8e5, -7.7e5, -1.02e6, -1.27e6, -1.54e6, -1.83e6, -2.14e6},
    {-2.7e5, -5.4e5, -7.3e5, -9.0e5, -8.0e5, -1.05e6, -1.30e6, -1.57e6, -1.86e6, -2.17e6},
    {-2.7e5, -5.5e5, -7.4e5, -9.2e5, -8.3e5, -1.09e6, -1.34e6, -1.61e6, -1.90e6, -2.22e6},
    {-2.8e5, -5.6e5, -7.6e5, -9.5e5, -8.8e5, -1.14e6, -1.40e6, -1.68e6, -1.98e6, -2.30e6},
    {-2.9e5, -5.8e5, -7.9e5, -9.9e5, -9.5e5, -1.20e6, -1.47e6, -1.76e6, -2.07e6, -2.40e6},
};

// =============================================================================
// TCHEM TABLO LOOKUP FONKSİYONLARI
// =============================================================================

// Bg için alt-üst bracket bul
void find_bg_bracket(double Bg, int& ilo, int& ihi, double& frac)
{
    if (Bg <= tbl_Bg[0])       { ilo=0;      ihi=0;      frac=0.0; return; }
    if (Bg >= tbl_Bg[N_BG-1])  { ilo=N_BG-1; ihi=N_BG-1; frac=0.0; return; }
    for (int i = 0; i < N_BG-1; i++)
        if (tbl_Bg[i+1] >= Bg) {
            ilo = i; ihi = i+1;
            frac = (Bg - tbl_Bg[i]) / (tbl_Bg[i+1] - tbl_Bg[i]);
            return;
        }
}

// Verilen 2D tablodan (Bg, L) interpolasyon
double lookup_val(double arr[][N_ROWS], double Bg, double L)
{
    int ilo, ihi; double frac;
    find_bg_bracket(Bg, ilo, ihi, frac);

    // L: 1..N_ROWS arası gerçek sayı → 0-tabanlı index
    double idx = L - 1.0;
    idx = max(0.0, min((double)(N_ROWS-1), idx));
    int jlo = (int)idx;
    int jhi = min(jlo+1, N_ROWS-1);
    double t = idx - jlo;

    double v_lo = arr[ilo][jlo] + t*(arr[ilo][jhi] - arr[ilo][jlo]);
    double v_hi = arr[ihi][jlo] + t*(arr[ihi][jhi] - arr[ihi][jlo]);
    return v_lo + frac*(v_hi - v_lo);
}

double lookup_Tw(double Bg, double L)    { return lookup_val(tbl_Tw,    Bg, L); }
double lookup_Bc(double Bg, double L)    { return lookup_val(tbl_Bc,    Bg, L); }
double lookup_Tchem(double Bg, double L) { return lookup_val(tbl_Tchem, Bg, L); }

// =============================================================================
// MAIN
// =============================================================================
int main()
{
    cout << string(80, '=') << "\n";
    cout << "ABLASYON ISI TRANSFERI - EWING TCHEM TABLOSU\n";
    cout << string(80, '=') << "\n";

    // =========================================================================
    // PARAMETRELER
    // =========================================================================
    const double t_end = 500.0;
    const double dt    = 0.1;
    const int    nstep = (int)round(t_end / dt);

    const int    N_nodes  = 200;
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

    const double k_v   = 0.17,   k_c   = 3.0;
    const double cp_v  = 1200.0, cp_c  = 1800.0;
    const double Q_p   = 0.5e7;

    const double phi_v      = 0.10,  phi_c      = 0.80;
    const double Gamma_perm = 1e-11, mu_g       = 3e-5;
    const double cp_g       = 1000.0;

    const double T_back = 300.0, P_surf = 101325.0, P_back = 101325.0;

    const double h_0_external = 200.0;
    const double T_recovery   = 8000.0;
    const double rho_e        = 1.2,  u_e = 1500.0;

    const double sigma_SB     = 5.67e-8;
    const double T_surr       = 300.0, emissivity = 0.85;

    // iterasyon parametreleri
    const int    max_iter  = 40;
    const double eps_T     = 0.1;
    const double eps_sdot  = 1e-9;

    printf("\n[INIT] dt=%.4f s, t_end=%.1f s, adim=%d\n", dt, t_end, nstep);
    printf("  rho_virgin=%.1f  rho_char=%.1f kg/m3\n", rho_virgin, rho_char_total);
    printf("  dx=%.3f mm, N=%d\n", (x[1]-x[0])*1e3, N_nodes);

    // =========================================================================
    // BASLANGIC KOSULLARI
    // =========================================================================
    vector<double>         T_old(N_nodes, T_back);
    vector<double>         T_prev(N_nodes, T_back);
    vector<double>         P_old(N_nodes, P_back);
    vector<double>         rho_solid_old(N_nodes, rho_virgin);
    vector<vector<double>> alpha_old(N_comp, vector<double>(N_nodes, 0.0));

    double recession_total = 0.0;
    double T_wall          = T_back;
    double m_dot_g         = 0.0;
    double k_surf;
    double h_eff           = h_0_external;
    double m_dot_surface   = 0.0;
    double sdot            = 0.0;

    // Ewing: NR başlangıç tahmini (line number)
    double L_prev = 1.0;

    // =========================================================================
    // CSV
    // =========================================================================
    ofstream fout("ablation_history.csv");
    fout << "time,Twall,T1,P0,mdot,heff,sdot_mm_s,recession_mm,thickness_mm,"
            "Bg,Bc,L_nr,iter\n";
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
                pyr_I  = B_arr[c] * exp(-E_over_R[c] / T_old[i]) * dt;
                pyr_ai = alpha_old[c][i];
                if (Psi_arr[c] == 1.0)
                    pyr_anew = 1.0 - (1.0 - pyr_ai) * exp(-pyr_I);
                else {
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

            k_node[i] = k_v*(1.0-alpha_eff[i]) + k_c*alpha_eff[i];
        }

        // =====================================================================
        // 2. dT/dt
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
            dT_dt[i] = (T_old[i] - T_prev[i]) / dt;

        // =====================================================================
        // 3. BASINC
        // =====================================================================
        for (int i = 0; i < N_nodes; i++)
            a_P[i]=b_P[i]=c_P[i]=d_P[i]=0.0;

        a_P[0] = 1.0;         d_P[0]         = P_surf;
        a_P[N_nodes-1] = 1.0; d_P[N_nodes-1] = P_back;

        for (int i = 1; i < N_nodes-1; i++)
        {
            p_dxl = x[i]   - x[i-1];
            p_dxr = x[i+1] - x[i];
            p_dxi = 0.5*(p_dxl + p_dxr);

            p_phi   = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];
            p_atime = p_phi / (R_univ*T_old[i]) * p_dxi / dt;
            p_K     = P_old[i] / (R_univ*T_old[i]) * Gamma_perm/mu_g;
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

        P_new = thomas_patankar(a_P, b_P, c_P, d_P, N_nodes);

        // =====================================================================
        // 4. BLOWING (ilk tahmin)
        // =====================================================================
        m_dot_surface = (P_new[0]/(R_univ*T_old[0]))
                      * (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf;
        blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);
        m_dot_g = m_dot_surface;

        // =====================================================================
        // 5. ITERATIF YUZEY + SICAKLIK COZUMU
        // =====================================================================
        k_surf = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);

        double T1        = T_old[1];
        double sdot_iter = sdot;
        int    iter_count = 0;

        // Ewing: diagnostik
        double Bg_now = 0.0, Bc_now = 0.0, L_cur = L_prev;

        for (int it = 0; it < max_iter; it++)
        {
            iter_count = it + 1;

            // --- (i) blowing güncelle ---
            m_dot_surface = (P_new[0]/(R_univ*T_wall))
                          * (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf;
            blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);
            m_dot_g = m_dot_surface;

            // --- a) EWING TCHEM RECONCILIATION (Eq.83, 86, 90) ---

            // Half-pass katsayıları: q_cond = hp_A*Tw + hp_B
            double hp_A = k_surf / dx_surf;
            double hp_B = -(k_surf / dx_surf) * T1;

            // B'g = m_dot_g / (rho_e*u_e*C_M),  C_M = h_0/( rho_e*u_e*cp_g )
            Bg_now = m_dot_g * cp_g / h_0_external;
            Bg_now = max(0.0, Bg_now);

            // NR on L (Ewing Eq.90):
            // f(L) = hp_A*Tw(L) + eps*sigma*Tw(L)^4 + hp_B
            //       - h_eff*T_recovery - (h_eff/cp_g)*Tchem(L)
            //       - eps*sigma*T_surr^4  =  0
            L_cur = L_prev;
            for (int nr = 0; nr < 50; nr++)
            {
                double Tw_L  = lookup_Tw(Bg_now, L_cur);
                double Tc_L  = lookup_Tchem(Bg_now, L_cur);

                double f = hp_A*Tw_L + emissivity*sigma_SB*pow(Tw_L,4) + hp_B
                         - h_eff*T_recovery - (h_eff/cp_g)*Tc_L
                         - emissivity*sigma_SB*pow(T_surr,4);

                // numerik türev
                double dL_fd = 0.01;
                double Tw_L2 = lookup_Tw(Bg_now, L_cur + dL_fd);
                double Tc_L2 = lookup_Tchem(Bg_now, L_cur + dL_fd);
                double f2    = hp_A*Tw_L2 + emissivity*sigma_SB*pow(Tw_L2,4) + hp_B
                             - h_eff*T_recovery - (h_eff/cp_g)*Tc_L2
                             - emissivity*sigma_SB*pow(T_surr,4);

                double dfdL = (f2 - f) / dL_fd;
                if (fabs(dfdL) < 1e-30) break;

                double step = -f / dfdL;
                step  = max(-10.0, min(10.0, step));   // salınım önle
                L_cur = max(1.0, min((double)N_ROWS, L_cur + step));

                if (fabs(step) < 1e-6) break;
            }

            // Tablodan sonuçları al
            double T_wall_new = lookup_Tw(Bg_now, L_cur);
            Bc_now = lookup_Bc(Bg_now, L_cur);

            // sdot = B'c * rho_e*u_e*C_M / rho_char   (Ewing Eq.47)
            // rho_e*u_e*C_M = h_0_external / cp_g
            double sdot_new = Bc_now * (h_0_external / cp_g) / rho_char_total;
            sdot = max(0.0, sdot_new);

            L_prev = L_cur;

            // --- (ii) sdot convergence ---
            double sdot_old_iter = sdot_iter;
            sdot_iter = sdot;

            // --- (iii) blowing'i yeni Twall ile tekrar bağla ---
            m_dot_surface = (P_new[0]/(R_univ*T_wall_new))
                          * (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf;
            blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);
            m_dot_g = m_dot_surface;

            // --- b) SICAKLIK TDMA ---
            for (int i = 0; i < N_nodes; i++)
                a_T[i]=b_T[i]=c_T[i]=d_T[i]=0.0;

            a_T[0] = 1.0; d_T[0] = T_wall_new;

            for (int i = 1; i < N_nodes-1; i++)
            {
                t_dxl = x[i]   - x[i-1];
                t_dxr = x[i+1] - x[i];
                t_dxi = 0.5*(t_dxl + t_dxr);

                t_cp  = cp_v*(1-alpha_eff[i]) + cp_c*alpha_eff[i];
                t_phi = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];

                t_ke = 2.0*k_node[i]*k_node[i+1] / (k_node[i]+k_node[i+1]);
                t_kw = 2.0*k_node[i]*k_node[i-1] / (k_node[i]+k_node[i-1]);

                t_rhogas = P_new[i] / (R_univ*T_old[i]);
                t_rhoc   = rho_solid_new[i]*t_cp + t_rhogas*t_phi*cp_g;
                t_atime  = t_rhoc * t_dxi / dt;
                t_ae     = t_ke / t_dxr;
                t_aw     = t_kw / t_dxl;

                t_mdotgas = t_rhogas * (Gamma_perm/mu_g)
                          * (P_new[i+1]-P_new[i-1]) / (t_dxl+t_dxr);
                t_hgrad   = cp_g * (T_old[i+1]-T_old[i-1]) / (t_dxl+t_dxr);

                t_hsolid = t_cp * T_old[i];
                t_hpyro  = t_hsolid + rho_solid_new[i]*(cp_v-cp_c)*T_old[i]
                                     / (rho_virgin-rho_char_total);
                t_Spyro  = -(Q_p - t_hpyro + cp_g*T_old[i]) * drho_dt[i];
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
                t_cp     = cp_v*(1-alpha_eff[j]) + cp_c*alpha_eff[j];
                t_phi    = phi_v*(1-alpha_eff[j]) + phi_c*alpha_eff[j];
                t_rhogas = P_new[j] / (R_univ*T_old[j]);
                t_rhoc   = rho_solid_new[j]*t_cp + t_rhogas*t_phi*cp_g;
                t_atime  = t_rhoc * (0.5*t_dxl) / dt;
                t_kw     = 2.0*k_node[j]*k_node[j-1] / (k_node[j]+k_node[j-1]);
                a_T[j]   = t_atime + t_kw/t_dxl;
                b_T[j]   = 0.0;
                c_T[j]   = t_kw/t_dxl;
                d_T[j]   = t_atime * T_old[j];
            }

            T_new = thomas_patankar(a_T, b_T, c_T, d_T, N_nodes);

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
            for (int j = 0; j < N_nodes; j++)
            {
                rho_solid_new[j] =
                    Gamma * ( (rho_v_comp[0]-(rho_v_comp[0]-rho_c_comp[0])*alpha_new[0][j])
                            +(rho_v_comp[1]-(rho_v_comp[1]-rho_c_comp[1])*alpha_new[1][j]) )
                + (1-Gamma)*( rho_v_comp[2]-(rho_v_comp[2]-rho_c_comp[2])*alpha_new[2][j] );
                alpha_eff[j] = (rho_virgin - rho_solid_new[j]) / (rho_virgin - rho_char_total);
                k_node[j]    = k_v*(1.0-alpha_eff[j]) + k_c*alpha_eff[j];
            }
            k_surf = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);
        }

        // =====================================================================
        // 7. UPDATE
        // =====================================================================
        T_prev = T_old;

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
                 << iter_count << "\n";
        }

        if ((n % (save_every*10)) == 0 || n == 1 || n == nstep)
        {
            printf("t=%7.1fs | Tw=%7.1fK | T1=%7.1fK | heff=%6.1f | "
                   "Bg=%.3f Bc=%.3f L=%.2f | "
                   "sdot=%8.4fmm/s | erim=%8.4fmm | kalan=%7.3fmm | it=%d\n",
                   time, T_wall, T_old[1], h_eff,
                   Bg_now, Bc_now, L_prev,
                   sdot*1e3, recession_total*1e3, thickness_now*1e3,
                   iter_count);
        }
    }

    fout.close();

    printf("\n%s\n", string(80, '=').c_str());
    printf("TAMAMLANDI\n");
    printf("Final Twall       : %.2f K\n", T_wall);
    printf("Final T[1]        : %.2f K\n", T_old[1]);
    printf("Final h_eff       : %.2f W/m2K\n", h_eff);
    printf("Final s_dot       : %.4f mm/s\n", sdot*1e3);
    printf("Toplam erime      : %.4f mm\n", recession_total*1e3);
    printf("Kalan kalinlik    : %.4f mm\n", (L_domain-recession_total)*1e3);

    double q_conv = h_eff*(T_recovery - T_wall);
    double q_rad  = emissivity*sigma_SB*(pow(T_surr,4) - pow(T_wall,4));
    double q_gas  = m_dot_g*(cp_g*T_old[1] - cp_g*T_wall);
    double q_cond = k_surf*(T_wall - T_old[1])/(x[1]-x[0]);
    double resid  = q_conv + q_rad + q_gas - q_cond;

    printf("\nYuzey isi akisi dengesi (kW/m^2):\n");
    printf("  q_conv = %+12.3f kW/m2\n", q_conv/1000.0);
    printf("  q_rad  = %+12.3f kW/m2\n", q_rad/1000.0);
    printf("  q_gas  = %+12.3f kW/m2\n", q_gas/1000.0);
    printf("  q_cond = %+12.3f kW/m2\n", q_cond/1000.0);
    printf("  resid  = %+12.3f kW/m2\n", resid/1000.0);
    printf("  Convergence Rate = %.6f %%\n", (resid/q_conv)*100.0);
    printf("%s\n", string(80, '=').c_str());

    return 0;
}