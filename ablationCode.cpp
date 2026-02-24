#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

// =============================================================================
// YARDIMCI: linspace
// =============================================================================
vector<double> linspace(double start, double end, int n)
{
    vector<double> v(n);
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++)
        v[i] = start + i * step;
    return v;
}

// =============================================================================
// YARDIMCI: 1D linear interpolasyon (tablo → değer)
// bounds_error yok: sınır dışı → uç değer döner
// =============================================================================
double interp1d(const vector<double>& xp, const vector<double>& fp, double x)
{
    int n = (int)xp.size();
    if (x <= xp[0])  return fp[0];
    if (x >= xp[n-1]) return fp[n-1];
    // binary search
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (xp[mid] <= x) lo = mid; else hi = mid;
    }
    double t = (x - xp[lo]) / (xp[hi] - xp[lo]);
    return fp[lo] + t * (fp[hi] - fp[lo]);
}

// =============================================================================
// YARDIMCI: grid'i yeniden oluştur ve alanları taşı (remap)
// x_old → x_new, 1D vektörü interpolasyonla doldurur
// =============================================================================
vector<double> remap1d(const vector<double>& x_old,
                        const vector<double>& x_new,
                        const vector<double>& f_old)
{
    int n = (int)x_new.size();
    vector<double> f_new(n);
    for (int i = 0; i < n; i++)
        f_new[i] = interp1d(x_old, f_old, x_new[i]);
    return f_new;
}

// 2D: (N_comp x N_nodes) için her satırı ayrı remap
vector<vector<double>> remap2d(const vector<double>& x_old,
                                const vector<double>& x_new,
                                const vector<vector<double>>& f_old)
{
    int ncomp = (int)f_old.size();
    int nnew  = (int)x_new.size();
    vector<vector<double>> f_new(ncomp, vector<double>(nnew));
    for (int c = 0; c < ncomp; c++)
        f_new[c] = remap1d(x_old, x_new, f_old[c]);
    return f_new;
}

// =============================================================================
// BLOWING FACTOR
// Omega = phi / (exp(phi) - 1)
// phi ~ 0 ise Taylor serisi kullan
// =============================================================================
double blowing_factor(double m_dot, double rho_e, double u_e,
                      double h0, double cp_g, double& h_eff_out)
{
    double C_h0  = h0 / (rho_e * u_e * cp_g);
    double denom = rho_e * u_e * C_h0;
    if (fabs(denom) < 1e-30) { h_eff_out = h0; return 1.0; }

    double phi = 2.0 * 0.5 * m_dot / denom;
    double Omega;
    if (fabs(phi) < 1e-4)
        Omega = 1.0 - phi/2.0 + phi*phi/12.0;
    else if (phi > 50.0)
        Omega = phi * exp(-phi);
    else {
        double ep = exp(phi) - 1.0;
        Omega = (fabs(ep) > 1e-30) ? phi / ep : 1.0;
    }
    Omega = max(0.0, min(1.0, Omega));
    h_eff_out = Omega * h0;
    return Omega;
}

// =============================================================================
// PATANKAR TDMA
// a[i]*T[i] = b[i]*T[i+1] + c[i]*T[i-1] + d[i]
// Patankar (1980) Ch.3 Eqs 3.16-3.19
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
    for (int i = 1; i < N; i++) {
        double den = a[i] - c[i] * P[i-1];
        if (fabs(den) < 1e-30) den = 1e-30;
        P[i] = b[i] / den;
        Q[i] = (d[i] + c[i] * Q[i-1]) / den;
    }
    T[N-1] = Q[N-1];
    for (int i = N-2; i >= 0; i--)
        T[i] = P[i] * T[i+1] + Q[i];
    return T;
}

// =============================================================================
// NEWTON-RAPHSON: yüzey sıcaklığı (stationary)
// f(Tw) = A + B*Tw + C*Tw^4 = 0
// Alanyalioglu (2019) Eq. 3.33-3.35
// =============================================================================
double solve_Twall_NR(double h_eff, double m_dot_g, double k_surf,
                      double T1, double Tw_guess,
                      double T_recovery, double emissivity,
                      double sigma_SB, double T_surr,
                      double dx_surf, double cp_g)
{
    double delta_n = dx_surf;
    double h_g_in  = cp_g * T1;

    double A = h_eff*T_recovery + (k_surf/delta_n)*T1
             + emissivity*sigma_SB*T_surr*T_surr*T_surr*T_surr
             + m_dot_g*h_g_in;
    double B = -h_eff - k_surf/delta_n - m_dot_g*cp_g;
    double C = -emissivity*sigma_SB;

    double Tw = Tw_guess;
    for (int iter = 0; iter < 100; iter++) {
        double Tw3 = Tw*Tw*Tw;
        double f   = A + B*Tw + C*Tw3*Tw;
        double df  = B + 4.0*C*Tw3;
        if (fabs(df) < 1e-20) break;
        double dTw = -f / df;
        // adım sınırı
        dTw = max(-500.0, min(500.0, dTw));
        Tw += dTw;
        Tw  = max(Tw, 100.0);
        if (fabs(dTw) < 1e-6) break;
    }
    return Tw;
}

// =============================================================================
// ABLASYON HIZI (s_dot)
// Alanyalioglu (2019) Eq. 3.37
// denominator = rho_char * Delta_H_melt  (sadece latent heat)
// =============================================================================
double compute_sdot(double h_eff, double m_dot_g, double k_surf,
                    double T1, double T_ablation,
                    double emissivity, double sigma_SB, double T_surr,
                    double dx_surf, double cp_g,
                    double rho_char, double Delta_H_melt,
                    double T_recovery)
{
    double Tw      = T_ablation;
    double delta_n = dx_surf;
    double Tw4     = Tw*Tw*Tw*Tw;

    double numerator = h_eff*(T_recovery - Tw)
                     + m_dot_g*(cp_g*T1 - cp_g*Tw)
                     + emissivity*sigma_SB*(T_surr*T_surr*T_surr*T_surr - Tw4)
                     - k_surf*(Tw - T1)/delta_n;

    double denominator = rho_char * Delta_H_melt;
    if (fabs(denominator) < 1e-12) return 0.0;
    return max(numerator / denominator, 0.0);
}

// =============================================================================
// MAIN
// =============================================================================
int main()
{
    cout << string(80, '=') << endl;
    cout << "ABLASYON ISI TRANSFERI v9 (C++) - MOVING BOUNDARY" << endl;
    cout << string(80, '=') << endl;

    // =========================================================================
    // ZAMAN PARAMETRELERİ
    // =========================================================================
    const double t_end = 200.0;
    const double dt    = 0.01;
    const int    nstep = (int)round(t_end / dt);
    printf("\n[INIT] dt=%.4f s, t_end=%.1f s, adim=%d\n", dt, t_end, nstep);

    // =========================================================================
    // FİZİKSEL PARAMETRELER
    // =========================================================================
    const int    N_nodes  = 100;
    const double L_domain = 0.05;   // m
    vector<double> x = linspace(0.0, L_domain, N_nodes);

    // --- Piroliz ---
    const int N_comp = 3;
    vector<double> B_arr   = {1.40e4,  9.75e8, 0.0};
    vector<double> Psi_arr = {3.0,     3.0,    0.0};
    vector<double> E_arr   = {71.14e6, 169.98e6, 0.0};  // J/kmol
    const double R_univ    = 8.314;
    vector<double> E_over_R(N_comp);
    for (int c = 0; c < N_comp; c++)
        E_over_R[c] = (E_arr[c] / 1000.0) / R_univ;   // K

    vector<double> rho_v_comp = {325.015,  973.926, 2066.380};
    vector<double> rho_c_comp = {0.0,      518.998, 2066.380};
    const double Gamma = 0.422;

    double rho_virgin = Gamma*(rho_v_comp[0]+rho_v_comp[1])
                      + (1-Gamma)*rho_v_comp[2];
    double rho_char_total = Gamma*(rho_c_comp[0]+rho_c_comp[1])
                          + (1-Gamma)*rho_c_comp[2];

    // --- Termal ---
    const double k_v   = 0.17;
    const double k_c   = 3.0;
    const double cp_v  = 1200.0;
    const double cp_c  = 1800.0;
    const double Q_p   = 0.5e7;

    // --- Gaz / gözenek ---
    const double phi_v      = 0.10;
    const double phi_c      = 0.80;
    const double Gamma_perm = 1e-11;
    const double mu_g       = 3e-5;
    const double cp_g       = 1000.0;

    // --- Sınır koşulları ---
    const double T_back = 300.0;
    const double P_surf = 101325.0;
    const double P_back = 101325.0;

    // --- Dış sınır tabakası ---
    const double h_0_external = 200.0;
    const double T_recovery   = 8000.0;
    const double rho_e        = 1.2;
    const double u_e          = 1500.0;

    // --- Radyasyon + ablasyon ---
    const double sigma_SB     = 5.67e-8;
    const double T_surr       = 300.0;
    const double emissivity   = 0.85;
    const double T_ablation   = 1996.0;
    const double Delta_H_melt = 160000.0;

    printf("  rho_virgin=%.1f  rho_char=%.1f kg/m3\n",
           rho_virgin, rho_char_total);
    printf("  Baslangic dx=%.3f mm, N=%d\n",
           (x[1]-x[0])*1e3, N_nodes);

    // =========================================================================
    // BAŞLANGIÇ KOŞULLARI
    // =========================================================================
    vector<double> T_old(N_nodes, T_back);
    vector<double> T_prev(N_nodes, T_back);
    vector<double> P_old(N_nodes, P_back);
    vector<double> rho_solid_old(N_nodes, rho_virgin);
    vector<vector<double>> alpha_old(N_comp, vector<double>(N_nodes, 0.0));

    // =========================================================================
    // GERİLEME TAKİBİ
    // =========================================================================
    double recession_total = 0.0;
    double x_surface       = 0.0;

    // =========================================================================
    // DURUM
    // =========================================================================
    bool   ablation_active = false;
    double T_wall = T_back;
    double s_dot  = 0.0;
    double m_dot_g = 0.0;
    double k_surf  = k_v;
    double h_eff   = h_0_external;
    double m_dot_surface = 0.0;

    // =========================================================================
    // HISTORY (CSV için)
    // =========================================================================
    ofstream fout("ablation_history.csv");
    fout << "time,Twall,T1,P0,mdot,heff,sdot_um_s,"
            "recession_mm,thickness_mm,mode\n";

    const int save_every = 2;

    printf("\n[LOOP] Time loop basliyor...\n\n");
    printf("%8s | %4s | %8s | %8s | %10s | %10s | %8s | %8s\n",
           "t(s)", "MODE", "Tw(K)", "T1(K)", "mdot", "sdot(um/s)",
           "erim(mm)", "kalan(mm)");
    printf("%s\n", string(80, '-').c_str());

    // =========================================================================
    // TIME LOOP
    // =========================================================================
    for (int n = 1; n <= nstep; n++)
    {
        double time = n * dt;

        double dx_surf = x[1] - x[0];

        // =====================================================================
        // 1. PİROLİZ
        // =====================================================================
        vector<vector<double>> alpha_new(N_comp, vector<double>(N_nodes, 0.0));

        for (int c = 0; c < N_comp; c++) {
            if (B_arr[c] == 0.0) {
                alpha_new[c] = alpha_old[c];
                continue;
            }
            for (int i = 0; i < N_nodes; i++) {
                double T_i = max(T_old[i], 1.0);
                double I   = B_arr[c] * exp(-E_over_R[c] / T_i) * dt;
                double a_i = alpha_old[c][i];
                double a_new;
                if (Psi_arr[c] == 1.0) {
                    a_new = 1.0 - (1.0 - a_i) * exp(-I);
                } else {
                    double exp_ = 1.0 / (1.0 - Psi_arr[c]);
                    double bkt  = I*(Psi_arr[c]-1.0)
                                + pow(max(1.0-a_i, 1e-30), 1.0-Psi_arr[c]);
                    bkt   = max(bkt, 1e-30);
                    a_new = 1.0 - pow(bkt, exp_);
                }
                alpha_new[c][i] = max(0.0, min(1.0, a_new));
            }
        }

        // rho solid
        vector<double> rho_solid_new(N_nodes);
        for (int i = 0; i < N_nodes; i++) {
            double r = 0.0;
            for (int c = 0; c < N_comp; c++)
                r += (c < 2 ? Gamma : (1-Gamma))
                   * (rho_v_comp[c] - (rho_v_comp[c]-rho_c_comp[c])*alpha_new[c][i]);
            // Düzeltme: Gamma sadece c=0,1 için, (1-Gamma) c=2 için
            rho_solid_new[i] = Gamma*(
                (rho_v_comp[0]-(rho_v_comp[0]-rho_c_comp[0])*alpha_new[0][i])
               +(rho_v_comp[1]-(rho_v_comp[1]-rho_c_comp[1])*alpha_new[1][i])
            ) + (1-Gamma)*(
                (rho_v_comp[2]-(rho_v_comp[2]-rho_c_comp[2])*alpha_new[2][i])
            );
        }

        // alpha_eff, drho_dt, k_node
        vector<double> alpha_eff(N_nodes), drho_dt(N_nodes, 0.0), k_node(N_nodes);
        for (int i = 0; i < N_nodes; i++) {
            alpha_eff[i] = max(0.0, min(1.0,
                (rho_virgin - rho_solid_new[i]) / (rho_virgin - rho_char_total)));
            for (int c = 0; c < N_comp; c++)
                drho_dt[i] += (rho_v_comp[c]-rho_c_comp[c])
                            * (alpha_new[c][i]-alpha_old[c][i]) / dt;
            k_node[i] = k_v*(1.0-alpha_eff[i]) + k_c*alpha_eff[i];
        }

        // =====================================================================
        // 2. dT/dt — BACKWARD DIFFERENCE
        // =====================================================================
        vector<double> dT_dt(N_nodes);
        for (int i = 0; i < N_nodes; i++)
            dT_dt[i] = (T_old[i] - T_prev[i]) / dt;

        // =====================================================================
        // 3. BASINÇ DENKLEMİ
        // =====================================================================
        vector<double> a_P(N_nodes,0), b_P(N_nodes,0),
                       c_P(N_nodes,0), d_P(N_nodes,0);
        a_P[0] = 1.0; d_P[0] = P_surf;
        a_P[N_nodes-1] = 1.0; d_P[N_nodes-1] = P_back;

        for (int i = 1; i < N_nodes-1; i++) {
            double dxl  = x[i]   - x[i-1];
            double dxr  = x[i+1] - x[i];
            double dx_i = 0.5*(dxl+dxr);

            double phi_i  = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];
            double C_cap  = phi_i / (R_univ * T_old[i]);
            double a_time = C_cap * dx_i / dt;

            double K      = P_old[i]/(R_univ*T_old[i]) * Gamma_perm/mu_g;
            double a_east = K / dxr;
            double a_west = K / dxl;

            double alpha_eff_old_i = max(0.0, min(1.0,
                (rho_virgin-rho_solid_old[i])/(rho_virgin-rho_char_total)));
            double dalpha_eff_dt_i = (alpha_eff[i]-alpha_eff_old_i) / dt;

            double S_P_temp = phi_i/(R_univ*T_old[i]*T_old[i]) * dT_dt[i];
            double S_P_poro = -(phi_v-phi_c)/(R_univ*T_old[i]) * dalpha_eff_dt_i;
            double S_linear = S_P_temp + S_P_poro;

            a_P[i] = a_time + a_east + a_west - S_linear*dx_i;
            b_P[i] = a_east;
            c_P[i] = a_west;
            d_P[i] = drho_dt[i]*dx_i + a_time*P_old[i];
        }

        vector<double> P_new = thomas_patankar(a_P, b_P, c_P, d_P);
        for (auto& p : P_new) p = max(p, 1.0);

        // =====================================================================
        // 4. BLOWING
        // =====================================================================
        double v_D_surf   = (Gamma_perm/mu_g)*(P_new[1]-P_new[0])/dx_surf;
        double rho_g_surf = P_new[0] / (R_univ * max(T_old[0], 1.0));
        m_dot_surface     = rho_g_surf * v_D_surf;
        blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);

        // =====================================================================
        // 5. YÜZEY ENERJİ DENGESİ — STATE MACHINE
        // =====================================================================
        m_dot_g = max(0.0, m_dot_surface);
        k_surf  = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1]);
        double T1 = T_old[1];

        if (!ablation_active) {
            double Tw_stat = solve_Twall_NR(
                h_eff, m_dot_g, k_surf, T1, T_wall,
                T_recovery, emissivity, sigma_SB, T_surr, dx_surf, cp_g);
            if (!isfinite(Tw_stat) || Tw_stat < 100.0) Tw_stat = T_wall;

            if (Tw_stat >= T_ablation) {
                double sdot_try = compute_sdot(
                    h_eff, m_dot_g, k_surf, T1, T_ablation,
                    emissivity, sigma_SB, T_surr, dx_surf, cp_g,
                    rho_char_total, Delta_H_melt, T_recovery);
                ablation_active = sdot_try > 0.0;
                T_wall = T_ablation;
                s_dot  = ablation_active ? sdot_try : 0.0;
            } else {
                T_wall = Tw_stat;
                s_dot  = 0.0;
            }
        } else {
            double sdot_try = compute_sdot(
                h_eff, m_dot_g, k_surf, T1, T_ablation,
                emissivity, sigma_SB, T_surr, dx_surf, cp_g,
                rho_char_total, Delta_H_melt, T_recovery);
            if (sdot_try > 0.0) {
                T_wall = T_ablation;
                s_dot  = sdot_try;
            } else {
                ablation_active = false;
                double Tw_stat = solve_Twall_NR(
                    h_eff, m_dot_g, k_surf, T1, T_ablation-1.0,
                    T_recovery, emissivity, sigma_SB, T_surr, dx_surf, cp_g);
                T_wall = min(isfinite(Tw_stat) ? Tw_stat : T_ablation-1.0,
                             T_ablation);
                s_dot = 0.0;
            }
        }

        // =====================================================================
        // 6. SICAKLIK DENKLEMİ
        // =====================================================================
        vector<double> a_T(N_nodes,0), b_T(N_nodes,0),
                       c_T(N_nodes,0), d_T(N_nodes,0);

        // Yüzey: Dirichlet
        a_T[0] = 1.0; d_T[0] = T_wall;

        for (int i = 1; i < N_nodes-1; i++) {
            double dxl  = x[i]   - x[i-1];
            double dxr  = x[i+1] - x[i];
            double dx_i = 0.5*(dxl+dxr);

            double cp_local = cp_v*(1-alpha_eff[i]) + cp_c*alpha_eff[i];
            double phi_i    = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];

            double k_e = 2.0*k_node[i]*k_node[i+1]/(k_node[i]+k_node[i+1]);
            double k_w = 2.0*k_node[i]*k_node[i-1]/(k_node[i]+k_node[i-1]);

            double rho_gas = P_new[i]/(R_univ*max(T_old[i],1.0));
            double rho_c_i = rho_solid_new[i]*cp_local + rho_gas*phi_i*cp_g;
            double a_time  = rho_c_i * dx_i / dt;
            double a_east  = k_e / dxr;
            double a_west  = k_w / dxl;

            double m_dot_gas  = (P_new[i]/(R_univ*max(T_old[i],1.0)))
                              * (Gamma_perm/mu_g)
                              * (P_new[i+1]-P_new[i-1])/(dxl+dxr);
            double h_gas_grad = cp_g*(T_old[i+1]-T_old[i-1])/(dxl+dxr);

            double h_solid    = cp_local*T_old[i];
            double h_eff_pyro = h_solid
                              + rho_solid_new[i]*(cp_v-cp_c)*T_old[i]
                                / (rho_virgin-rho_char_total);
            double S_pyro  = -(Q_p - h_eff_pyro + cp_g*T_old[i]) * drho_dt[i];
            double S_total = S_pyro - m_dot_gas*h_gas_grad;

            a_T[i] = a_time + a_east + a_west;
            b_T[i] = a_east;
            c_T[i] = a_west;
            d_T[i] = S_total*dx_i + a_time*T_old[i];
        }

        // Arka yüzey: adiabatik
        {
            int i = N_nodes - 1;
            double dxl_last = x[i] - x[i-1];
            double cp_last  = cp_v*(1-alpha_eff[i]) + cp_c*alpha_eff[i];
            double phi_last = phi_v*(1-alpha_eff[i]) + phi_c*alpha_eff[i];
            double rho_g_l  = P_new[i]/(R_univ*max(T_old[i],1.0));
            double rho_c_l  = rho_solid_new[i]*cp_last + rho_g_l*phi_last*cp_g;
            double a_time_l = rho_c_l*(0.5*dxl_last)/dt;
            double k_w_last = 2.0*k_node[i]*k_node[i-1]/(k_node[i]+k_node[i-1]);
            a_T[i] = a_time_l + k_w_last/dxl_last;
            b_T[i] = 0.0;
            c_T[i] = k_w_last/dxl_last;
            d_T[i] = a_time_l * T_old[i];
        }

        vector<double> T_new = thomas_patankar(a_T, b_T, c_T, d_T);

        // Fiziksel kontrol
        bool bad = false;
        for (int i = 0; i < N_nodes; i++)
            if (!isfinite(T_new[i]) || T_new[i] < 0.0) { bad = true; break; }
        if (bad) {
            printf("  WARNING t=%.4f: T_new fiziksel degil, onceki adimda kal\n", time);
            T_new = T_old;
        }

        // =====================================================================
        // 7. YÜZEY GERİLEMESİ — MOVING BOUNDARY
        // =====================================================================
        double recession_step = s_dot * dt;
        recession_total += recession_step;
        x_surface       += recession_step;

        if (recession_step > 0.0) {
            double thickness_current = L_domain - recession_total;
            if (thickness_current < 2.0*dx_surf) {
                printf("\nWARNING t=%.2f: Malzeme tamamen eridi! Kalan=%.3f mm\n",
                       time, thickness_current*1e3);
                break;
            }

            vector<double> x_new = linspace(x_surface, L_domain, N_nodes);

            // Remap
            T_new         = remap1d(x, x_new, T_new);
            P_new         = remap1d(x, x_new, P_new);
            rho_solid_new = remap1d(x, x_new, rho_solid_new);
            T_prev        = remap1d(x, x_new, T_prev);
            alpha_new     = remap2d(x, x_new, alpha_new);

            x = x_new;

            // Yüzey BC yeniden uygula
            T_new[0] = T_wall;
            P_new[0] = P_surf;
        }

        // =====================================================================
        // 8. UPDATE
        // =====================================================================
        if (recession_step == 0.0)
            T_prev = T_old;

        T_old         = T_new;
        P_old         = P_new;
        alpha_old     = alpha_new;
        rho_solid_old = rho_solid_new;

        // =====================================================================
        // HISTORY + PRINT
        // =====================================================================
        double thickness_now = L_domain - recession_total;

        if ((n % save_every) == 0 || n == 1 || n == nstep) {
            string mode_str = ablation_active ? "ABL" : "STA";
            fout << time << ","
                 << T_wall << ","
                 << T_old[1] << ","
                 << P_new[0] << ","
                 << m_dot_surface << ","
                 << h_eff << ","
                 << s_dot*1e6 << ","
                 << recession_total*1e3 << ","
                 << thickness_now*1e3 << ","
                 << mode_str << "\n";
        }

        if ((n % (save_every*200)) == 0 || n == 1 || n == nstep) {
            printf("t=%7.1fs | %s | Tw=%7.1fK | T1=%7.1fK | "
                   "sdot=%8.3fum/s | erim=%7.3fmm | kalan=%7.3fmm\n",
                   time,
                   ablation_active ? "ABL" : "STA",
                   T_wall, T_old[1],
                   s_dot*1e6,
                   recession_total*1e3, thickness_now*1e3);
        }
    }

    // =========================================================================
    // BİTİŞ
    // =========================================================================
    fout.close();
    printf("\n%s\n", string(80, '=').c_str());
    printf("TAMAMLANDI\n");
    printf("Final Twall       : %.2f K\n", T_wall);
    printf("Final T[0]        : %.2f K\n", T_old[0]);
    printf("Final T[1]        : %.2f K\n", T_old[1]);
    printf("Final P[0]        : %.3f kPa\n", P_old[0]/1000.0);
    printf("Final mdot        : %.3e kg/m2s\n", m_dot_surface);
    printf("Final h_eff       : %.2f W/m2K\n", h_eff);
    printf("Final s_dot       : %.4f um/s\n", s_dot*1e6);
    printf("Ablation active   : %s\n", ablation_active ? "EVET" : "HAYIR");
    printf("Toplam erime      : %.4f mm\n", recession_total*1e3);
    printf("Kalan kalinlik    : %.4f mm\n", (L_domain-recession_total)*1e3);

    double q_conv = h_eff*(T_recovery - T_wall);
    double q_rad  = emissivity*sigma_SB*(T_surr*T_surr*T_surr*T_surr
                  - T_wall*T_wall*T_wall*T_wall);
    double q_gas  = m_dot_g*(cp_g*T_old[1] - cp_g*T_wall);
    double q_cond = k_surf*(T_wall - T_old[1])/(x[1]-x[0]);
    double resid  = q_conv + q_rad + q_gas - q_cond;

    printf("\nYuzey isi akisi dengesi:\n");
    printf("  q_conv = %+9.3f kW/m2\n", q_conv/1e3);
    printf("  q_rad  = %+9.3f kW/m2\n", q_rad/1e3);
    printf("  q_gas  = %+9.3f kW/m2\n", q_gas/1e3);
    printf("  q_cond = %+9.3f kW/m2\n", q_cond/1e3);
    printf("  resid  = %+9.3f kW/m2\n", resid/1e3);
    printf("%s\n", string(80, '=').c_str());
    printf("Sonuclar: ablation_history.csv\n");

    return 0;
}