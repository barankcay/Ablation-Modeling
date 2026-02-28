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

// solve_Twall_NR
double nr_A, nr_B, nr_C;
double nr_Tw, nr_Tw3, nr_f, nr_df, nr_dTw;

// compute_sdot
double sd_Tw4, sd_num, sd_den;

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
// NEWTON-RAPHSON: yüzey sicakligi
// f(Tw) = A + B*Tw + C*Tw^4 = 0
// =============================================================================
double solve_Twall_NR(double h_eff, double m_dot_g, double k_surf,
                      double T1, double Tw_guess,
                      double T_recovery, double emissivity,
                      double sigma_SB, double T_surr,
                      double dx_surf, double cp_g)
{
    nr_A = h_eff*T_recovery + (k_surf/dx_surf)*T1
         + emissivity*sigma_SB*pow(T_surr,4)
         + m_dot_g*cp_g*T1;

    nr_B = -h_eff - k_surf/dx_surf - m_dot_g*cp_g;
    nr_C = -emissivity*sigma_SB;

    nr_Tw = Tw_guess;

    for (int iter = 0; iter < 100; iter++)
    {
        nr_Tw3 = nr_Tw*nr_Tw*nr_Tw;
        nr_f   = nr_A + nr_B*nr_Tw + nr_C*nr_Tw3*nr_Tw;
        nr_df  = nr_B + 4.0*nr_C*nr_Tw3;

        if (fabs(nr_df) < 1e-20) break;

        nr_dTw = -nr_f / nr_df;
        nr_Tw  = nr_Tw + nr_dTw;

        if (nr_Tw < 1.0) nr_Tw = 1.0;   // güvenlik
        if (fabs(nr_dTw) < 1e-6) break;
    }
    return nr_Tw;
}

// =============================================================================
// ABLASYON HIZI
// =============================================================================
double compute_sdot(double h_eff, double m_dot_g, double k_surf,
                    double T1, double T_ablation,
                    double emissivity, double sigma_SB, double T_surr,
                    double dx_surf, double cp_g,
                    double rho_char, double Delta_H_melt,
                    double T_recovery)
{
    sd_Tw4 = pow(T_ablation,4);

    sd_num = h_eff*(T_recovery - T_ablation)
           + m_dot_g*cp_g*(T1 - T_ablation)
           + emissivity*sigma_SB*(pow(T_surr,4) - sd_Tw4)
           - k_surf*(T_ablation - T1)/dx_surf;

    sd_den = rho_char * Delta_H_melt;

    double sdot_raw = sd_num / sd_den;

    const double sdot_min = 1e-7;  // m/s  (0.1 µm/s)
    return (sdot_raw > sdot_min) ? sdot_raw : 0.0;
}

// =============================================================================
// MAIN
// =============================================================================
int main()
{
    cout << string(80, '=') << "\n";
    cout << "ABLASYON ISI TRANSFERI v12.2 (C++) - ITERATIF YUZEY (STABILIZE PATCH)\n";
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
    const double T_surr       = 300.0, emissivity   = 0.85;
    const double T_ablation   = 1996.0, Delta_H_melt = 160000.0;

    // iterasyon parametreleri
    const int    max_iter  = 40;
    const double eps_T     = 0.1;   // K
    const double eps_sdot  = 1e-9;   // m/s

    // ---- STABILIZE PATCH PARAMS ----
    const double dT_hys  = 3.0;   // K (ABL aç/kapa histerezisi)
    const double ur_sdot = 0.2;   // sdot under-relax (0..1)
    // --------------------------------

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

    bool   ablation_active = false;
    double T_wall        = T_back;
    double s_dot         = 0.0;
    double m_dot_g       = 0.0;
    double k_surf        = k_v;
    double h_eff         = h_0_external;
    double m_dot_surface = 0.0;

    // =========================================================================
    // CSV
    // =========================================================================
    ofstream fout("ablation_history.csv");
    // debug kolonları eklendi: Tw_stat, sdot_raw, dP01, dx_surf, k_surf
    fout << "time,Twall,T1,P0,mdot,heff,sdot_mm_s,recession_mm,thickness_mm,mode,iter,"
            "Tw_stat,sdot_raw_mm_s,dP01_Pa,dx_surf_m,k_surf\n";
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

    // debug cache
    double dbg_Tw_stat = T_back;
    double dbg_sdot_raw = 0.0;

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

            k_node[i] = k_v*(1.0-alpha_eff[i]) + k_c*alpha_eff[i];
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

        P_new = thomas_patankar(a_P, b_P, c_P, d_P);

        // =====================================================================
        // 4. BLOWING (dış iterasyon dışında ilk tahmin)
        // =====================================================================
        m_dot_surface = (P_new[0]/(R_univ*T_old[0]))
                      * (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf;
        blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);
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
        double sdot_iter  = s_dot;
        int    iter_count = 0;

        for (int it = 0; it < max_iter; it++)
        {
            iter_count = it + 1;

            // --- (i) iterasyon içi: mdot/heff güncelle (Twall’a bağlı yüzey yoğunluğu) ---
            // önce Twall tahmini olarak mevcut T_wall kullan
            double Tw_for_rho = T_wall;
            m_dot_surface = (P_new[0]/(R_univ*Tw_for_rho))
                          * (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf;
            blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);
            m_dot_g = m_dot_surface;

            // --- a) STATE MACHINE ---
            double T_wall_new = T_wall;
            double sdot_new_raw = 0.0;

            // statik Tw (ablasyon yok varsayımı)
            dbg_Tw_stat = solve_Twall_NR(
                h_eff, m_dot_g, k_surf, T1, T_wall,
                T_recovery, emissivity, sigma_SB, T_surr, dx_surf, cp_g);

            if (!ablation_active)
            {
                // ABL açma şartı: Tablation + histerezis
                if (dbg_Tw_stat >= (T_ablation + dT_hys))
                {
                    sdot_new_raw = compute_sdot(
                        h_eff, m_dot_g, k_surf, T1, T_ablation,
                        emissivity, sigma_SB, T_surr, dx_surf, cp_g,
                        rho_char_total, Delta_H_melt, T_recovery);

                    if (sdot_new_raw > 0.0)
                    {
                        ablation_active = true;
                        T_wall_new = T_ablation;
                    }
                    else
                    {
                        T_wall_new = dbg_Tw_stat;
                        sdot_new_raw = 0.0;
                    }
                }
                else
                {
                    T_wall_new = dbg_Tw_stat;
                    sdot_new_raw = 0.0;
                }
            }
            else
            {
                // ABL modunda sdot
                sdot_new_raw = compute_sdot(
                    h_eff, m_dot_g, k_surf, T1, T_ablation,
                    emissivity, sigma_SB, T_surr, dx_surf, cp_g,
                    rho_char_total, Delta_H_melt, T_recovery);

                if (sdot_new_raw > 0.0)
                {
                    T_wall_new = T_ablation;
                }
                else
                {
                    // Kapatma şartı: Tw_stat <= Tablation - histerezis
                    if (dbg_Tw_stat <= (T_ablation - dT_hys))
                    {
                        ablation_active = false;
                        T_wall_new = dbg_Tw_stat;
                        sdot_new_raw = 0.0;
                    }
                    else
                    {
                        // histerezis bandında: ABL'yi koru, sdot=0 (chatter kır)
                        T_wall_new = T_ablation;
                        sdot_new_raw = 0.0;
                    }
                }
            }

            dbg_sdot_raw = sdot_new_raw;

            // --- (ii) sdot under-relax ---
            double sdot_old_iter = sdot_iter;
            sdot_iter = (1.0 - ur_sdot) * sdot_iter + ur_sdot * sdot_new_raw;
            if (sdot_iter < 0.0) sdot_iter = 0.0;

            // --- (iii) mdot/heff’i yeni Twall ile tekrar bağla (çok pahalı değil) ---
            Tw_for_rho = T_wall_new;
            m_dot_surface = (P_new[0]/(R_univ*Tw_for_rho))
                          * (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf;
            blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g, h_eff);
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

        s_dot = sdot_iter;

        // =====================================================================
        // 6. MOVING BOUNDARY + REMAP
        // Patch:
        //  - remap distance fabs+clamp (inverseAverage zaten bunu yapıyor)
        // =====================================================================
        double recession_step = s_dot * dt;
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
                T_prev[j]        = inverseAverage(T_prev[j],        d1, T_prev[j+1],        d2);

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
                k_node[j]    = k_v*(1.0-alpha_eff[j]) + k_c*alpha_eff[j];
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
                 << s_dot*1e3 << "," << recession_total*1e3 << ","
                 << thickness_now*1e3 << ","
                 << (ablation_active ? "ABL" : "STA") << ","
                 << iter_count << ","
                 << dbg_Tw_stat << "," << dbg_sdot_raw*1e3 << ","
                 << (P_new[1]-P_new[0]) << "," << dx_surf << "," << k_surf
                 << "\n";
        }

        if ((n % (save_every*10)) == 0 || n == 1 || n == nstep)
        {
            printf("t=%7.1fs | %s | Tw=%7.1fK | T1=%7.1fK | heff=%6.1f | "
                   "sdot=%8.3fmm/s | erim=%8.4fmm | kalan=%7.3fmm | it=%d\n",
                   time, ablation_active ? "ABL" : "STA",
                   T_wall, T_old[1], h_eff,
                   s_dot*1e3, recession_total*1e3, thickness_now*1e3,
                   iter_count);
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
    printf("Final s_dot       : %.4f mm/s\n", s_dot*1e3);
    printf("Ablation active   : %s\n", ablation_active ? "EVET" : "HAYIR");
    printf("Toplam erime      : %.4f mm\n", recession_total*1e3);
    printf("Kalan kalinlik    : %.4f mm\n", (L_domain-recession_total)*1e3);

    // =============================================================================
    // Yüzey enerji dengesi (W/m^2)
    // Not: burada /1e3 YOK, doğrudan W/m^2 basıyoruz (kafa karışmasın)
    // =============================================================================
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
    printf("%s\n", string(80, '=').c_str());


    


    return 0;
}