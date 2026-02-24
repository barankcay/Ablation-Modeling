import numpy as np
from scipy.interpolate import interp1d

print("="*80)
print("ABLASYON ISI TRANSFERİ v9 - MOVING BOUNDARY")
print("Yüzey gerilemesi: s_dot → x[0] kayması → grid güncelleme")
print("="*80)

# =============================================================================
# ZAMAN PARAMETRELERİ
# =============================================================================
t_end = 200.0
dt    = 0.01
nstep = int(np.round(t_end / dt))
print(f"\n[INIT] dt={dt} s, t_end={t_end} s, adım={nstep}")

# =============================================================================
# FİZİKSEL PARAMETRELER
# =============================================================================
N_nodes  = 100
L_domain = 0.05          # başlangıç kalınlığı [m]
x        = np.linspace(0, L_domain, N_nodes)   # başlangıç grid'i

# Piroliz
N_comp     = 3
B          = np.array([1.40e4,  9.75e8, 0.0])
Psi        = np.array([3.0,     3.0,    0.0])
E          = np.array([71.14e6, 169.98e6, 0.0])   # J/kmol
R_univ     = 8.314
E_over_R   = (E / 1000.0) / R_univ               # K

rho_v_comp = np.array([325.015,  973.926, 2066.380])
rho_c_comp = np.array([0.0,      518.998, 2066.380])
Gamma      = 0.422

rho_virgin = Gamma*(rho_v_comp[0]+rho_v_comp[1]) + (1-Gamma)*rho_v_comp[2]
rho_char   = Gamma*(rho_c_comp[0]+rho_c_comp[1]) + (1-Gamma)*rho_c_comp[2]

# Termal
k_v  = 0.17;   k_c  = 3
cp_v = 1200.0; cp_c = 1800.0
Q_p  = 0.5e7

# Gaz / gözenek
phi_v      = 0.10; phi_c = 0.80
Gamma_perm = 1e-11
mu_g       = 3e-5
cp_g       = 1000.0

# Sınır koşulları
T_back = 300.0
P_surf = 101325.0
P_back = 101325.0

# Dış sınır tabakası
h_0_external = 200.0
T_recovery   = 8000.0
rho_e        = 1.2
u_e          = 1500.0

# Radyasyon + ablasyon
sigma_SB     = 5.67e-8
T_surr       = 300.0
emissivity   = 0.85
T_ablation   = 1996.0
Delta_H_melt = 160000.0

print(f"  rho_virgin={rho_virgin:.1f}  rho_char={rho_char:.1f} kg/m³")
print(f"  Başlangıç dx={np.diff(x)[0]*1e3:.3f} mm, N={N_nodes}")

# =============================================================================
# BAŞLANGIÇ KOŞULLARI
# =============================================================================
T_old         = np.full(N_nodes, T_back)
T_prev        = T_old.copy()
alpha_old     = np.zeros((N_comp, N_nodes))
P_old         = np.full(N_nodes, P_back)
rho_solid_old = np.full(N_nodes, rho_virgin)

# =============================================================================
# GERİLEME TAKİBİ
# =============================================================================
recession_total = 0.0        # toplam eriyen kalınlık [m]
x_surface       = 0.0        # yüzey konumu (başlangıçta 0)

# =============================================================================
# YARDIMCI FONKSİYONLAR
# =============================================================================

def blowing_factor(m_dot, rho_e, u_e, h0, cp_g, lam=0.5):
    C_h0  = h0 / (rho_e * u_e * cp_g)
    denom = rho_e * u_e * C_h0
    if abs(denom) < 1e-30:
        return 1.0, h0
    phi = 2.0 * lam * m_dot / denom
    if abs(phi) < 1e-4:
        Omega = 1.0 - phi/2.0 + phi**2/12.0
    elif phi > 50.0:
        Omega = phi * np.exp(-phi)
    else:
        ep    = np.exp(phi) - 1.0
        Omega = phi / ep if abs(ep) > 1e-30 else 1.0
    Omega = np.clip(Omega, 0.0, 1.0)
    return Omega, Omega * h0


def thomas_patankar(a, b, c, d):
    N = len(a)
    P = np.zeros(N); Q = np.zeros(N)
    P[0] = b[0] / a[0]
    Q[0] = d[0] / a[0]
    for i in range(1, N):
        denom = a[i] - c[i]*P[i-1]
        if abs(denom) < 1e-30:
            denom = 1e-30
        P[i] = b[i] / denom
        Q[i] = (d[i] + c[i]*Q[i-1]) / denom
    T = np.zeros(N)
    T[-1] = Q[-1]
    for i in range(N-2, -1, -1):
        T[i] = P[i]*T[i+1] + Q[i]
    return T


def solve_Twall_NR(h_eff, m_dot_g, k_surf, T1, Tw_guess,
                   T_recovery, emissivity, sigma_SB, T_surr, dx_surf, cp_g):
    delta_n = dx_surf
    h_g_in  = cp_g * T1
    A = (h_eff*T_recovery + (k_surf/delta_n)*T1
         + emissivity*sigma_SB*T_surr**4 + m_dot_g*h_g_in)
    B = -h_eff - k_surf/delta_n - m_dot_g*cp_g
    C = -emissivity*sigma_SB
    Tw = float(Tw_guess)
    for _ in range(100):
        f  = A + B*Tw + C*Tw**4
        df = B + 4.0*C*Tw**3
        if abs(df) < 1e-20:
            break
        dTw = np.clip(-f/df, -500.0, 500.0)
        Tw  = max(Tw + dTw, 100.0)
        if abs(dTw) < 1e-6:
            break
    return Tw


def compute_sdot(h_eff, m_dot_g, k_surf, T1, T_ablation,
                 emissivity, sigma_SB, T_surr, dx_surf, cp_g,
                 rho_char, Delta_H_melt, T_recovery):
    Tw      = T_ablation
    delta_n = dx_surf
    numerator = (h_eff*(T_recovery - Tw)
                 + m_dot_g*(cp_g*T1 - cp_g*Tw)
                 + emissivity*sigma_SB*(T_surr**4 - Tw**4)
                 - k_surf*(Tw - T1)/delta_n)
    denominator = rho_char * Delta_H_melt
    if abs(denominator) < 1e-12:
        return 0.0
    return max(numerator / denominator, 0.0)


def remap_to_new_grid(x_old, x_new, fields):
    """
    Eski grid x_old üzerindeki alanları yeni grid x_new'e interpolasyon ile taşı.
    fields: dict of {name: array}  — 1D veya (N_comp, N_nodes) shape
    Dönüş: aynı yapıda yeni dict
    """
    result = {}
    for name, arr in fields.items():
        if arr.ndim == 1:
            f   = interp1d(x_old, arr, kind='linear',
                           bounds_error=False,
                           fill_value=(arr[0], arr[-1]))
            result[name] = f(x_new)
        else:
            # (N_comp, N_nodes)
            new_arr = np.zeros((arr.shape[0], len(x_new)))
            for comp in range(arr.shape[0]):
                f = interp1d(x_old, arr[comp], kind='linear',
                             bounds_error=False,
                             fill_value=(arr[comp, 0], arr[comp, -1]))
                new_arr[comp] = f(x_new)
            result[name] = new_arr
    return result

# =============================================================================
# DURUM
# =============================================================================
ablation_active = False
T_wall = T_old[0]
s_dot  = 0.0

# History
time_hist      = []; Twall_hist    = []; T0_hist    = []
P0_hist        = []; mdot_hist     = []; heff_hist  = []
sdot_hist      = []; mode_hist     = []
recession_hist = []; thickness_hist = []
save_every = 10

print("\n[LOOP] Time loop başlıyor...\n")
print(f"{'t':>8s} | {'MODE':4s} | {'Tw':>7s}K | {'T1':>7s}K | "
      f"{'mdot':>10s} | {'h_eff':>6s} | "
      f"{'sdot':>10s} | {'erim':>8s} | {'kalan':>8s}")
print("-"*100)

for n in range(1, nstep + 1):
    time = n * dt

    # Grid spacing — yüzey geriledikçe x[0] kayar ama uniform grid koruyoruz
    # x dizisi her adımda güncel: x[0] = x_surface (yüzey konumu)
    dx_arr   = np.diff(x)                     # her aralık farklı olabilir
    dx_surf  = x[1] - x[0]                    # yüzey-node1 mesafesi

    # =========================================================================
    # 1. PİROLİZ
    # =========================================================================
    alpha_new = np.zeros((N_comp, N_nodes))
    for comp in range(N_comp):
        if B[comp] == 0.0:
            alpha_new[comp, :] = alpha_old[comp, :]
            continue
        I_arr = B[comp] * np.exp(-E_over_R[comp] / np.maximum(T_old, 1.0)) * dt
        if Psi[comp] == 1.0:
            alpha_new[comp, :] = 1.0 - (1.0 - alpha_old[comp, :]) * np.exp(-I_arr)
        else:
            exp_ = 1.0 / (1.0 - Psi[comp])
            bkt  = I_arr*(Psi[comp]-1.0) + (1.0-alpha_old[comp, :])**(1.0-Psi[comp])
            alpha_new[comp, :] = 1.0 - np.maximum(bkt, 1e-30)**exp_
        alpha_new[comp, :] = np.clip(alpha_new[comp, :], 0.0, 1.0)

    rho_comp = np.zeros((N_comp, N_nodes))
    for comp in range(N_comp):
        rho_comp[comp] = (rho_v_comp[comp]
                          - (rho_v_comp[comp]-rho_c_comp[comp])*alpha_new[comp])
    rho_solid_new = Gamma*(rho_comp[0]+rho_comp[1]) + (1-Gamma)*rho_comp[2]
    alpha_eff     = np.clip((rho_virgin-rho_solid_new)/(rho_virgin-rho_char), 0.0, 1.0)

    dalpha_dt = (alpha_new - alpha_old) / dt
    drho_dt   = np.zeros(N_nodes)
    for comp in range(N_comp):
        drho_dt += (rho_v_comp[comp]-rho_c_comp[comp]) * dalpha_dt[comp]

    k_node = k_v*(1.0-alpha_eff) + k_c*alpha_eff

    # =========================================================================
    # 2. dT/dt — BACKWARD DIFFERENCE
    # =========================================================================
    dT_dt = (T_old - T_prev) / dt

    # =========================================================================
    # 3. BASINÇ DENKLEMİ
    # =========================================================================
    a_P = np.zeros(N_nodes); b_P = np.zeros(N_nodes)
    c_P = np.zeros(N_nodes); d_P = np.zeros(N_nodes)
    a_P[0]  = 1.0; d_P[0]  = P_surf
    a_P[-1] = 1.0; d_P[-1] = P_back

    for i in range(1, N_nodes-1):
        dxi    = x[i+1] - x[i-1]   # 2*dx eğer uniform
        dxl    = x[i]   - x[i-1]
        dxr    = x[i+1] - x[i]
        dx_i   = 0.5 * dxi          # kontrol hacmi genişliği

        phi_i  = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i]
        C_cap  = phi_i / (R_univ * T_old[i])
        a_time = C_cap * dx_i / dt

        K      = P_old[i] / (R_univ * T_old[i]) * Gamma_perm / mu_g
        a_east = K / dxr
        a_west = K / dxl

        S_const  = drho_dt[i]
        S_P_temp = phi_i / (R_univ * T_old[i]**2) * dT_dt[i]
        alpha_eff_old_i  = np.clip(
            (rho_virgin-rho_solid_old[i])/(rho_virgin-rho_char), 0.0, 1.0)
        dalpha_eff_dt_i  = (alpha_eff[i] - alpha_eff_old_i) / dt
        S_P_poro = -(phi_v-phi_c) / (R_univ*T_old[i]) * dalpha_eff_dt_i
        S_linear = S_P_temp + S_P_poro

        a_P[i] = a_time + a_east + a_west - S_linear*dx_i
        b_P[i] = a_east
        c_P[i] = a_west
        d_P[i] = S_const*dx_i + a_time*P_old[i]

    P_new = np.maximum(thomas_patankar(a_P, b_P, c_P, d_P), 1.0)

    # =========================================================================
    # 4. BLOWING
    # =========================================================================
    v_D_surf      = (Gamma_perm/mu_g) * (P_new[1]-P_new[0]) / dx_surf
    rho_g_surf    = P_new[0] / (R_univ * max(T_old[0], 1.0))
    m_dot_surface = rho_g_surf * v_D_surf
    _, h_eff      = blowing_factor(m_dot_surface, rho_e, u_e, h_0_external, cp_g)

    # =========================================================================
    # 5. YÜZEY ENERJİ DENGESİ — STATE MACHINE
    # =========================================================================
    m_dot_g = max(0.0, m_dot_surface)
    k_surf  = 2.0*k_node[0]*k_node[1] / (k_node[0]+k_node[1])
    T1      = T_old[1]

    if not ablation_active:
        Tw_stat = solve_Twall_NR(
            h_eff, m_dot_g, k_surf, T1, Tw_guess=T_wall,
            T_recovery=T_recovery, emissivity=emissivity,
            sigma_SB=sigma_SB, T_surr=T_surr, dx_surf=dx_surf, cp_g=cp_g)
        if not np.isfinite(Tw_stat) or Tw_stat < 100.0:
            Tw_stat = T_wall

        if Tw_stat >= T_ablation:
            sdot_try = compute_sdot(
                h_eff, m_dot_g, k_surf, T1, T_ablation,
                emissivity, sigma_SB, T_surr, dx_surf, cp_g,
                rho_char, Delta_H_melt, T_recovery)
            ablation_active = sdot_try > 0.0
            T_wall = T_ablation
            s_dot  = sdot_try if ablation_active else 0.0
        else:
            T_wall = Tw_stat
            s_dot  = 0.0
    else:
        sdot_try = compute_sdot(
            h_eff, m_dot_g, k_surf, T1, T_ablation,
            emissivity, sigma_SB, T_surr, dx_surf, cp_g,
            rho_char, Delta_H_melt, T_recovery)
        if sdot_try > 0.0:
            T_wall = T_ablation
            s_dot  = sdot_try
        else:
            ablation_active = False
            Tw_stat = solve_Twall_NR(
                h_eff, m_dot_g, k_surf, T1, Tw_guess=T_ablation-1.0,
                T_recovery=T_recovery, emissivity=emissivity,
                sigma_SB=sigma_SB, T_surr=T_surr, dx_surf=dx_surf, cp_g=cp_g)
            T_wall = min(Tw_stat if np.isfinite(Tw_stat) else T_ablation-1.0,
                         T_ablation)
            s_dot = 0.0

    # =========================================================================
    # 6. SICAKLIK DENKLEMİ
    # =========================================================================
    a_T = np.zeros(N_nodes); b_T = np.zeros(N_nodes)
    c_T = np.zeros(N_nodes); d_T = np.zeros(N_nodes)
    a_T[0] = 1.0; d_T[0] = T_wall

    for i in range(1, N_nodes-1):
        dxl = x[i]   - x[i-1]
        dxr = x[i+1] - x[i]
        dx_i = 0.5*(dxl + dxr)

        cp_local = cp_v*(1.0-alpha_eff[i]) + cp_c*alpha_eff[i]
        phi_i    = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i]

        k_e = 2.0*k_node[i]*k_node[i+1] / (k_node[i]+k_node[i+1])
        k_w = 2.0*k_node[i]*k_node[i-1] / (k_node[i]+k_node[i-1])

        rho_gas = P_new[i] / (R_univ * max(T_old[i], 1.0))
        rho_c_i = rho_solid_new[i]*cp_local + rho_gas*phi_i*cp_g
        a_time  = rho_c_i * dx_i / dt
        a_east  = k_e / dxr
        a_west  = k_w / dxl

        m_dot_gas  = (P_new[i]/(R_univ*max(T_old[i],1.0))) * \
                     (Gamma_perm/mu_g) * (P_new[i+1]-P_new[i-1])/(dxl+dxr)
        h_gas_grad = cp_g * (T_old[i+1]-T_old[i-1]) / (dxl+dxr)

        h_solid    = cp_local * T_old[i]
        h_eff_pyro = (h_solid
                      + rho_solid_new[i]*(cp_v-cp_c)*T_old[i]
                        / (rho_virgin-rho_char))
        S_pyro  = -(Q_p - h_eff_pyro + cp_g*T_old[i]) * drho_dt[i] 
        #cp_g*Told gas entalpisini veriyo. sıcaklığa bağlı verildiği hali de var ona bi bak
        S_total = S_pyro - m_dot_gas*h_gas_grad

        a_T[i] = a_time + a_east + a_west
        b_T[i] = a_east
        c_T[i] = a_west
        d_T[i] = S_total*dx_i + a_time*T_old[i]

    # Arka yüzey: adiabatik
    i           = N_nodes - 1
    dxl_last    = x[i] - x[i-1]
    cp_last     = cp_v*(1.0-alpha_eff[i]) + cp_c*alpha_eff[i]
    phi_last    = phi_v*(1.0-alpha_eff[i]) + phi_c*alpha_eff[i]
    rho_g_last  = P_new[i] / (R_univ*max(T_old[i],1.0))
    rho_c_last  = rho_solid_new[i]*cp_last + rho_g_last*phi_last*cp_g
    a_time_last = rho_c_last * (0.5*dxl_last) / dt
    k_w_last    = 2.0*k_node[i]*k_node[i-1] / (k_node[i]+k_node[i-1])
    a_T[-1]     = a_time_last + k_w_last/dxl_last
    b_T[-1]     = 0.0
    c_T[-1]     = k_w_last/dxl_last
    d_T[-1]     = a_time_last * T_old[-1]

    T_new = thomas_patankar(a_T, b_T, c_T, d_T)
    if not np.all(np.isfinite(T_new)) or np.any(T_new < 0.0):
        print(f"  ⚠ t={time:.4f}: T_new fiziksel değil")
        T_new = T_old.copy()

    # =========================================================================
    # 7. YÜZEY GERİLEMESİ — MOVING BOUNDARY
    # =========================================================================
    recession_step   = s_dot * dt           # bu adımda eriyen mesafe [m]
    recession_total += recession_step
    x_surface       += recession_step       # yüzey konumu ilerliyor

    if recession_step > 0.0:
        # Yeni grid: yüzey x_surface'ten başlıyor, arka yüzey sabit kalıyor
        # → etkin kalınlık azalıyor
        thickness_current = L_domain - recession_total
        if thickness_current < 2*dx_surf:
            print(f"\n⚠ t={time:.2f}s: Malzeme tamamen eridi! "
                  f"Kalan={thickness_current*1e3:.3f} mm")
            break

        x_new = np.linspace(x_surface, L_domain, N_nodes)

        # Mevcut alanları eski grid'den yeni grid'e taşı
        fields_old = {
            'T':     T_new,
            'P':     P_new,
            'alpha': alpha_new,
            'rho_s': rho_solid_new,
            'T_prev': T_prev,
        }
        fields_new = remap_to_new_grid(x, x_new, fields_old)

        # Grid ve alanları güncelle
        x             = x_new
        T_new         = fields_new['T']
        P_new         = fields_new['P']
        alpha_new     = fields_new['alpha']
        rho_solid_new = fields_new['rho_s']
        T_prev        = fields_new['T_prev']

        # Yüzey BC yeniden uygula (interpolasyon bozabilir)
        T_new[0]  = T_wall
        P_new[0]  = P_surf

    # =========================================================================
    # 8. UPDATE
    # =========================================================================
    if recession_step == 0.0:
        T_prev = T_old.copy()
    # (recession varsa T_prev zaten remap'te güncellendi)

    T_old         = T_new.copy()
    P_old         = P_new.copy()
    alpha_old     = alpha_new.copy()
    rho_solid_old = rho_solid_new.copy()

    # =========================================================================
    # HISTORY + PRINT
    # =========================================================================
    thickness_now = L_domain - recession_total

    if (n % save_every) == 0 or n == 1 or n == nstep:
        time_hist.append(time)
        Twall_hist.append(T_wall)
        T0_hist.append(T_new[0])
        P0_hist.append(P_new[0])
        mdot_hist.append(m_dot_surface)
        heff_hist.append(h_eff)
        sdot_hist.append(s_dot)
        mode_hist.append("ABL" if ablation_active else "STA")
        recession_hist.append(recession_total * 1e3)     # mm
        thickness_hist.append(thickness_now * 1e3)       # mm

    if (n % (save_every*20) == 0) or n == 1 or n == nstep:
        mode = "ABL" if ablation_active else "STA"
        print(f"t={time:7.1f}s | {mode} | "
              f"Tw={T_wall:7.1f}K | T1={T_old[1]:7.1f}K | "
              f"mdot={m_dot_surface:+.2e} | "
              f"sdot={s_dot*1e6:8.3f}um/s | "
              f"erim={recession_total*1e3:7.3f}mm | "
              f"kalan={thickness_now*1e3:7.3f}mm")

# =============================================================================
# BİTİŞ
# =============================================================================
print("\n" + "="*80)
print("✅ TAMAMLANDI")
print(f"Final Twall          : {T_wall:.2f} K")
print(f"Final T[0]           : {T_old[0]:.2f} K")
print(f"Final T[1]           : {T_old[1]:.2f} K")
print(f"Final P[0]           : {P_old[0]/1000:.3f} kPa")
print(f"Final mdot           : {m_dot_surface:.3e} kg/m²s")
print(f"Final h_eff          : {h_eff:.2f} W/m²K")
print(f"Final s_dot          : {s_dot*1e6:.4f} μm/s")
print(f"Ablation active      : {ablation_active}")
print(f"Toplam erime         : {recession_total*1e3:.4f} mm")
print(f"Kalan kalınlık       : {(L_domain-recession_total)*1e3:.4f} mm")
print(f"Başlangıç kalınlığı  : {L_domain*1e3:.1f} mm")

q_conv = h_eff*(T_recovery - T_wall)
q_rad  = emissivity*sigma_SB*(T_surr**4 - T_wall**4)
q_gas  = m_dot_g*(cp_g*T_old[1] - cp_g*T_wall)
q_cond = k_surf*(T_wall - T_old[1]) / (x[1]-x[0])
resid  = q_conv + q_rad + q_gas - q_cond

print("\nYüzey ısı akısı dengesi:")
print(f"  q_conv = {q_conv/1e3:+9.3f} kW/m²")
print(f"  q_rad  = {q_rad /1e3:+9.3f} kW/m²")
print(f"  q_gas  = {q_gas /1e3:+9.3f} kW/m²")
print(f"  q_cond = {q_cond/1e3:+9.3f} kW/m² (içe)")
print(f"  resid  = {resid /1e3:+9.3f} kW/m²")
print("="*80)
print("\nHistory dizileri:")
print("  time_hist, Twall_hist, T0_hist, P0_hist, mdot_hist,")
print("  heff_hist, sdot_hist, mode_hist, recession_hist, thickness_hist")