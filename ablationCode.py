import numpy as np

print("="*80)
print("ABLASYON ISI TRANSFERİ - TEK ZAMAN ADIMI")
print("3 Komponent Piroliz + Convective BC with Blowing")
print("Kaynak: Ewing, Laker, Walker (2013) - J. Thermophysics and Heat Transfer")
print("="*80)

# =============================================================================
# BÖLÜM 1: FİZİKSEL PARAMETRELERİN TANIMI
# =============================================================================
print("\n[1/8] Fiziksel parametreler tanımlanıyor...")

# --- UZAMSAL DİSKRETİZASYON ---
N_nodes = 100
L_domain = 0.05
dx = L_domain / (N_nodes - 1)
x = np.linspace(0, L_domain, N_nodes)

print(f"  ✓ Grid: {N_nodes} düğüm, Δy = {dx*1000:.3f} mm")

# --- ZAMANSAL DİSKRETİZASYON ---
dt = 0.01

print(f"  ✓ Zaman adımı: Δt = {dt} s")

# =============================================================================
# MALZEME ÖZELLİKLERİ - 3 KOMPONENTLİ PİROLİZ
# =============================================================================

N_comp = 3

# Arrhenius parametreleri
B = np.array([1.40e4, 4.48e9, 0.0])
Psi = np.array([3.0, 3.0, 0.0])
E_over_R = np.array([8555.6, 20444.4, 0.0])

# Yoğunluklar
rho_v_comp = np.array([229.0, 972.0, 160.0])
rho_c_comp = np.array([0.0, 792.0, 160.0])

# Hacim fraksiyonu
Gamma = 0.1

# Toplam yoğunluklar
rho_virgin = Gamma * (rho_v_comp[0] + rho_v_comp[1]) + (1 - Gamma) * rho_v_comp[2]
rho_char = Gamma * (rho_c_comp[0] + rho_c_comp[1]) + (1 - Gamma) * rho_c_comp[2]

print(f"\n  Komponent A: B={B[0]:.2e}, Ψ={Psi[0]}, E/R={E_over_R[0]:.1f}K")
print(f"  Komponent B: B={B[1]:.2e}, Ψ={Psi[1]}, E/R={E_over_R[1]:.1f}K")
print(f"  Komponent C: İnert")
print(f"  Toplam: ρ_virgin={rho_virgin:.1f}, ρ_char={rho_char:.1f} kg/m³")

# --- TERMOFİZİKSEL ÖZELLİKLER ---
k_v = 0.17
k_c = 0.60
cp_v = 1200.0
cp_c = 1800.0
Q_p = 500000.0

print(f"  Termal: k_v={k_v}, k_c={k_c}, Q_p={Q_p/1000:.0f} kJ/kg")

# --- GAZ AKIŞI ---
phi_v = 0.10
phi_c = 0.80
Gamma_perm = 1e-11
mu_g = 3e-5
cp_g = 1000.0
R_univ = 8.314

print(f"  Gözeneklilik: φ_v={phi_v}, φ_c={phi_c}")

# --- SINIR KOŞULLARI ---
# NOT: T_surf artık kullanılmıyor - convective BC var!
T_back = 300.0
P_surf = 101325.0
P_back = 101325.0

# --- EXTERNAL BOUNDARY LAYER PARAMETERS ---
# Bu değerler dış HTc korelasyon kodundan gelecek
h_0_external = 5000.0       # Heat transfer coefficient without blowing (W/m²·K)
T_recovery = 1500.0         # Recovery temperature (K)
rho_e = 0.5                 # Boundary layer edge density (kg/m³)
u_e = 1500.0                # Boundary layer edge velocity (m/s)

print(f"\n  External BC (from HTc correlation):")
print(f"    h_0 = {h_0_external:.0f} W/m²·K (without blowing)")
print(f"    T_recovery = {T_recovery:.1f} K")
print(f"    ρ_e = {rho_e:.2f} kg/m³, u_e = {u_e:.1f} m/s")


# =============================================================================
# BÖLÜM 2: BAŞLANGIÇ KOŞULLARI
# =============================================================================
print("\n[2/8] Başlangıç koşulları...")

# Sıcaklık profili (recovery temperature'dan içe üstel azalma)
T_old = T_back + (T_recovery - T_back) * np.exp(-x / 0.01)

# Her komponent için α
alpha_old = np.zeros((N_comp, N_nodes))

# Basınç
P_old = np.full(N_nodes, P_back)

# Yoğunluk
rho_solid_old = np.full(N_nodes, rho_virgin)

print(f"  ✓ T: yüzey={T_old[0]:.0f}K, arka={T_old[-1]:.0f}K")
print(f"  ✓ α=0 (virgin), P={P_old[0]/1000:.1f}kPa")


# =============================================================================
# BÖLÜM 3: PİROLİZ HESAPLAMA
# =============================================================================
print("\n[3/8] Piroliz hesaplanıyor...")

alpha_new = np.zeros((N_comp, N_nodes))

for comp in range(N_comp):
    if B[comp] == 0.0:
        continue
    
    for i in range(N_nodes):
        I = B[comp] * np.exp(-E_over_R[comp] / T_old[i]) * dt
        
        if Psi[comp] == 1.0:
            alpha_new[comp, i] = 1.0 - (1.0 - alpha_old[comp, i]) / np.exp(I)
        else:
            exponent = 1.0 / (1.0 - Psi[comp])
            bracket = I * (Psi[comp] - 1.0) + (1.0 - alpha_old[comp, i])**(1.0 - Psi[comp])
            alpha_new[comp, i] = 1.0 - bracket**exponent
        
        alpha_new[comp, i] = np.clip(alpha_new[comp, i], 0.0, 1.0)

# Yoğunluklar
rho_comp = np.zeros((N_comp, N_nodes))
for comp in range(N_comp):
    rho_comp[comp, :] = rho_v_comp[comp] - (rho_v_comp[comp] - rho_c_comp[comp]) * alpha_new[comp, :]

rho_solid_new = Gamma * (rho_comp[0, :] + rho_comp[1, :]) + (1 - Gamma) * rho_comp[2, :]

# Efektif α
alpha_eff = (rho_virgin - rho_solid_new) / (rho_virgin - rho_char)

# Piroliz hızı
dalpha_dt = (alpha_new - alpha_old) / dt

# Toplam gaz üretimi
drho_dt = np.zeros(N_nodes)
for comp in range(N_comp):
    drho_dt += (rho_v_comp[comp] - rho_c_comp[comp]) * dalpha_dt[comp, :]

print(f"  ✓ Komponent A: max α={np.max(alpha_new[0,:]):.4f}")
print(f"  ✓ Komponent B: max α={np.max(alpha_new[1,:]):.4f}")
print(f"  ✓ Efektif: max α={np.max(alpha_eff):.4f}")


# =============================================================================
# BÖLÜM 4: ∂T/∂t TAHMİNİ
# =============================================================================
print("\n[4/8] ∂T/∂t tahmini...")

dT_dt = np.zeros(N_nodes)

for i in range(1, N_nodes - 1):
    k_local = k_v * (1.0 - alpha_eff[i]) + k_c * alpha_eff[i]
    cp_local = cp_v * (1.0 - alpha_eff[i]) + cp_c * alpha_eff[i]
    alpha_thermal = k_local / (rho_solid_old[i] * cp_local)
    
    d2T_dx2 = (T_old[i+1] - 2*T_old[i] + T_old[i-1]) / (dx**2)
    dT_dt[i] = alpha_thermal * d2T_dx2

print(f"  ✓ Max |∂T/∂t|={np.max(np.abs(dT_dt)):.2f} K/s")


# =============================================================================
# BÖLÜM 5: BASINÇ DENKLEMİ
# =============================================================================
print("\n[5/8] Basınç çözülüyor...")

a_P = np.zeros(N_nodes)
b_P = np.zeros(N_nodes)
c_P = np.zeros(N_nodes)
d_P = np.zeros(N_nodes)

# Sınır koşulları
a_P[0] = 1.0
d_P[0] = P_surf

a_P[-1] = 1.0
d_P[-1] = P_back

# İç düğümler
for i in range(1, N_nodes - 1):
    phi = phi_v * (1.0 - alpha_eff[i]) + phi_c * alpha_eff[i]
    
    C_cap = phi / (R_univ * T_old[i])
    a_time = C_cap * dx / dt
    
    K = P_old[i] / (R_univ * T_old[i]) * Gamma_perm / mu_g
    a_east = K / dx
    a_west = K / dx
    
    S_const = drho_dt[i]
    
    S_P_temp = phi / (R_univ * T_old[i]**2) * dT_dt[i]
    
    dalpha_eff_dt = (alpha_eff[i] - (rho_virgin - rho_solid_old[i])/(rho_virgin - rho_char)) / dt
    S_P_poro = -(phi_v - phi_c) / (R_univ * T_old[i]) * dalpha_eff_dt
    
    S_linear = S_P_temp + S_P_poro
    
    a_P[i] = a_time + a_east + a_west - S_linear * dx
    b_P[i] = a_east
    c_P[i] = a_west
    d_P[i] = S_const * dx + a_time * P_old[i]

# Thomas algoritması
for i in range(1, N_nodes):
    factor = c_P[i] / a_P[i-1]
    a_P[i] -= factor * b_P[i-1]
    d_P[i] -= factor * d_P[i-1]

P_new = np.zeros(N_nodes)
P_new[-1] = d_P[-1] / a_P[-1]
for i in range(N_nodes - 2, -1, -1):
    P_new[i] = (d_P[i] - b_P[i] * P_new[i+1]) / a_P[i]

print(f"  ✓ P: yüzey={P_new[0]/1000:.2f}kPa, arka={P_new[-1]/1000:.2f}kPa")


# =============================================================================
# BÖLÜM 6: BLOWING CORRECTION FUNCTIONS
# =============================================================================

def calculate_phi_parameter(m_dot_surface, rho_e, u_e, C_h0, lambda_param=0.5):
    """
    EQUATION 2.22 - Blowing rate parameter
    φ = 2*λ * (ṁ" / (ρ_e * u_e * C_h,0))
    """
    phi = 2.0 * lambda_param * m_dot_surface / (rho_e * u_e * C_h0)
    return phi


def calculate_blowing_factor(phi):
    """
    EQUATION 2.20 - Blowing reduction factor
    Ω_blowing = φ / (e^φ - 1)
    """
    if abs(phi) < 1e-6:
        # Taylor expansion: φ/(e^φ-1) ≈ 1 - φ/2
        Omega = 1.0 - phi/2.0
    elif phi > 20.0:
        # Asymptotic: φ/e^φ → 0
        Omega = phi * np.exp(-phi)
    else:
        # Standard formula
        Omega = phi / (np.exp(phi) - 1.0)
    
    return Omega


# =============================================================================
# BÖLÜM 7: SICAKLIK DENKLEMİ (CONVECTIVE BC with BLOWING)
# =============================================================================
print("\n[6/8] Sıcaklık çözülüyor (convective BC with blowing)...")

a_T = np.zeros(N_nodes)
b_T = np.zeros(N_nodes)
c_T = np.zeros(N_nodes)
d_T = np.zeros(N_nodes)

# === SURFACE NODE (i=0) - CONVECTIVE BC WITH BLOWING ===

i = 0

# Surface properties
k_surf = k_v * (1.0 - alpha_eff[i]) + k_c * alpha_eff[i]

# --- CALCULATE SURFACE BLOWING ---
# Darcy velocity at surface
v_D_surf = (Gamma_perm / mu_g) * (P_new[1] - P_new[0]) / dx

# Surface gas density (using old temperature for stability)
rho_g_surf = P_new[0] / (R_univ * T_old[0])

# Surface mass flux
m_dot_surface = rho_g_surf * v_D_surf

print(f"\n  Surface blowing:")
print(f"    ṁ = {m_dot_surface:.6f} kg/m²·s")

# --- BLOWING CORRECTION (Equations 2.20-2.22) ---

# Stanton number without blowing (Eq. 2.23: ρ*u*C_h = h/cp)
C_h0 = h_0_external / (rho_e * u_e * cp_g)

# Lambda parameter (0.4 laminar, 0.5 turbulent)
lambda_param = 0.5

# Equation 2.22: phi parameter
phi_blowing = calculate_phi_parameter(m_dot_surface, rho_e, u_e, C_h0, lambda_param)

# Equation 2.20: blowing factor
Omega_blowing = calculate_blowing_factor(phi_blowing)

# Equation 2.21: effective heat transfer coefficient
h_eff = Omega_blowing * h_0_external

print(f"  Blowing correction:")
print(f"    φ = {phi_blowing:.6f}")
print(f"    Ω_blowing = {Omega_blowing:.6f}")
print(f"    h_eff = {h_eff:.1f} W/m²·K")
print(f"    Reduction: {(1-Omega_blowing)*100:.1f}%")

# --- CONVECTIVE BOUNDARY CONDITION ---
# Energy balance: h_eff*(T_r - T_w) = k*(T[1] - T[0])/dx
# Rearrange: (k/dx + h_eff)*T[0] - (k/dx)*T[1] = h_eff*T_recovery

a_T[0] = k_surf/dx + h_eff
b_T[0] = k_surf/dx
c_T[0] = 0.0
d_T[0] = h_eff * T_recovery

print(f"  Surface BC coefficients:")
print(f"    a[0] = {a_T[0]:.2e}, b[0] = {b_T[0]:.2e}, d[0] = {d_T[0]:.2e}")

# === INTERIOR NODES (i=1 to N-2) ===

for i in range(1, N_nodes - 1):
    
    k_local = k_v * (1.0 - alpha_eff[i]) + k_c * alpha_eff[i]
    cp_local = cp_v * (1.0 - alpha_eff[i]) + cp_c * alpha_eff[i]
    phi = phi_v * (1.0 - alpha_eff[i]) + phi_c * alpha_eff[i]
    
    rho_gas = P_new[i] / (R_univ * T_old[i])
    rho_c = rho_solid_new[i] * cp_local + rho_gas * phi * cp_g
    a_time = rho_c * dx / dt
    
    a_east = k_local / dx
    a_west = k_local / dx
    
    # Gaz kütle akışı
    m_dot_gas = (P_new[i] / (R_univ * T_old[i])) * (Gamma_perm / mu_g) * \
                (P_new[i+1] - P_new[i-1]) / (2*dx)
    
    h_gas_grad = cp_g * (T_old[i+1] - T_old[i-1]) / (2*dx)
    advection = m_dot_gas * h_gas_grad
    
    # Piroliz kaynağı
    h_solid = cp_local * T_old[i]
    h_eff_pyro = h_solid + rho_solid_new[i] * (cp_v - cp_c) * T_old[i] / (rho_virgin - rho_char)
    h_gas = cp_g * T_old[i]
    
    S_pyro = -(Q_p - h_eff_pyro + h_gas) * drho_dt[i]
    S_adv = -advection
    S_total = S_pyro + S_adv
    
    a_T[i] = a_time + a_east + a_west
    b_T[i] = a_east
    c_T[i] = a_west
    d_T[i] = S_total * dx + a_time * T_old[i]

# === BACK SURFACE (i=N-1) ===

a_T[-1] = 1.0
b_T[-1] = 0.0
c_T[-1] = 0.0
d_T[-1] = T_back

# Thomas algoritması
for i in range(1, N_nodes):
    factor = c_T[i] / a_T[i-1]
    a_T[i] -= factor * b_T[i-1]
    d_T[i] -= factor * d_T[i-1]

T_new = np.zeros(N_nodes)
T_new[-1] = d_T[-1] / a_T[-1]
for i in range(N_nodes - 2, -1, -1):
    T_new[i] = (d_T[i] - b_T[i] * T_new[i+1]) / a_T[i]

print(f"\n  ✓ T: yüzey={T_new[0]:.0f}K, orta={T_new[N_nodes//2]:.0f}K, arka={T_new[-1]:.0f}K")


# =============================================================================
# BÖLÜM 8: ABLASYON HIZI
# =============================================================================
print("\n[7/8] Ablasyon hesaplanıyor...")

# Yüzey Darcy hızı (güncellenmiş basınç ile)
v_D_final = (Gamma_perm / mu_g) * (P_new[1] - P_new[0]) / dx

# Yüzey gaz yoğunluğu (güncellenmiş sıcaklık ile)
rho_gas_surf_final = P_new[0] / (R_univ * T_new[0])

# Yüzey gaz kütle akışı
m_dot_gas_final = rho_gas_surf_final * v_D_final

# Sınır tabakası parametreleri
C_M = 0.002

# Boyutsuz gaz akışı
B_prime_g = m_dot_gas_final / (rho_e * u_e * C_M)

# Char tüketimi (basitleştirilmiş)
B_prime_c = 0.01 * B_prime_g

# Ablasyon hızı
s_dot = B_prime_c * rho_e * u_e * C_M / rho_char

print(f"  ✓ Gaz akışı: {m_dot_gas_final*1e6:.2f} mg/m²·s")
print(f"  ✓ B'_g={B_prime_g:.4f}, B'_c={B_prime_c:.4f}")
print(f"  ✓ Ablasyon: {s_dot*1e6:.2f} μm/s = {s_dot*3600*1000:.2f} mm/saat")


# =============================================================================
# BÖLÜM 9: VERIFICATION - HEAT FLUX BALANCE
# =============================================================================
print("\n[8/8] Isı akısı dengesi kontrol...")

# Convective heat flux
q_conv = h_eff * (T_recovery - T_new[0])

# Conductive heat flux into material
q_cond = -k_surf * (T_new[1] - T_new[0]) / dx

# Balance check
error_percent = abs(q_conv - q_cond) / abs(q_conv) * 100

print(f"  Convective: q_conv = {q_conv/1e6:.3f} MW/m²")
print(f"  Conductive: q_cond = {q_cond/1e6:.3f} MW/m²")
print(f"  Balance error: {error_percent:.2f}%")


# =============================================================================
# ÖZET SONUÇLAR
# =============================================================================
print("\n" + "="*80)
print("SONUÇLAR ÖZETİ - TEK ZAMAN ADIMI (Δt = 0.01s)")
print("="*80)

print(f"\nSICAKLIK:")
print(f"  Yüzey: {T_new[0]:.1f} K ({T_new[0]-273:.0f}°C)")
print(f"  Orta: {T_new[N_nodes//2]:.1f} K ({T_new[N_nodes//2]-273:.0f}°C)")
print(f"  Arka: {T_new[-1]:.1f} K ({T_new[-1]-273:.0f}°C)")

print(f"\nPİROLİZ:")
print(f"  Komponent A: max α = {np.max(alpha_new[0,:]):.4f}")
print(f"  Komponent B: max α = {np.max(alpha_new[1,:]):.4f}")
print(f"  Efektif: max α = {np.max(alpha_eff):.4f}")

print(f"\nBASINÇ:")
print(f"  Yüzey: {P_new[0]/1000:.2f} kPa")
print(f"  Arka: {P_new[-1]/1000:.2f} kPa")

print(f"\nBLOWING EFFECT:")
print(f"  ṁ = {m_dot_surface:.6f} kg/m²·s")
print(f"  φ = {phi_blowing:.6f}")
print(f"  Ω = {Omega_blowing:.6f} ({(1-Omega_blowing)*100:.1f}% reduction)")
print(f"  h_0 = {h_0_external:.0f} W/m²·K → h_eff = {h_eff:.0f} W/m²·K")

print(f"\nABLASYON:")
print(f"  Hız: {s_dot*1e6:.2f} μm/s = {s_dot*3600*1000:.2f} mm/saat")

print("\n" + "="*80)
print("✅ TEK ZAMAN ADIMI TAMAMLANDI")
print("="*80)

print("\n📚 KAYNAKLAR:")
print("Ewing, D. J., Laker, T. S., and Walker, D. G.,")
print('"Numerical Modeling of Ablation Heat Transfer,"')
print("J. Thermophysics and Heat Transfer, Vol. 27, No. 4, 2013")
print("\nBlowing correction: Equations 2.20-2.22 (Semi-empirical approach)")