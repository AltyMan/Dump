import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

# ==== User & equipment parameters ====
# Optical source
P_in_dbm = -20            # Input optical power (dBm)
P_in_mw = 10**(P_in_dbm/10)
lambda0_nm = 780.0        # Central wavelength (nm)

# Modulator
Vpi = 4.0                 # Vπ (V)
V_RF = 1.827              # RF drive peak (V) for m≈1.4347
V_DC = 0.0                # Example DC offset (V) - adjust as needed
m = np.pi * V_RF / Vpi
phi_DC = np.pi * V_DC / Vpi

# RF frequency & OSA
f_RF = 20e9               # RF frequency (Hz)
c = 3e8                   # Speed of light (m/s)
delta_lambda_nm = (lambda0_nm*1e-9)**2 / c * f_RF * 1e9  # nm
resolution_nm = 0.02      # OSA resolution (nm)
sigma_nm = np.divide(resolution_nm, (2*np.sqrt(2*np.log(2))))
#sigma_nm = resolution_nm / 2.355

# Noise floor
P_noise_dbm = -50         # dBm
P_noise_mw = 10**(P_noise_dbm/10)

# Wavelength axis
wls = np.linspace(lambda0_nm-0.3, lambda0_nm+0.3, 6001)  # nm

# ==== Compute sideband intensities ====
N = 20
orders = np.arange(-N, N+1)
intensities = np.zeros_like(orders, dtype=float)
for i, n in enumerate(orders):
    if n == 0:
        # Carrier term amplitude ∝ J0(m) * cos(phi_DC)
        intensities[i] = np.abs((jv(0, m) * np.cos(phi_DC)))**2
    elif n % 2 == 0:
        # Even-order sideband: ∝ [J_n(m) * cos(phi_DC)]^2
        intensities[i] = np.abs((jv(n, m) * np.cos(phi_DC)))**2
    else:
        # Odd-order sideband: ∝ [J_n(m) * sin(phi_DC)]^2
        intensities[i] = np.abs((jv(n, m) * np.cos(phi_DC)))**2

# ==== Build simulated OSA trace ====
spectrum_mw = np.zeros_like(wls)
for n, frac in zip(orders, intensities):
    center = lambda0_nm + n * delta_lambda_nm
    spectrum_mw += P_in_mw * frac * np.exp(-((wls - center)**2) / (2 * sigma_nm**2))

# Add noise floor
spectrum_mw += P_noise_mw

# Convert to dBm
spectrum_dbm = 10 * np.log10(spectrum_mw)

# ==== Plot ====
plt.figure(figsize=(8, 4))
plt.plot(wls, spectrum_dbm, label="OSA Trace")
plt.axhline(P_noise_dbm, color='gray', linestyle='--', label="Noise Floor")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Power (dBm)")
plt.title(f"MZM Output with DC Offset (V_DC={V_DC} V)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
