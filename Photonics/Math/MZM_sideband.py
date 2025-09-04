import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv
from scipy.optimize import fsolve

# Given parameters from spec sheets
Vpi_RF = 4.0  # RF Port Vπ at 1 GHz from LNX7840A spec (typical)
Vpi_DC = 6.0 # DC Port Vπ at 1 kHz from LNX7840A spec (typical)
lam = 765e-9  # central wavelength (m)
dlam = 0.02e-9  # OSA resolution (m)
c = 3e8  # speed of light (m/s)
driver_gain_db = 26  # typical gain of DR-AN-20-HO

# Solve for modulation index m where J0(m) = J1(m) (so that the sidebands are equal in intensity)
m_solution = fsolve(lambda x: jv(0, x) - jv(2, x), 1.0)[0] # finds the root(s)/solution(s), take first root from solution (only one root anyways)
print(f"Modulation index m ≈ {m_solution:.4f}")

# Calculate the required RF peak voltage for that m
V_peak = m_solution * Vpi_RF / np.pi
print(f"Required RF peak voltage (Vₚ) ≈ {V_peak:.3f} V")

# Convert that to required RF power at modulator input (50-ohm load)
V_rms_mod = V_peak / np.sqrt(2) # RMS voltage
print(f"RMS voltage at modulator input (V_rms_mod) ≈ {V_rms_mod:.3f} V")
P_mod_mw = ((V_rms_mod**2) / 50) * 1000  # in mW
P_mod_dbm = 10 * np.log10(P_mod_mw)
print(f"Required RF power at modulator input ≈ {P_mod_dbm:.2f} dBm")

# Calculate required generator power before driver
P_gen_dbm = P_mod_dbm - driver_gain_db
print(f"Required RF generator power before driver ≈ {P_gen_dbm:.2f} dBm")

# Get the new m value for a set RF power from generator
P_gen_dbm = 0
P_mod_dbm = P_gen_dbm + driver_gain_db  # power at modulator input
P_mod_mw = 10**(P_mod_dbm / 10)  # convert dBm to mW
V_rms_gen = np.sqrt((P_mod_mw  / 1000) * 50) # RMS voltage at generator output
V_peak_gen = V_rms_gen * np.sqrt(2)  # peak voltage at generator output
m_gen = (V_peak_gen * np.pi) / Vpi_RF  # modulation index for this generator power
print(f"Modulation index m for generator power {P_gen_dbm} dBm ≈ {m_gen:.4f}")

# Determine the minimum RF frequency to see sidebands separated by Δλ
f_RF_min = (c * dlam) / (lam**2)
print(f"Minimum RF frequency f_RF ≈ {f_RF_min/1e9:.2f} GHz")

# Gaussian lineshape
def gaussian(wl, center, amp, sigma):
    return amp * np.exp(-((wl - center)**2) / (2 * sigma**2))

# Build total power spectrum (mW)
def compute_spectrum(wls, m, phi_DC, orders, lambda0, delta_lambda, P_in_mw, sigma):
    spectrum = np.zeros_like(wls)
    for n in orders:
        Jn = jv(n, m)  # Bessel function value for order n
        if n == 0 or n % 2 == 0:
            intensity = np.abs(Jn * np.cos(phi_DC))**2
        else:
            intensity = np.abs(Jn * np.sin(phi_DC))**2
        center_wl = lambda0 + n * delta_lambda
        g_spectrum = gaussian(wls, center_wl, intensity, sigma)
        spectrum += g_spectrum  # accumulate contributions from all orders
    spectrum_mw = P_in_mw * spectrum  # scale by input power
    return spectrum_mw

# Given parameters
m = 1.8412 # modulation index for π/4 DC bias and no DC bias

lambda0 = [765.0, 780.0]  # nm
f_RF = 1.5e10
delta_lambda = np.multiply((f_RF * np.pow(np.multiply(lambda0, 10**-9), 2)) / c, 1e9)  # nm
resolution = 0.02  # nm
sigma = np.divide(delta_lambda, (2*np.sqrt(2*np.log(2))))  # Gaussian sigma (standard deviation) for 0.02 nm FWHM

# Input optical power
P_in_dbm = -15 # dBm
P_in_mw = 10**(P_in_dbm / 10)  # mW

# Input DC bias voltage
V_DC = [0, 1.5, 3]  # V, example DC offset (can be adjusted)
phi_DC = np.divide(np.multiply(np.pi, V_DC), Vpi_DC)  # phase offset from DC bias

# Noise floor
P_noise_dbm = -50  # dBm
P_noise_mw = 10**(P_noise_dbm / 10)  # mW

# Calculate Bessel intensities for orders -N to N
N = 20  # include orders up to ±20
sine_orders = np.arange(-N, N + 1)
cosine_orders = sine_orders[::2]

# Wavelength axis
wls_765 = np.linspace(lambda0[0] - 0.3, lambda0[0] + 0.3, 40001)  # nm
wls_780 = np.linspace(lambda0[1] - 0.3, lambda0[1] + 0.3, 40001)  # nm

# Spectrum intensities for sine and cosine modulation
spectrum_765_none_mw = compute_spectrum(wls_765, m, phi_DC[0], cosine_orders, lambda0[0], delta_lambda[0], P_in_mw, sigma[0])
spectrum_765_quarter_mw = compute_spectrum(wls_765, m, phi_DC[1], sine_orders, lambda0[0], delta_lambda[0], P_in_mw, sigma[0])
spectrum_765_half_mw = compute_spectrum(wls_765, m, phi_DC[2], sine_orders, lambda0[0], delta_lambda[0], P_in_mw, sigma[0])
spectrum_780_none_mw = compute_spectrum(wls_780, m, phi_DC[0], cosine_orders, lambda0[1], delta_lambda[1], P_in_mw, sigma[1])
spectrum_780_quarter_mw = compute_spectrum(wls_780, m, phi_DC[1], sine_orders, lambda0[1], delta_lambda[1], P_in_mw, sigma[1])
spectrum_780_half_mw = compute_spectrum(wls_780, m, phi_DC[2], sine_orders, lambda0[1], delta_lambda[1], P_in_mw, sigma[1])

# Add noise floor
spectrum_765_none_mw += P_noise_mw
spectrum_765_quarter_mw += P_noise_mw
spectrum_765_half_mw += P_noise_mw
spectrum_780_none_mw += P_noise_mw
spectrum_780_quarter_mw += P_noise_mw
spectrum_780_half_mw += P_noise_mw

# Convert to dBm
spectrum_765_none_dbm = 10 * np.log10(spectrum_765_none_mw)
spectrum_765_quarter_dbm = 10 * np.log10(spectrum_765_quarter_mw)
spectrum_765_half_dbm = 10 * np.log10(spectrum_765_half_mw)
spectrum_780_none_dbm = 10 * np.log10(spectrum_780_none_mw)
spectrum_780_quarter_dbm = 10 * np.log10(spectrum_780_quarter_mw)
spectrum_780_half_dbm = 10 * np.log10(spectrum_780_half_mw)

# Plot
fig, axs = plt.subplots(1, 2, figsize=(16, 4))
axs[0].plot(wls_765, spectrum_765_none_dbm, label='No DC', color='blue')
axs[0].plot(wls_765, spectrum_765_quarter_dbm, label='Quarter DC', color='orange')
axs[0].plot(wls_765, spectrum_765_half_dbm, label='Half DC', color='green')
axs[0].axhline(P_noise_dbm, color='gray', linestyle='--', label='Noise Floor (-50 dBm)')
axs[0].axvline(764.9, color='gray', linestyle='--', alpha=0.5)
axs[0].axvline(765.1, color='gray', linestyle='--', alpha=0.5)
axs[0].set_xlabel("Wavelength (nm)")
axs[0].set_ylabel("Power (dBm)")
axs[0].set_title(f"{lambda0[0]:.0f}nm Modulation")
axs[0].legend()
axs[0].grid(True)
axs[1].plot(wls_780, spectrum_780_none_dbm, label='No DC', color='blue')
axs[1].plot(wls_780, spectrum_780_quarter_dbm, label='Quarter DC', color='orange')
axs[1].plot(wls_780, spectrum_780_half_dbm, label='Half DC', color='green')
axs[1].axhline(P_noise_dbm, color='gray', linestyle='--', label='Noise Floor (-50 dBm)')
axs[1].axvline(779.9, color='gray', linestyle='--', alpha=0.5)
axs[1].axvline(780.1, color='gray', linestyle='--', alpha=0.5)
axs[1].set_xlabel("Wavelength (nm)")
axs[1].set_ylabel("Power (dBm)")
axs[1].set_title(f"{lambda0[1]:.0f}nm Modulation")
axs[1].legend()
axs[1].grid(True)
fig.suptitle(f"Simulated OSA Trace at m = {m}, f = {f_RF/1e9:.0f}GHz\nLeft: {lambda0[0]:.0f}nm Modulation, Right: {lambda0[1]:.0f}nm Modulation")
fig.tight_layout()
plt.show()