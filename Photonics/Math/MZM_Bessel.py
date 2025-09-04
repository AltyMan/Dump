import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

# Parameters
m_range = np.linspace(0, 10, 500)
detuning = np.linspace(-10, 10, 40000)
mod_indices = [1.4347, 1.8412]  # Modulation indices to show in intensity plots

# Bessel Coefficients for sine and cosine (same for both)
orders = np.arange(0, 5)  # Bessel orders from -10 to 10
orders_even = np.arange(0, 5, 2) # Only even orders for cosine modulation
bessel_sine = {n: jv(n, m_range) for n in orders}
bessel_cosine = {n: jv(n, m_range) for n in orders_even}

# Function to compute intensity spectrum for sine modulation
def compute_sine_spectrum(m_val):
    intensity = np.zeros_like(detuning)
    for n in range(-4, 5):
        Jn = jv(n, m_val)
        intensity += (np.abs(Jn) ** 2) * np.exp(-((detuning - n) ** 2) / 0.0165)
    return intensity

# Function to compute intensity spectrum for cosine modulation
def compute_cosine_spectrum(m_val):
    intensity = np.zeros_like(detuning)
    for n in range(-4, 5, 2):  # Only even orders for cosine
        Jn = jv(n, m_val)
        magnitude_squared = np.abs(Jn) ** 2 #if n == 0 else (2 * Jn) ** 2
        intensity += magnitude_squared * np.exp(-((detuning - n) ** 2) / 0.0165)
        #if n != 0:
            #intensity += magnitude_squared * np.exp(-((detuning - n) ** 2) / 0.0165)
    return intensity

P_in_dbm = -20  # dBm
P_in_mw = 10**(P_in_dbm / 10)  # mW

# Plotting all in one figure
fig, axs = plt.subplots(2, 2, figsize=(14, 8.7))

# Top-left: Bessel coefficients (sine)
for n in orders:
    axs[0, 0].plot(m_range, np.abs(bessel_sine[n]), label=f"|J_{n}(m)|")
axs[0, 0].axvline(x=1.4347, color='gray', linestyle='--', label='m = 1.4347')
axs[0, 0].set_title("Bessel Coefficients (ϕ = m·sin(ωₘt))")
axs[0, 0].set_xlabel("Modulation Index m")
axs[0, 0].set_ylabel("Amplitude |Jₙ(m)|")
axs[0, 0].set_xticks(np.arange(0, 11, 1))
axs[0, 0].legend()
axs[0, 0].grid(True)

# Top-right: Intensity spectra (sine)
for m_val in mod_indices:
    intensity = compute_sine_spectrum(m_val)
    spectrum = intensity * P_in_mw  # Convert to mW
    spectrum_dbm = 10 * np.log10(spectrum + 10**(-50 / 10))  # Convert to dBm
    axs[0, 1].plot(detuning, intensity, label=f"m = {m_val}")
axs[0, 1].set_title("Optical Spectrum (ϕ = m·sin(ωₘt))")
axs[0, 1].set_xlabel("Detuning ( (ω - ω₀) / ωₘ )")
axs[0, 1].set_ylabel("Normalized Intensity")
axs[0, 1].legend()
axs[0, 1].grid(True)

# Bottom-left: Bessel coefficients (cosine)
for n in orders_even:
    axs[1, 0].plot(m_range, np.abs(bessel_cosine[n]), label=f"|J_{n}(m)")
axs[1, 0].axvline(x=1.8412, color='gray', linestyle='--', label='m = 1.8412')
axs[1, 0].set_title("Bessel Coefficients (ϕ = m·cos(ωₘt))")
axs[1, 0].set_xlabel("Modulation Index m")
axs[1, 0].set_ylabel("Amplitude |Jₙ(m)|")
axs[1, 0].set_xticks(np.arange(0, 11, 1))
axs[1, 0].legend()
axs[1, 0].grid(True)

# Bottom-right: Intensity spectra (cosine)
for m_val in mod_indices:
    intensity = compute_cosine_spectrum(m_val)
    spectrum = intensity * P_in_mw  # Convert to mW
    spectrum_dbm = 10 * np.log10(spectrum + 10**(-50 / 10))  # Convert to dBm
    axs[1, 1].plot(detuning, intensity, label=f"m = {m_val}")
axs[1, 1].set_title("Optical Spectrum (ϕ = m·cos(ωₘt))")
axs[1, 1].set_xlabel("Detuning ( (ω - ω₀) / ωₘ )")
axs[1, 1].set_ylabel("Normalized Intensity")
axs[1, 1].legend()
axs[1, 1].grid(True)

# Overall layout
fig.suptitle("Bessel Coefficient and Optical Spectrum Comparisons\nTop: Sine Modulation, Bottom: Cosine Modulation", fontsize=16)
plt.show()
