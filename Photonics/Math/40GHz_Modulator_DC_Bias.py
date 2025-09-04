# IQP 2025 (CQT @ NUS)

import numpy as np
import matplotlib.pyplot as plt

v_range = [x/10 for x in range(-50, 66)]  # Voltage range from 0 to 10V

# Note: Higher marker adds 0.2, lower removes 0.2
power_dbm = [-0.3, -0.6, -0.9, -1.2, -1.5, -1.9, -2.2, -2.6, -3.0, -3.5, -3.9, -4.4, -5.0, -5.5, -6.2, -6.8, -7.5, -8.2, -8.9, -9.6, -10.3, -10.9, -11.3, -11.6, -11.7, -11.5, -10.9, -10.4, -9.7, -9.0, -8.4, -7.8, -7.1, -6.4, -5.7, -5.2, -4.7, -4.2, -3.7, -3.2, -2.8, -2.4, -2.0, -1.6, -1.3, -1.0, -0.7, -0.4, -0.2, 0.08, 0.22, 0.5, 0.7, 0.88, 1.08, 1.28, 1.48, 1.6, 1.72, 1.88, 1.98, 2.08, 2.18, 2.28, 2.32, 2.4, 2.48, 2.48, 2.5, 2.52, 2.52, 2.52, 2.5, 2.48, 2.4, 2.4, 2.32, 2.28, 2.18, 2.1, 2.0, 1.9, 1.78, 1.68, 1.5, 1.38, 1.2, 1.1, 0.9, 0.7, 0.5, 0.3, 0.1, -0.1, -0.4, -0.6, -0.9, -1.2, -1.5, -1.9, -2.2, -2.6, -3.0, -3.5, -3.9, -4.4, -5.0, -5.6, -6.2, -7.0, -7.7, -8.6, -9.5, -10.6, -11.8, -13.1]

plt.figure(figsize=(10, 6))
plt.plot(v_range, power_dbm, marker='o', linestyle='-', color='b', markersize=4)
plt.title("40GHz Modulator DC Bias vs. Power (dBm)")
plt.xlabel("DC Bias Voltage (V)")
plt.ylabel("Output Power (dBm)")
plt.grid(True)
plt.xticks(np.arange(-5, 7, 1))
plt.yticks(np.arange(-15, 4, 1))
plt.xlim(-5, 6.5)
plt.ylim(-15, 4)
plt.legend()
plt.tight_layout()
plt.show()