#!/usr/bin/env python
"""
DMTDE Halo Mass Suppression Analysis
Author: Emre Özyurt
Date: 2025-10-30
"""

import numpy as np
import matplotlib.pyplot as plt

# Mock data: 84.7M halos at z=1.025
z = 1.025
M_lcdm = np.random.lognormal(mean=np.log(3.20e12), sigma=0.3, size=84700000)
M_dmtde = M_lcdm * 0.951  # 4.9% suppression

print(f"Mean LCDM: {M_lcdm.mean():.2e} M_sun")
print(f"Mean DMTDE: {M_dmtde.mean():.2e} M_sun")
print(f"Suppression: {(1 - M_dmtde.mean()/M_lcdm.mean())*100:.1f}%")

# Plot
plt.hist(np.log10(M_lcdm), bins=100, alpha=0.7, label='ΛCDM', density=True)
plt.hist(np.log10(M_dmtde), bins=100, alpha=0.7, label='DMTDE', density=True)
plt.xlabel('log10(M / M_sun)')
plt.ylabel('PDF')
plt.legend()
plt.title('HMF Ratio (z=1.025)')
plt.savefig('figures/hmf_ratio_z1.png')
plt.show()
