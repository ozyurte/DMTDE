#!/usr/bin/env python3
"""
AbacusSummit Verilerini MCMC İçin Hazırlama
===========================================
Tablo I'deki 27 redshift'i CSV/NPY formatına çevirir
"""

import numpy as np
import pandas as pd

# =============================================================================
# Tablo I'deki veriler (makalenizden)
# =============================================================================

# Redshift değerleri
redshifts = np.array([
    0.300, 0.350, 0.400, 0.450, 0.500, 0.575, 0.650, 0.725, 0.800,
    0.875, 0.950, 1.025, 1.100, 1.175, 1.250, 1.325, 1.400, 1.475,
    1.550, 1.625, 1.700, 1.850, 2.000, 2.250, 2.500, 2.750, 3.000
])

# Halo sayıları (milyon cinsinden)
N_halos_millions = np.array([
    92.6, 92.6, 92.5, 92.4, 92.1, 91.5, 90.8, 89.9, 88.9,
    87.6, 86.3, 84.8, 83.1, 81.4, 79.6, 77.6, 75.5, 73.4,
    71.3, 69.0, 66.8, 62.1, 57.4, 49.6, 42.1, 35.0, 28.7
])

# ΛCDM ortalama kütle (10^12 M_sun)
M_LCDM_1e12 = np.array([
    4.53, 4.42, 4.31, 4.20, 4.09, 3.94, 3.80, 3.67, 3.54,
    3.42, 3.31, 3.20, 3.10, 3.00, 2.91, 2.83, 2.74, 2.67,
    2.59, 2.52, 2.46, 2.33, 2.22, 2.06, 1.93, 1.81, 1.71
])

# DMTDE ortalama kütle (10^12 M_sun)
M_DMTDE_1e12 = np.array([
    4.42, 4.31, 4.20, 4.09, 3.99, 3.84, 3.71, 3.58, 3.45,
    3.34, 3.22, 3.12, 3.02, 2.93, 2.84, 2.75, 2.67, 2.60,
    2.53, 2.46, 2.40, 2.28, 2.17, 2.01, 1.88, 1.76, 1.66
])

# Baskılama (%)
suppression_percent = np.full(27, 4.9)

# İstatistiksel belirsizlik (%)
sigma_stat_percent = np.array([
    0.12, 0.12, 0.12, 0.13, 0.13, 0.13, 0.13, 0.14, 0.14,
    0.14, 0.14, 0.15, 0.15, 0.15, 0.16, 0.16, 0.16, 0.17,
    0.17, 0.17, 0.18, 0.18, 0.19, 0.20, 0.22, 0.24, 0.27
])

# =============================================================================
# DataFrame oluşturma
# =============================================================================

df = pd.DataFrame({
    'redshift': redshifts,
    'N_halos': N_halos_millions * 1e6,  # Tam sayıya çevir
    'mean_mass_LCDM': M_LCDM_1e12 * 1e12,  # M_sun cinsine
    'mean_mass_DMTDE': M_DMTDE_1e12 * 1e12,
    'suppression_percent': suppression_percent,
    'sigma_stat_percent': sigma_stat_percent,
    'sigma_mass_LCDM': (M_LCDM_1e12 * 1e12) / np.sqrt(N_halos_millions * 1e6),
    'sigma_mass_DMTDE': (M_DMTDE_1e12 * 1e12) / np.sqrt(N_halos_millions * 1e6)
})

# =============================================================================
# Kaydetme
# =============================================================================

# CSV formatı (insanlar için okunabilir)
df.to_csv('abacus_halos_27z.csv', index=False, float_format='%.6e')
print("✓ CSV kaydedildi: abacus_halos_27z.csv")
print(f"  Boyut: {df.shape}")
print(f"  Sütunlar: {list(df.columns)}")

# NumPy formatı (hızlı yükleme için)
np.save('abacus_halos_27z.npy', df[['redshift', 'mean_mass_LCDM', 
                                     'mean_mass_DMTDE', 'N_halos']].values)
print("✓ NPY kaydedildi: abacus_halos_27z.npy")

# İlk 5 satırı göster
print("\nÖrnek veriler:")
print(df.head())

# İstatistikler
print("\n" + "="*70)
print("VERİ İSTATİSTİKLERİ")
print("="*70)
print(f"Redshift aralığı: {redshifts.min():.3f} - {redshifts.max():.3f}")
print(f"Toplam halo: {df['N_halos'].sum():.2e}")
print(f"Ortalama baskılama: {df['suppression_percent'].mean():.2f}%")
print(f"Baskılama std: {df['suppression_percent'].std():.4f}%")
print("="*70)

# Baskılama faktörünü doğrula
measured_suppression = 100 * (1 - df['mean_mass_DMTDE'] / df['mean_mass_LCDM'])
print(f"\nÖlçülen baskılama (hesaplanan): {measured_suppression.mean():.2f}% ± {measured_suppression.std():.4f}%")
print(f"Teorik beklenti: 4.90%")
print(f"Uyum: ✓" if np.abs(measured_suppression.mean() - 4.9) < 0.01 else "⚠️ Uyumsuzluk!")
