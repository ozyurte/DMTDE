#!/usr/bin/env python3
import numpy as np
from abacusnbody.halos import read_halo_catalog  # Düzeltilmiş import (halos modülü)
from pylians import power_spectrum
import os

# Veri yolları (z=0.5 için)
data_dir = '/home/emreozyurt/abacussummit/AbacusSummit_data/Abacus_base_c000_ph000'
sim_name = 'base_c000_ph000'
redshift = 2.5  # z2.500 klasörü için (0.5 değil, Abacus redshift formatı)

print("Halo verisi yükleniyor...")
try:
    # Haloları oku
    halos = read_halo_catalog(data_dir, sim_name=sim_name, redshift=redshift)
    # Pozisyonları al (alan isimleri değişebilir)
    if 'x' in halos.dtype.names:
        pos = np.column_stack((halos['x'], halos['y'], halos['z']))
    elif 'Position' in halos.dtype.names:
        pos = halos['Position']
    else:
        raise ValueError("Pozisyon alanı bulunamadı!")
    if pos.shape[0] > 100000:  # Subsample
        pos = pos[:100000]
    print(f"Pozisyon yüklendi: {pos.shape}")
except Exception as e:
    print(f"Ver
