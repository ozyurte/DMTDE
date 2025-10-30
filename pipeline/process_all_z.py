#!/usr/bin/env python3
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import numpy as np
import os
import json

base_dir = "/home/emreozyurt/abacussummit/AbacusSummit_data/Abacus_base_c000_ph000/halos"
output_file = "all_z_results.txt"
summary_json = "z_summary.json"
particle_mass = 2.109e+10  # Msun

# Otomatik z klasörlerini bul
z_folders = sorted([d for d in os.listdir(base_dir) 
                   if d.startswith('z') and os.path.isdir(os.path.join(base_dir, d))])

print(f"BULUNAN {len(z_folders)} Z-SHIFT: {z_folders}")

results = []
detailed = []

for z_folder in z_folders:
    folder = os.path.join(base_dir, z_folder, "halo_info")
    if not os.path.exists(folder):
        detailed.append(f"{z_folder}: halo_info YOK")
        continue
    
    files = sorted([f for f in os.listdir(folder) if f.endswith('.asdf')])
    if not files:
        detailed.append(f"{z_folder}: .asdf DOSYA YOK")
        continue
    
    total_halos = 0
    mean_masses = []
    z_val = None
    
    print(f"\n=== {z_folder} ({len(files)} dosya) ===")
    for f in files:
        path = os.path.join(folder, f)
        try:
            cat = CompaSOHaloCatalog(path, cleaned=False)
            z_val = cat.header["Redshift"]
            N = cat.halos["N"]
            mass = N * particle_mass
            num_halos = len(mass)
            mean_mass = np.mean(mass)
            mean_dmtde = mean_mass * 0.975
            
            total_halos += num_halos
            mean_masses.append(mean_mass)
            
            print(f"OKUNDU: {f} → {num_halos:,} halo")
            del cat, N, mass
        except Exception as e:
            print(f"{f}: HATA → {e}")
    
    if mean_masses:
        mean_all = np.mean(mean_masses)
        mean_dmtde = mean_all * 0.975
        suppression = 100 * (1 - 0.975**2)
        z_clean = z_folder[1:].replace('_', '.')
        summary = f"{z_folder}: z≈{z_val:.3f}, Halos={total_halos:,}, <M>={mean_all:.2e}, <M_dmtde>={mean_dmtde:.2e}, Bastırma={suppression:.1f}%"
        results.append({
            "z_folder": z_folder,
            "z": z_val,
            "halos": total_halos,
            "mean_mass": mean_all,
            "mean_dmtde": mean_dmtde,
            "suppression": suppression
        })
        detailed.append(f"\n--- {z_folder} ---")
        detailed.append(summary)
    else:
        results.append({"z_folder": z_folder, "error": "VERİ YOK"})

# Dosyaya yaz
with open(output_file, 'w') as f:
    f.write("\n".join(detailed))
with open(summary_json, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTÜM {len(results)} Z-SHIFT TAMAM! → {output_file} + {summary_json}")
print("Şimdi: python3 make_pdf.py → GRAFİKLİ PDF!")
