import asdf
import h5py
import numpy as np
import os

# Blosc sıkıştırmasını manuel olarak tanımla
try:
    import blosc
    asdf.config.AsdfConfig().add_extension(
        asdf.extension.CompressionExtension(
            compression_types=['blsc'],
            compress=blosc.compress,
            decompress=blosc.decompress,
            types=[np.ndarray]
        )
    )
except ImportError as e:
    print(f"Blosc kayıt hatası: {e}. Lütfen 'pip3 install blosc' çalıştırın.")
    exit(1)

# Veri klasörü
folder = "/home/emreozyurt/abacussummit/AbacusSummit_data/Abacus_base_c000_ph000/halos/z1.025/halo_info"
output_folder = "/home/emreozyurt/abacussummit/AbacusSummit_data/Abacus_base_c000_ph000/halos/z1.025/halo_info_hdf5"

# Çıkış klasörünü oluştur
os.makedirs(output_folder, exist_ok=True)

for file in os.listdir(folder):
    if file.endswith('.asdf'):
        file_path = os.path.join(folder, file)
        try:
            with asdf.open(file_path, lazy_load=True) as af:
                z = af['header']['redshift']
                mass = af['data']['halos']['mass']
                dgamma2 = 0.975  # z=1.025 tahmini
                mass_dmtde = mass * dgamma2

                # HDF5'e kaydet
                output_file = os.path.join(output_folder, file.replace('.asdf', '_dmtde.h5'))
                with h5py.File(output_file, 'w') as hf:
                    hf.create_dataset('redshift', data=z, compression='gzip', compression_opts=4)
                    hf.create_dataset('mass', data=mass, compression='gzip', compression_opts=4)
                    hf.create_dataset('mass_dmtde', data=mass_dmtde, compression='gzip', compression_opts=4)
                    print(f"{file}: z={z}, Mean Halo Mass (DMTDE): {np.mean(mass_dmtde):.2e} Msun")
        except Exception as e:
            print(f"Hata: {file} işlenirken sorun oluştu: {str(e)}")
