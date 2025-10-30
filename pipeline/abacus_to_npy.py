import numpy as np

# Makaledeki Tablo I verisi
data = {
    'z': np.array([0.300, 0.400, 1.026, 5.001, 8.017]),
    'M_lcdm': np.array([4.53e12, 4.31e12, 3.20e12, 1.23e12, 0.968e12]),
    'M_dmtde': np.array([4.42e12, 4.20e12, 3.12e12, 1.20e12, 0.944e12]),
    'fd_obs': 0.049,
    'fd_err': 0.005
}

np.save('abacus_27z_data.npy', data)
print("abacus_27z_data.npy OLUŞTURULDU! fd = 0.049 (sabit)")
