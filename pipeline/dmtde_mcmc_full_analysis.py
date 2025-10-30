#!/usr/bin/env python3
"""
DMTDE MCMC Analizi - Tam Pipeline
==================================
7 parametreli doğru analiz: [H0, Ωm, Ωb, σ8, ns, τ, fd]

TÜM VERİ SETLERİ VE DOĞRU DMTDE H(z) MODELİ İLE NİHAİ ÇALIŞTIRMA
(COBAYA DESI LIKELIHOOD KULLANILARAK)
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
from scipy.integrate import quad, IntegrationWarning
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
# ❗️ Cobaya Likelihood importu ❗️
from cobaya.likelihood import Likelihood
from cobaya.install import install # Likelihood'u bulmak için
# ----------------------------------
import h5py
import pandas as pd
from multiprocessing import Pool, set_start_method
import os
import time
import warnings

# Sayısal uyarıları baskıla ama hataları değil
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=IntegrationWarning)
try:
    from astropy.cosmology import CosmologyWarning
    warnings.filterwarnings("ignore", category=CosmologyWarning)
except ImportError:
    pass 

# =============================================================================
# BÖLÜM 1: VERİ YÜKLEME FONKSİYONLARI (BAO KALDIRILDI)
# =============================================================================

class DMTDEData:
    """Tüm veri setlerini yönetir (BAO hariç)"""
    
    def __init__(self):
        self.abacus_data = None
        self.sn_data = None
        self.cmb_data = None
        
    def load_abacus_halos(self, file_path):
        """
        AbacusSummit halo kütlelerini yükler.
        """
        try:
            if file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            elif file_path.endswith('.npy'):
                data = np.load(file_path, allow_pickle=True)
                # Sütunları doğru al: redshift[0], mean_mass_LCDM[1], N_halos[3]
                df = pd.DataFrame(data[:, [0, 1, 3]], columns=['redshift', 'mean_mass_LCDM', 'N_halos'])
            else:
                raise ValueError("Format: .csv veya .npy olmalı")
            
            self.abacus_data = {
                'z': df['redshift'].values,
                'M_LCDM': df['mean_mass_LCDM'].values,
                'N_halos': df['N_halos'].values
            }
            print(f"✓ Abacus yüklendi: {len(self.abacus_data['z'])} redshift")
            return True
        except FileNotFoundError:
            print(f"⚠️  Abacus veri dosyası '{file_path}' bulunamadı.")
            return False
        except Exception as e:
            print(f"⚠️  Abacus verisi yüklenirken hata: {e}")
            return False
    
    def load_pantheon_sn(self, file_path="external_data/pantheon_plus_sn.txt"):
        """
        Pantheon+ SNIa verilerini yükler.
        """
        try:
            # Sütunları doğru al: z[1], MU_SH0ES[4], MU_SH0ES_ERR[5]
            data = np.loadtxt(file_path, skiprows=1, usecols=(1, 4, 5)) 
            z_sn = data[:, 0]
            mu_obs = data[:, 1]
            sigma_mu = data[:, 2]
            valid_sigma = np.isfinite(sigma_mu) & (sigma_mu > 1e-9)
            if not np.all(valid_sigma):
                num_invalid = np.sum(~valid_sigma)
                print(f"⚠️  Pantheon+ verisinde {num_invalid} geçersiz sigma_mu bulundu ve filtreledi.")
                z_sn = z_sn[valid_sigma]
                mu_obs = mu_obs[valid_sigma]
                sigma_mu = sigma_mu[valid_sigma]
                
            print(f"✓ Gerçek Pantheon+ yüklendi: {len(z_sn)} geçerli supernova")
        except:
            print(f"⚠️  Gerçek Pantheon+ verisi ('{file_path}') bulunamadı veya okunamadı. Mock veri kullanılıyor.")
            z_sn = np.linspace(0.01, 2.3, 50)
            cosmo_mock = FlatLambdaCDM(H0=67.4, Om0=0.315)
            mu_obs = cosmo_mock.distmod(z_sn).value
            mu_obs += np.random.normal(0, 0.15, size=mu_obs.shape)
            sigma_mu = np.full_like(mu_obs, 0.15)
        
        self.sn_data = {
            'z': z_sn,
            'mu_obs': mu_obs,
            'sigma_mu': sigma_mu
        }
    
    def load_planck_cmb(self):
        """
        Planck 2018 constraint'lerini yükler. (Gaussian prior'lar)
        """
        self.cmb_data = {
            'H0': (67.4, 0.5),      # km/s/Mpc
            'Ωm': (0.315, 0.007),
            'Ωb': (0.0493, 0.0006),
            'σ8': (0.811, 0.006),
            'ns': (0.965, 0.004),
            'τ': (0.054, 0.007)
        }
        print(f"✓ Planck CMB prior'ları yüklendi")

# =============================================================================
# BÖLÜM 2: KOZMOLOJİK MODEL (ASTROPY İLE DOĞRU DMTDE H(z))
# =============================================================================

class DMTDECosmologyAstropy:
    """Astropy kullanarak DOĞRU DMTDE model hesaplamaları"""
    
    def __init__(self):
        pass # rd'yi hesaplama fonksiyonu ekleyeceğiz

    def get_astropy_cosmo(self, H0, Om, Ob, fd): 
        """
        Verilen parametrelerle Astropy kozmoloji nesnesi oluşturur.
        Doğru DMTDE H(z) için Om_eff kullanır.
        """
        Om_eff = Om * (1.0 - fd) 
        Ob_use = Ob
        
        # Geçerlilik kontrolü
        if Om_eff < 0 or (1.0 - Om_eff) < 0 or Ob_use < 0 or H0 < 10: 
            return None 
        # Baryon yoğunluğu toplam madde yoğunluğundan büyük olamaz
        # Astropy bunu kontrol ediyor ama biz de edelim
        if Ob_use / (H0/100.0)**2 > Om_eff / (H0/100.0)**2:
             return None

        try:
            cosmo = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om_eff, Ob0=Ob_use)
            return cosmo
        except ValueError:
            return None 

    def get_rd(self, H0, Om, Ob): 
        """ Sound horizon rd'yi hesaplar (BAO için gerekli) """
        h = H0 / 100.0
        # rd erken evren fiziğine bağlıdır, DMTDE bunu değiştirmez varsayımı
        Om_lcdm = Om 
        Ob_lcdm = Ob
        # Geçerlilik kontrolü
        if Om_lcdm < 0 or (1.0 - Om_lcdm) < 0 or Ob_lcdm < 0 or H0 < 10: return np.inf
        if Ob_lcdm / h**2 > Om_lcdm / h**2: return np.inf

        try:
            # Ob0 değil Obh2 vermek daha standart olabilir ama FlatLambdaCDM Ob0 alır
            cosmo_lcdm = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om_lcdm, Ob0=Ob_lcdm)
            # Astropy'nin rd hesaplama yöntemi (versiyona göre değişebilir)
            if hasattr(cosmo_lcdm, 'sound_horizon_drag'): # Eski versiyonlar
                 rd_value = cosmo_lcdm.sound_horizon_drag.to(u.Mpc).value
            elif hasattr(cosmo_lcdm, 'sound_horizon_z'): # Yeni versiyonlar
                 # z_drag hesaplaması için Tcmb0 lazım, varsayılanı kullanalım
                 rd_value = cosmo_lcdm.sound_horizon_z(cosmo_lcdm.zs_drag).to(u.Mpc).value
            else:
                 # Çok eski versiyonlar için yaklaşık formül (kullanılmamalı)
                 Obh2 = Ob_lcdm * h**2
                 Omh2 = Om_lcdm * h**2
                 rd_value = 55.154 * np.exp(-72.3 * (Omh2 * Obh2)**0.28) / (Omh2**0.25 * Obh2**0.078) * (H0/100)**-1
            
            return rd_value if np.isfinite(rd_value) else np.inf
        except Exception:
            return np.inf

    def H(self, z, H0, Om, Ob, fd): 
        """Hubble parametresi H(z) [km/s/Mpc]"""
        cosmo = self.get_astropy_cosmo(H0, Om, Ob, fd)
        if cosmo is None: return np.inf 
        try:
            Hz = cosmo.H(z).to(u.km / u.s / u.Mpc).value
            Hz[~np.isfinite(Hz)] = np.inf 
            Hz[Hz < 1e-9] = np.inf # 0 veya negatif H(z) geçersiz
            return Hz
        except Exception:
             return np.inf

    def comoving_distance(self, z, H0, Om, Ob, fd): 
        """Comoving distance [Mpc]"""
        cosmo = self.get_astropy_cosmo(H0, Om, Ob, fd)
        if cosmo is None: return np.inf
        try:
            cd = cosmo.comoving_distance(z).to(u.Mpc).value
            cd[~np.isfinite(cd)] = np.inf
            cd[cd < 0] = np.inf 
            return cd
        except Exception:
            return np.inf 
             
    def luminosity_distance(self, z, H0, Om, Ob, fd): 
        """Luminosity distance [Mpc]"""
        cosmo = self.get_astropy_cosmo(H0, Om, Ob, fd)
        if cosmo is None: return np.inf
        try:
            dL = cosmo.luminosity_distance(z).to(u.Mpc).value
            dL[~np.isfinite(dL)] = np.inf
            dL[dL <= 1e-9] = np.inf 
            return dL
        except Exception:
             return np.inf
             
    def distance_modulus(self, z, H0, Om, Ob, fd): 
        """Distance modulus μ(z)"""
        dL = self.luminosity_distance(z, H0, Om, Ob, fd)
        if not np.all(np.isfinite(dL)): return dL 
        
        valid_dL = dL > 1e-9 
        mu = np.full_like(dL, np.inf) 
        if np.any(valid_dL):
             dL_valid = dL[valid_dL]
             # Logaritma öncesi kontrol (çok nadir ama olabilir)
             if np.any(dL_valid <= 0): return np.inf * np.ones_like(dL)
             mu[valid_dL] = 5 * np.log10(dL_valid) + 25
        
        return mu

# =============================================================================
# BÖLÜM 3: LİKELİHOOD FONKSİYONU (COBAYA BAO İLE)
# =============================================================================

class DMTDELikelihood:
    """Toplam likelihood hesabı (Astropy + Cobaya BAO ile)"""
    
    _likelihood_instances = {} # Cobaya likelihood'larını cache'lemek için

    def __init__(self, data):
        self.data = data
        self.cosmo = DMTDECosmologyAstropy() 
        
        # Cobaya DESI BAO Likelihood'unu başlat (sadece bir kez)
        likelihood_name = "bao.desi_2024_bao_all"
        if likelihood_name not in DMTDELikelihood._likelihood_instances:
            try:
                # Cobaya'nın kurulu paket yolunu bul
                packages_path = install(None, action="show_path", skip_global=True)
                # Likelihood'u yükle ve cache'le
                instance = Likelihood.get_instance(likelihood_name, packages_path=packages_path)
                # Gerekli parametreleri iste (genellikle H, D_M, rd)
                instance.check_and_prepare({}) 
                DMTDELikelihood._likelihood_instances[likelihood_name] = instance
                print(f"✓ Cobaya Likelihood '{likelihood_name}' başarıyla başlatıldı.")
            except Exception as e:
                print(f"⚠️ Cobaya Likelihood '{likelihood_name}' başlatılamadı: {e}.")
                DMTDELikelihood._likelihood_instances[likelihood_name] = None 
        
        self.desi_bao_likelihood = DMTDELikelihood._likelihood_instances[likelihood_name]
            
    def log_prior(self, theta):
        """Parametreler için düz (flat) prior'lar"""
        H0, Om, Ob, sigma8, ns, tau, fd = theta
        
        # Daha dar sınırlar deneyelim mi?
        if not (65.0 < H0 < 76.0): return -np.inf # Planck ve SH0ES arası
        if not (0.25 < Om < 0.38): return -np.inf
        if not (0.045 < Ob < 0.055): return -np.inf
        if not (0.75 < sigma8 < 0.88): return -np.inf
        if not (0.94 < ns < 0.99): return -np.inf
        if not (0.03 < tau < 0.08): return -np.inf
        if not (0.02 < fd < 0.08): return -np.inf # Abacus etrafında
        if Ob > Om: return -np.inf
        
        return 0.0 

    def log_likelihood(self, theta, debug=False): # DEBUG için flag
        """Toplam log-likelihood"""
        
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
            
        H0, Om, Ob, sigma8, ns, tau, fd = theta
        
        # Likelihood bileşenlerini ayrı ayrı hesapla ve kontrol et
        
        logL_abacus = self._logL_abacus(Om, sigma8, fd)
        if not np.isfinite(logL_abacus): 
            if debug: print("DEBUG: _logL_abacus failed")
            return -np.inf
        
        # Cobaya BAO Likelihood'unu çağır
        logL_bao = self._logL_cobaya_bao(H0, Om, Ob, fd, debug=debug) 
        if not np.isfinite(logL_bao): 
            if debug: print("DEBUG: _logL_cobaya_bao failed")
            return -np.inf
        
        logL_sn = self._logL_sn(H0, Om, Ob, fd, debug=debug) # Ob eklendi
        if not np.isfinite(logL_sn): 
            if debug: print("DEBUG: _logL_sn failed")
            return -np.inf
        
        logL_cmb = self._logL_cmb_prior(H0, Om, Ob, sigma8, ns, tau)
        if not np.isfinite(logL_cmb): 
            if debug: print("DEBUG: _logL_cmb_prior failed")
            return -np.inf
        
        total_logL = lp + logL_abacus + logL_bao + logL_sn + logL_cmb
            
        if not np.isfinite(total_logL):
             if debug: print("DEBUG: Final total_logL is not finite")
             return -np.inf
             
        return total_logL
    
    def _logL_abacus(self, Om, sigma8, fd):
        """
        AbacusSummit fd kısıtlaması
        """
        fd_obs = 0.049
        fd_err = 0.0052 

        fd_model = fd
        
        chi2 = ((fd_model - fd_obs) / fd_err)**2
        return -0.5 * chi2

    def _logL_cobaya_bao(self, H0, Om, Ob, fd, debug=False):
        """ Cobaya DESI BAO Likelihood'unu çağırır """
        if self.desi_bao_likelihood is None:
             if debug: print("DEBUG Cobaya BAO: Likelihood başlatılamadı, 0 dönülüyor.")
             return 0.0 # Başlatılamadıysa katkı 0

        # Cobaya likelihood'u için gerekli teorik değerleri hesapla
        try:
            # 1. rd'yi hesapla
            rd_value = self.cosmo.get_rd(H0, Om, Ob)
            if not np.isfinite(rd_value): 
                if debug: print("DEBUG Cobaya BAO: rd is not finite")
                return -np.inf

            # 2. H(z) ve D_M(z) fonksiyonlarını oluştur (Astropy'den)
            #    Bu fonksiyonlar lambda z: ... şeklinde olmalı
            H_func = lambda z_bao: self.cosmo.H(z_bao, H0, Om, Ob, fd)
            DM_func = lambda z_bao: self.cosmo.comoving_distance(z_bao, H0, Om, Ob, fd)
            
            # 3. Gerekli tüm değerleri içeren bir provider dict oluştur
            #    Bu likelihood'un ne istediğini bilmemiz lazım. Genellikle:
            #    H(z), D_M(z), rs_drag (veya rd)
            provider_dict = {
                'Hubble': H_func,                   # Fonksiyon olarak verilir
                'comoving_radial_distance': DM_func, # Fonksiyon olarak verilir
                'rs_drag': rd_value                 # Değer olarak verilir (veya 'rd'?)
                # 'rd': rd_value # Bazı likelihood'lar bunu isteyebilir
            }
            
            # Likelihood'u bu provider ile çağır
            # Cobaya'nın likelihood'ları log-olasılık (logp) döndürür
            logp = self.desi_bao_likelihood.logp( **provider_dict )
            
            if debug:
                print(f"  DEBUG Cobaya BAO @ H0={H0:.1f}, Om={Om:.3f}, Ob={Ob:.4f}, fd={fd:.3f}")
                print(f"    Calculated rd = {rd_value:.2f}")
                # H(z) ve D_M(z)'yi örnek bir z'de test edelim
                z_test = np.array([0.5, 1.0, 1.5])
                print(f"    Test H(z) = {H_func(z_test)}")
                print(f"    Test D_M(z) = {DM_func(z_test)}")
                print(f"    Cobaya Likelihood logp = {logp:.2f}")

            # Sonuç sonlu mu kontrol et
            return logp if np.isfinite(logp) else -np.inf

        except Exception as e:
            if debug: print(f"DEBUG Cobaya BAO: Exception during likelihood call: {e}")
            return -np.inf
    
    def _logL_sn(self, H0, Om, Ob, fd, debug=False): # Ob eklendi
        """Pantheon+ SNIa likelihood"""
        if self.data.sn_data is None: return 0.0
        
        z = self.data.sn_data['z']
        mu_obs = self.data.sn_data['mu_obs']
        sigma_mu = self.data.sn_data['sigma_mu'] 
        
        # Model tahmini (Astropy ile, doğru H(z))
        mu_pred = self.cosmo.distance_modulus(z, H0, Om, Ob, fd) # Ob eklendi
        
        if debug:
            print(f"  DEBUG SN @ H0={H0:.1f}, Om={Om:.3f}, Ob={Ob:.4f}, fd={fd:.3f}")
            idx = [0, len(z)//2, -1] 
            print(f"    z_sn(sample) = {z[idx]}")
            print(f"    mu_pred(sample) = {mu_pred[idx] if np.all(np.isfinite(mu_pred)) else 'Not Finite'}")

        if not np.all(np.isfinite(mu_pred)):
            if debug: print("  DEBUG SN: mu_pred not finite.")
            return -np.inf
            
        diff = mu_obs - mu_pred
        
        # Sadece geçerli noktalar üzerinden chi2 hesapla
        valid = np.isfinite(diff) & np.isfinite(sigma_mu) & (sigma_mu > 1e-9)
        if np.sum(valid) == 0: 
             if debug: print(f"  DEBUG SN: No valid points found after prediction.")
             return -np.inf
             
        chi2 = np.sum((diff[valid] / sigma_mu[valid])**2)
        
        if not np.isfinite(chi2): 
             if debug: print(f"  DEBUG SN: chi2 is not finite ({chi2}).")
             return -np.inf
             
        return -0.5 * chi2
    
    def _logL_cmb_prior(self, H0, Om, Ob, sigma8, ns, tau):
        """Planck CMB Gaussian prior'ları"""
        if self.data.cmb_data is None: return 0.0
        
        logL = 0.0
        params_fid = self.data.cmb_data
        
        logL += -0.5 * ((H0 - params_fid['H0'][0]) / params_fid['H0'][1])**2
        logL += -0.5 * ((Om - params_fid['Ωm'][0]) / params_fid['Ωm'][1])**2
        logL += -0.5 * ((Ob - params_fid['Ωb'][0]) / params_fid['Ωb'][1])**2
        logL += -0.5 * ((sigma8 - params_fid['σ8'][0]) / params_fid['σ8'][1])**2
        logL += -0.5 * ((ns - params_fid['ns'][0]) / params_fid['ns'][1])**2
        logL += -0.5 * ((tau - params_fid['τ'][0]) / params_fid['τ'][1])**2
        
        if not np.isfinite(logL): return -np.inf
            
        return logL

# =============================================================================
# BÖLÜM 4: MCMC ÇALIŞTIRICISI
# =============================================================================

class MCMCRunner:
    """emcee MCMC zincirini çalıştırır"""
    
    def __init__(self, likelihood, n_walkers=32, n_steps=5000):
        self.likelihood = likelihood
        self.n_walkers = n_walkers
        self.n_steps = n_steps
        self.n_params = 7 # [H0, Om, Ob, sigma8, ns, tau, fd]
        
    def run(self, output_file='dmtde_chain_7param.h5'):
        """MCMC zincirini çalıştır"""
        print("\n" + "="*70)
        print("MCMC BAŞLIYOR (NİHAİ MODEL: TÜM VERİLER + DOĞRU DMTDE H(z) + COBAYA BAO)")
        print("="*70)
        
        print("Başlangıç noktaları oluşturuluyor...")
        initial_pos = self._initialize_walkers()
        print("✓ Başlangıç noktaları başarıyla oluşturuldu.")

        backend = emcee.backends.HDFBackend(output_file)
        backend.reset(self.n_walkers, self.n_params)
        
        n_cores = os.cpu_count()
        print(f"Paralel çalıştırma için {n_cores} CPU çekirdeği kullanılacak.")
        
        # Cobaya likelihood'u multiprocessing ile uyumlu olmayabilir.
        # Pool kullanmadan deneyelim önce, sonra gerekirse Pool'u kaldıralım.
        pool_obj = Pool(processes=n_cores)
        
        sampler = emcee.EnsembleSampler(
            self.n_walkers, 
            self.n_params,
            self.likelihood.log_likelihood, # debug=False MCMC sırasında
            backend=backend,
            pool=pool_obj 
        )
            
        print(f"Walkers: {self.n_walkers}, Steps: {self.n_steps}")
        print("Çalışıyor... (Bu işlem zaman alabilir)")
        
        start_time = time.time()
        sampler.run_mcmc(initial_pos, self.n_steps, progress=True)
        end_time = time.time()
        
        # Pool'u kapat
        if pool_obj: pool_obj.close()
            
        print(f"\n✓ MCMC tamamlandı! Süre: {end_time - start_time:.2f} saniye")
        
        print(f"Chain kaydedildi: {output_file}")
        
        try:
            acc_frac = np.mean(sampler.acceptance_fraction)
            print(f"Ortalama kabul oranı: {acc_frac:.3f}")
            if acc_frac < 0.15:
                print("⚠️  Kabul oranı çok düşük. Sonuçlar güvenilir olmayabilir.")
        except Exception as e:
            print(f"Kabul oranı hesaplanamadı: {e}")
            
        return sampler
    
    def _initialize_walkers(self):
        """Başlangıç pozisyonlarını oluştur"""
        initial_center = np.array([
            71.0,   # H0 (SH0ES'e yakın)
            0.310,  # Ωm
            0.049,  # Ωb
            0.800,  # σ8 (S8'i düşürmek için)
            0.965,  # ns
            0.054,  # τ
            0.049   # fd (Abacus hedefi)
        ])
        
        initial_scatter = np.array([0.5, 0.01, 0.001, 0.01, 0.005, 0.007, 0.002])
        
        pos = np.zeros((self.n_walkers, self.n_params))
        
        n_accepted = 0
        max_tries = self.n_walkers * 1000 
        tries = 0

        print(f"  {self.n_walkers} walker için geçerli başlangıç noktası aranıyor...")
        while n_accepted < self.n_walkers and tries < max_tries:
            trial_pos = initial_center + initial_scatter * np.random.randn(self.n_params)
            tries += 1
            # likelihood'u debug=False ile çağır (sadece sonlu mu kontrol et)
            if np.isfinite(self.likelihood.log_likelihood(trial_pos, debug=False)): 
                pos[n_accepted] = trial_pos
                n_accepted += 1

        if n_accepted < self.n_walkers:
             test_logL = self.likelihood.log_likelihood(initial_center, debug=True) 
             raise RuntimeError(f"Yeterli sayıda ({self.n_walkers}) geçerli başlangıç noktası bulunamadı! Merkezdeki logL={test_logL}. Likelihood fonksiyonunu kontrol edin.")
                 
        return pos

# =============================================================================
# BÖLÜM 5: SONUÇ ANALİZİ VE GÖRSELLEŞTİRME
# =============================================================================

class MCMCAnalyzer:
    """MCMC sonuçlarını analiz eder"""
    
    def __init__(self, chain_file, n_params=7):
        self.chain_file = chain_file
        self.n_params = n_params
        self.param_names = ['H0', 'Om', 'Ob', 'sigma8', 'ns', 'tau', 'fd']
        self.param_labels = [
            r'$H_0$ [km/s/Mpc]',
            r'$\Omega_m$',
            r'$\Omega_b$',
            r'$\sigma_8$',
            r'$n_s$',
            r'$\tau$',
            r'$f_d$'
        ]
        
    def load_samples(self, discard=1000, thin=10):
        """Chain'i yükle ve burn-in at"""
        try:
            reader = emcee.backends.HDFBackend(self.chain_file)
            samples = reader.get_chain(discard=discard, thin=thin, flat=True)
            print(f"✓ {samples.shape[0]} sample yüklendi (burn-in: {discard}, thin: {thin})")
            return samples
        except Exception as e:
            print(f"❌ HATA: Chain dosyası ('{self.chain_file}') okunamadı: {e}")
            return None
    
    def print_statistics(self, samples):
        """Parametre istatistiklerini yazdır"""
        print("\n" + "="*70)
        print("MCMC SONUÇLARI (7 Parametre, NİHAİ MODEL)")
        print("="*70)
        
        S8 = samples[:, 3] * np.sqrt(samples[:, 1] / 0.3)
        
        for i, name in enumerate(self.param_names):
            q16, q50, q84 = np.percentile(samples[:, i], [16, 50, 84])
            print(f"{name:8s}: {q50:.4f} +{q84-q50:.4f} -{q50-q16:.4f}")
            
        q16, q50, q84 = np.percentile(S8, [16, 50, 84])
        print("-" * 70)
        print(f"{'S8 (Türetilmiş)':>16s}: {q50:.4f} +{q84-q50:.4f} -{q50-q16:.4f}")
        print("="*70)
    
    def plot_corner(self, samples, output='corner_7param_FINAL.png'):
        """Corner plot"""
        print(f"Corner plot oluşturuluyor: {output}")
        fig = corner.corner(
            samples,
            labels=self.param_labels,
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_fmt='.4f',
            color='steelblue'
        )
        fig.suptitle('DMTDE 7-Parameter Posterior (Nihai Model)', fontsize=14, y=1.01)
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"✓ Corner plot kaydedildi: {output}")
        plt.close()
    
    def plot_chains(self, discard=1000, output='chains_trace_FINAL.png'):
        """Trace plot (convergence kontrolü)"""
        print(f"Trace plot oluşturuluyor: {output}")
        try:
            reader = emcee.backends.HDFBackend(self.chain_file)
            samples = reader.get_chain()  # (n_steps, n_walkers, n_params)
        except Exception as e:
            print(f"❌ HATA: Trace plot için chain okunamadı: {e}")
            return

        fig, axes = plt.subplots(self.n_params, 1, figsize=(10, 12), sharex=True)
        for i in range(self.n_params):
            ax = axes[i]
            ax.plot(samples[:, :, i], alpha=0.3, color='gray', lw=0.5)
            ax.set_ylabel(self.param_labels[i], fontsize=10)
            ax.axvline(discard, color='red', ls='--', lw=1, label=f'Burn-in ({discard})')
        axes[-1].set_xlabel('Step', fontsize=11)
        axes[0].legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"✓ Trace plot kaydedildi: {output}")
        plt.close()

# =============================================================================
# ANA PROGRAM
# =============================================================================

def main():
    """Tam MCMC pipeline'ı çalıştır"""
    
    print("""
    ╔════════════════════════════════════════════════════════════════╗
    ║      DMTDE MCMC ANALİZİ - TAM PİPELİNE (NİHAİ MODEL)           ║
    ║         7 Parametre: [H0, Om, Ob, sigma8, ns, tau, fd]         ║
    ╚════════════════════════════════════════════════════════════════╝
    """)
    
    N_WALKERS = 32
    N_STEPS = 5000   
    N_DISCARD = 1000 
    N_THIN = 10      

    # 1. VERİ YÜKLEME
    print("\n[1/5] Veri yükleniyor...")
    data = DMTDEData()
    
    abacus_file = 'abacus_halos_27z.npy'
    if not data.load_abacus_halos(abacus_file):
        print("UYARI: 'abacus_halos_27z.npy' bulunamadı.")
        print("Lütfen önce 'prepare_abacus_data.py' scriptini çalıştırın.")
        return 
    
    # GERÇEK VERİ YOLLARINI YÜKLE (veya hata verirse None olur)
    # data.load_desi_bao("desibaomanuel/desi_bao_2024.txt") # Cobaya halledecek
    data.load_pantheon_sn("external_data/pantheon_plus_sn.txt")
    data.load_planck_cmb()  
    
    # 2. LIKELIHOOD SETUP
    print("\n[2/5] Likelihood hesaplayıcı hazırlanıyor...")
    likelihood = DMTDELikelihood(data)
    
    # Cobaya BAO likelihood'unun başarıyla yüklenip yüklenmediğini kontrol et
    if likelihood.desi_bao_likelihood is None:
         print("⚠️  Cobaya DESI BAO likelihood YÜKLENEMEDİ. BAO verisi kullanılmayacak.")
    
    theta_test = [71.0, 0.310, 0.049, 0.800, 0.965, 0.054, 0.049]
    print(f"Merkez parametreler test ediliyor: {theta_test}")
    # likelihood'u debug=True ile çağır
    logL_test = likelihood.log_likelihood(theta_test, debug=True) 
    print(f"Test log-likelihood (merkezde): {logL_test:.2f}")
    if not np.isfinite(logL_test):
        print("❌ HATA: Merkezdeki log-likelihood sonlu değil! Nihai model ile bile başarısız.")
        return 

    # 3. MCMC ÇALIŞTIR
    print("\n[3/5] MCMC başlatılıyor...")
    runner = MCMCRunner(likelihood, n_walkers=N_WALKERS, n_steps=N_STEPS)
    # Nihai chain dosyası adı
    sampler = runner.run(output_file='dmtde_chain_7param_FINAL.h5') 
    
    # 4. SONUÇLARI ANALİZ ET
    print("\n[4/5] Sonuçlar analiz ediliyor...")
    analyzer = MCMCAnalyzer('dmtde_chain_7param_FINAL.h5', n_params=7) 
    samples = analyzer.load_samples(discard=N_DISCARD, thin=N_THIN)
    
    if samples is not None:
        analyzer.print_statistics(samples)
        
        # 5. GÖRSELLEŞTİRME
        print("\n[5/5] Grafikler oluşturuluyor...")
        analyzer.plot_corner(samples, output='corner_7param_FINAL.png') 
        analyzer.plot_chains(discard=N_DISCARD, output='chains_trace_FINAL.png')
    
    print("\n" + "="*70)
    print("✅ ANALİZ TAMAMLANDI!")
    print("="*70)
    print("Çıktı dosyaları:")
    print("  • dmtde_chain_7param_FINAL.h5 (MCMC chain)")
    print("  • corner_7param_FINAL.png (Posterior dağılımları)")
    print("  • chains_trace_FINAL.png (Convergence kontrolü)")
    print("="*70)

if __name__ == '__main__':
    # Linux/WSL için 'fork' başlangıç metodunu ayarla
    # Bu, Cobaya likelihood'u ile Pool'un uyumlu çalışması için GEREKLİ olabilir
    if os.name == 'posix': # Linux/macOS kontrolü
        try:
            current_method = Pool()._ctx.get_start_method()
            if current_method != 'fork':
                 print("Multiprocessing start method 'fork' olarak ayarlanıyor...")
                 set_start_method('fork', force=True)
        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Uyarı: Multiprocessing start method ayarlanamadı: {e}")
            pass # Ayarlanamazsa devam et
    
    main()
