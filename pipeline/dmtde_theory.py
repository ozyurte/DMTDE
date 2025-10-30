# dmtde_theory.py
# Cobaya için özel DMTDE teorik modeli (Astropy kullanarak)

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from cobaya.theory import Theory # Cobaya'nın Theory sınıfını miras alıyoruz
import warnings
try:
    from astropy.cosmology import CosmologyWarning
    warnings.filterwarnings("ignore", category=CosmologyWarning)
except ImportError:
    pass

class DMTDETheory(Theory):
    """
    Cobaya için DMTDE teorik modeli.
    Astropy kullanarak doğru DMTDE H(z) ve mesafeleri hesaplar.
    Ayrıca BAO için rs_drag (rd) değerini hesaplar.
    
    Parametreler: H0, Omega_m, Omega_b, fd 
                  (sigma8, ns, tau bu sınıf tarafından kullanılmaz)
    """
    # Cobaya'ya bu teorinin hangi parametreleri kullandığını söyleyelim
    # (YAML'daki params ile eşleşmeli)
    params = {'H0': None, 'Omega_m': None, 'Omega_b': None, 'fd': None}

    def initialize(self):
        """Teori başlatıldığında bir kez çağrılır."""
        # Başlatma sırasında özel bir şey yapmamıza gerek yok şimdilik
        pass

    def get_requirements(self):
        """
        Bu teorinin hangi çıktıları üreteceğini ve 
        hangi parametrelere ihtiyaç duyduğunu tanımlar.
        Likelihood'lar bu isimleri kullanarak değerleri ister.
        """
        return {
            # Fonksiyonlar: Likelihood'lar z değeri vererek çağırabilir
            'Hubble': {'z': np.linspace(0, 5, 100)}, # Örnek z aralığı, likelihoodlar kendi z'lerini kullanır
            'comoving_radial_distance': {'z': np.linspace(0, 5, 100)},
            'luminosity_distance': {'z': np.linspace(0, 5, 100)},
            'angular_diameter_distance': {'z': np.linspace(0, 5, 100)},
            'distmod': {'z': np.linspace(0, 5, 100)}, # Mesafe modülü
            # Değerler: Tek bir sayısal değer
            'rs_drag': None # BAO için sound horizon at drag epoch
            # Not: Cobaya'da rd yerine genelde rs_drag kullanılır
        }

    def must_provide(self, **requirements):
        """
        Likelihood'ların istediği 'requirements'leri hesaplamak için çağrılır.
        """
        # Gelen parametreleri alalım
        H0 = self.provider.get_param('H0')
        Om = self.provider.get_param('Omega_m')
        Ob = self.provider.get_param('Omega_b')
        fd = self.provider.get_param('fd')

        # Astropy kozmoloji nesnesini oluştur
        Om_eff = Om * (1.0 - fd)
        Ob_use = Ob # Baryon yoğunluğunun değişmediğini varsayıyoruz

        # Geçerlilik kontrolü
        if Om_eff < 0 or (1.0 - Om_eff) < 0 or Ob_use < 0 or H0 < 10 or Ob_use > Om_eff:
            # Geçersiz parametreler, Cobaya'ya likelihood'un -inf olması gerektiğini söyle
            self.provider.set_aborted() # Bu adımı geçersiz say
            return False # Hesaplama yapılamadı

        try:
            # Ana DMTDE kozmolojisi (H(z) için)
            cosmo_dmtde = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om_eff, Ob0=Ob_use)
            # rs_drag hesaplaması için LCDM kozmolojisi (erken evren fd=0)
            cosmo_lcdm = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om, Ob0=Ob)
        except ValueError:
            self.provider.set_aborted() 
            return False 

        # İstenen hesaplamaları yap ve state sözlüğüne ekle
        state = {}
        # rs_drag'ı hesapla (LCDM kullanarak)
        try:
            # Astropy'nin rd hesaplama yöntemi
             if hasattr(cosmo_lcdm, 'sound_horizon_drag'): # Eski versiyonlar
                 rd_value = cosmo_lcdm.sound_horizon_drag.to(u.Mpc).value
             elif hasattr(cosmo_lcdm, 'sound_horizon_z'): # Yeni versiyonlar
                 rd_value = cosmo_lcdm.sound_horizon_z(cosmo_lcdm.zs_drag).to(u.Mpc).value
             else: # Yaklaşık formül (kullanılmamalı)
                 h = H0 / 100.0; Obh2 = Ob * h**2; Omh2 = Om * h**2
                 rd_value = 55.154 * np.exp(-72.3 * (Omh2 * Obh2)**0.28) / (Omh2**0.25 * Obh2**0.078) * h**-1
                 
            if not np.isfinite(rd_value): raise ValueError("rs_drag sonlu değil")
            state['rs_drag'] = rd_value
        except Exception:
             self.provider.set_aborted()
             return False

        # Fonksiyonları tanımla (DMTDE kozmolojisini kullanarak)
        # Bu fonksiyonlar z değeri alıp ilgili miktarı hesaplamalı
        # Hata kontrolü ekleyelim
        get_val = lambda func, z: func(z).value if np.all(np.isfinite(func(z).value)) else np.inf
        
        state['Hubble'] = \
            lambda z: cosmo_dmtde.H(z).to(u.km / u.s / u.Mpc).value
        state['comoving_radial_distance'] = \
            lambda z: cosmo_dmtde.comoving_distance(z).to(u.Mpc).value
        state['luminosity_distance'] = \
            lambda z: cosmo_dmtde.luminosity_distance(z).to(u.Mpc).value
        state['angular_diameter_distance'] = \
            lambda z: cosmo_dmtde.angular_diameter_distance(z).to(u.Mpc).value
        state['distmod'] = \
            lambda z: cosmo_dmtde.distmod(z).value

        # Hesaplanan değerleri Cobaya'ya bildir
        self.provider.set_current_results(state)
        return True # Hesaplama başarılı

    def get_Hubble(self, z):
        """H(z) değerini döndürür (km/s/Mpc)"""
        # Bu fonksiyon aslında doğrudan çağrılmaz, state içindeki lambda kullanılır
        # Ama Cobaya bazen bu isimde fonksiyonlar arayabilir
        return self._current_state['Hubble'](z)

    def get_comoving_radial_distance(self, z):
        """Comoving distance D_M(z) değerini döndürür (Mpc)"""
        return self._current_state['comoving_radial_distance'](z)
        
    def get_luminosity_distance(self, z):
        """Luminosity distance D_L(z) değerini döndürür (Mpc)"""
        return self._current_state['luminosity_distance'](z)

    def get_angular_diameter_distance(self, z):
        """Angular diameter distance D_A(z) değerini döndürür (Mpc)"""
        return self._current_state['angular_diameter_distance'](z)
        
    def get_distmod(self, z):
        """Distance modulus mu(z) değerini döndürür"""
        # Güvenlik: log10(0) hatasını önlemek için
        dL = self.get_luminosity_distance(z)
        valid = np.isfinite(dL) & (dL > 1e-9)
        mu = np.full_like(np.atleast_1d(z), np.inf, dtype=float)
        if np.any(valid):
             mu[valid] = 5 * np.log10(dL[valid]) + 25
        return mu[0] if mu.shape==(1,) else mu # Skalerse skaler döndür

    def get_rs_drag(self):
        """Sound horizon rs_drag değerini döndürür (Mpc)"""
        return self._current_state['rs_drag']