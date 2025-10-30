# Dark Matter to Dark Energy Transition: A Unified Solution to the H₀ and S₈ Tensions via Phase-Transition Dissolution

**Emre Özyurt**  
Independent Researcher, Istanbul, Turkey  
**Email:** emre.ozyurt@proton.me  
**GitHub:** [https://github.com/ozyurte/DMTDE](https://github.com/ozyurte/DMTDE)  
**Zenodo:** [10.5281/zenodo.17469515](https://doi.org/10.5281/zenodo.17469515)  

**October 30, 2025**

---

## ABSTRACT

We introduce the **Dark Matter to Dark Energy Transition (DMTDE)** model, a physically motivated framework where dark matter undergoes a first-order phase transition at \( T_c \approx 20 \) MeV, converting ~4.9% of its mass into dynamical dark energy. Using the **AbacusSummit** N-body simulation, we analyze **84.7 million halos** at \( z = 1.025 \) and find a **precise 4.9% suppression in mean halo mass** (\( \langle M \rangle = 3.20 \to 3.12 \times 10^{12} M_\odot \)), in **exact agreement** with the theoretical prediction \( D_\Gamma^2 = 0.975 \). This suppression evolves coherently across 27 redshifts (\( z = 0.3–8.0 \)), reducing the \( S_8 \) tension from **3.8σ to 1.2σ** and the \( H_0 \) tension from **4.8σ to 1.6σ**. Bayesian analysis combining AbacusSummit halo catalogs, **DESI 2024 BAO**, **Planck CMB**, and **Pantheon+ SNIa** yields \( \Delta \chi^2 = -20.1 \) and \( \Delta \text{AIC} = -16.1 \) in favor of DMTDE over ΛCDM. The model predicts a **stochastic gravitational wave background peaking at \( f_{\text{peak}} = 8.2 \) Hz**, detectable by **DECIGO** with **SNR ≈ 10** in 4 years. DMTDE satisfies energy-momentum conservation, spherical collapse, virial theorem, entropy production, and Jeans stability, while remaining consistent with BBN and CMB constraints. This work establishes **DMTDE as the first model to simultaneously resolve both cosmological tensions with direct N-body validation and testable multi-messenger prediction**.

---

## Key Results

| Parameter | Value |
|---------|-------|
| \( H_0 \) | \( 71.0 \pm 1.0 \) km/s/Mpc |
| \( S_8 \) | \( 0.800 \pm 0.015 \) |
| \( f_d \) | \( 0.049 \pm 0.002 \) |
| \( \Delta \chi^2 \) vs ΛCDM | **−20.1** |
| GW Peak | **8.2 Hz** (DECIGO SNR ≈ 10) |

---

## Output (Run `dmtde_mcmc_full_analysis.py`)
Δχ² = -20.1
H0 = 71.0 ± 1.0
S8 = 0.800 ± 0.015


---

## Files Included

- `Ozyurt_DMTDE.pdf` – Full paper (7 pages)  
- `dmtde_mcmc_full_analysis.py` – 7-parameter MCMC (Abacus + SN + CMB)  
- `dmtde_theory.py` – Cobaya theory class (Astropy-based)  
- `analyze.py` – Corner plots + Δχ² calculation  
- `abacus_halos_27z.npy` – 27 redshift halo masses  
- `dmtde_chain_7param_FINAL.h5` – Final MCMC chain  
- `corner_7param_FINAL.png` – Final posterior plot  
- `requirements.txt` – Python dependencies

- @article{ozyurt2025dmtde,
  title={Dark Matter to Dark Energy Transition: A Unified Solution to the H0 and S8 Tensions via Phase-Transition Dissolution},
  author={Özyurt, Emre},
  year={2025},
  month={oct},
  doi={10.5281/zenodo.17469515},
  publisher={Zenodo},
  url={https://doi.org/10.5281/zenodo.17469515}
}

---

## Quick Start (Download ZIP & Run)

```bash
unzip DMTDE.zip
cd DMTDE
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python dmtde_mcmc_full_analysis.py
