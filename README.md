# Dark Matter to Dark Energy Transition: A Unified Solution to the H₀ and S₈ Tensions via Phase-Transition Dissolution

**Emre Özyurt**  
Independent Researcher, Istanbul, Turkey  
**Email:** emre.ozyurt@proton.me  
**GitHub:** [https://github.com/ozyurte/DMTDE](https://github.com/ozyurte/DMTDE)  
**Zenodo:** [10.5281/zenodo.17469515](https://doi.org/10.5281/zenodo.17469515)  

**October 31, 2025 (v2.0)**

---

## ABSTRACT

We introduce the **Dark Matter to Dark Energy Transition (DMTDE)** model, a physically motivated framework where dark matter undergoes a first-order phase transition at Tc ≈ 20 MeV, converting ~4.9% of its mass into dynamical dark energy. Using a **semi-analytical approach validated against the AbacusSummit** N-body simulation suite, we analyze **84.7 million halos** at z = 1.025 and find a **precise 4.9% suppression in mean halo mass** (⟨M⟩ = 3.20 → 3.12 × 10¹² M☉), in **exact agreement** with the theoretical prediction D²ᵧ = 0.975. This suppression evolves coherently across 27 redshifts (z = 0.3–8.0), reducing the S₈ tension from **3.8σ to 1.2σ** and the H₀ tension from **4.8σ to 1.6σ**. Bayesian analysis combining AbacusSummit halo catalogs, **DESI 2024 BAO**, **Planck CMB**, and **Pantheon+ SNIa** yields Δχ² = -20.1 and ΔAIC = -16.1 in favor of DMTDE over ΛCDM. The model predicts a **stochastic gravitational wave background peaking at f_peak = 8.2 Hz**, detectable by **DECIGO** with **SNR ≈ 10** in 4 years. DMTDE satisfies energy-momentum conservation, spherical collapse, virial theorem, entropy production, and Jeans stability, while remaining consistent with BBN and CMB constraints. This work establishes **DMTDE as the first model to simultaneously resolve both cosmological tensions through large-scale semi-analytical validation and testable multi-messenger prediction**.

---

## Key Results

| Parameter | Value |
|-----------|-------|
| H₀ | 71.5 ± 0.8 km/s/Mpc |
| S₈ | 0.792 ± 0.010 |
| f_d | 0.049 ± 0.002 |
| Δχ² vs ΛCDM | **−20.1** |
| GW Peak | **8.2 Hz** (DECIGO SNR ≈ 10) |

---

## Output (Run `pipeline/dmtde_mcmc_full_analysis.py`)
Δχ² = -20.1 H0 = 71.5 ± 0.8 km/s/Mpc S8 = 0.792 ± 0.010 fd = 0.049 ± 0.002
---

## Version History

### v2.0 (October 31, 2025) - **CURRENT**
- **Methodology clarification:** Explicitly described semi-analytical approach using AbacusSummit ΛCDM simulations
- **Technical fixes:** Resolved all LaTeX compilation issues (eliminated question marks in references)
- **Enhanced transparency:** Improved documentation for independent researcher constraints
- **New content:** Added Figure 3 (cosmological tension reduction visualization)
- **BibTeX corrections:** Fixed Pantheon+ [12], barreira2014, caprini2020 entries
- **Scientific results:** Unchanged from v1.0 (all numerical findings identical)

### v1.0 (October 29, 2025)
- Initial submission to Physical Review D, SCPMA, CPL, Nature
- Established DMTDE framework with AbacusSummit validation
- Demonstrated H₀ and S₈ tension resolution to <2σ

---

## Repository Structure

| Path | Description |
|------|-------------|
| **📄 Main Files** | |
| `Ozyurt_DMTDE.pdf` | Main paper (7 pages, v2.0) |
| `Ozyurt_DMTDE.tex` | LaTeX source |
| `references.bib` | BibTeX bibliography |
| **📁 pipeline/** | Analysis scripts |
| `├─ dmtde_mcmc_full_analysis.py` | Main MCMC analysis |
| `├─ dmtde_theory.py` | Cobaya theory class |
| `├─ prepare_abacus_data.py` | Data preprocessing |
| `├─ run_halo_suppression.py` | Suppression analysis |
| `└─ requirements.txt` | Python dependencies |
| **📁 data/** | Input datasets |
| `├─ abacus_halos_27z.npy` | 27 redshift halo masses |
| `├─ abacus_halos_27z.csv` | Same in CSV format |
| `├─ pantheon_plus_sn.txt` | Pantheon+ SNIa data |
| `└─ planck_2018_baseline.txt` | Planck CMB data |
| **📁 chains/** | MCMC outputs |
| `├─ dmtde_final.1.txt` | Final chain |
| `├─ dmtde_final.covmat` | Covariance matrix |
| `└─ dmtde_final.input.yaml` | Configuration |
| **📁 figures/** | All figures |
| `├─ dmtde_corner.png` | Posterior distributions |
| `├─ dmtde_mass_evolution.png` | Halo mass evolution |
| `├─ dmtde_hmf_ratio.png` | HMF ratio plot |
| `├─ dmtde_tension_reduction.pdf` | H₀/S₈ tensions (Fig 3) |
| `└─ dmtde_gw_spectrum.png` | GW prediction |
| **📁 supplementary/** | Supplementary materials |
| `└─ supplementary.pdf` | Extended derivations |

---

## Citation

```bibtex
@article{ozyurt2025dmtde,
  title={Dark Matter to Dark Energy Transition: A Unified Solution to the H0 and S8 Tensions via Phase-Transition Dissolution},
  author={Özyurt, Emre},
  year={2025},
  month={oct},
  doi={10.5281/zenodo.17469515},
  publisher={Zenodo},
  url={https://doi.org/10.5281/zenodo.17469515}
}

Data Availability
AbacusSummit simulations: https://abacussummit.readthedocs.io
Pantheon+ SNIa: https://pantheonplussh0es.github.io
Planck 2018: https://pla.esac.esa.int
DESI 2024 BAO: https://data.desi.lbl.gov
Contact
For questions or collaborations:
Emre Özyurt | emre.ozyurt@proton.me | Istanbul, Turkey

