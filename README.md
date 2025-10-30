# DMTDE: Dark Matter to Dark Energy Transition

**Paper:** [PDF (Zenodo)](https://doi.org/10.5281/zenodo.17469515)  
**Author:** Emre Özyurt  
**Email:** emre.ozyurt@proton.me  
**Date:** 2025-10-30

## Abstract
We introduce the **Dark Matter to Dark Energy Transition (DMTDE)** model... *(tam abstract yukarıda)*

## Repository Contents
- `main.tex`: LaTeX source
- `DMTDE_paper.pdf`: Manuscript
- `figures/`: All plots
- `supplementary.pdf`: Extended tables, derivations

## Reproducibility

## Pipeline Usage
1. `process_all_z.py`  
   - Input: AbacusSummit halo catalogs (`z*` folders)  
   - Output: `z_summary.json`, `all_z_results.txt`  
   - Computes mean halo mass and applies DMTDE suppression (D_Γ² = 0.975)

2. `plot_highz.py`  
   - Reads `z_summary.json`  
   - Generates `figures/highz_suppression.png`
