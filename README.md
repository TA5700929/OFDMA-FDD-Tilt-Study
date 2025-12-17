# Downlink Throughput vs. Antenna Tilt (OFDMA/FDD)

This repository contains the MATLAB simulation data and figures for the IEEE Access submission  
**“The Effect of Sectorized Antenna Tilt on the Throughput and Range Characteristics of the OFDMA/FDD Radio Interface in the Downlink”**.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17619548.svg)](https://doi.org/10.5281/zenodo.17619548)

## Folders
- `data/` – four CSV files produced by the simulation.  
- `figures/` – 12 plots used in the paper.  
- `scripts/` – MATLAB code (run `lte_tilt_validation.m` to reproduce results).

## How to reproduce
1. Install MATLAB (R2020b or newer).  
2. Clone/download this repo.  
3. Open MATLAB, `cd` to the `scripts` folder.  
4. Run `lte_tilt_validation` (just type the name without `.m`).  
5. When it finishes, the same CSV files and figures will appear again.

## License
MIT – see LICENSE file.

## Citation
If you use the dataset, please cite:  
Tanvir Ahmed, “Simulation dataset for ‘The effect of sectorized antenna tilt on the throughput and range characteristics of the OFDMA/FDD radio interface in the downlink’,” Zenodo, 2025. https://doi.org/10.5281/zenodo.17619548
