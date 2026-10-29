# ReadMe

The code in this repository was used in the following study:

> Arai, T., and Kawamura, Y. *Data-driven Reconstruction of the Phase Dynamics for the El Niño–Southern Oscillation Reveals Seasonality and Statistical Synchronization with the Annual Cycle*.

The MATLAB analysis is adapted from the code for the following study in https://github.com/dg227/NLSA. 
The NLSA library is included in `nlsa/` with minor local modifications, including support for ERSSTv5.

The Koopman analysis of sea surface temperature (SST) in this repository is based on the original NLSA implementation accompanying the following study:

> Froyland, G., Giannakis, D., Lintner, B. R., Pike, M., and Slawinska, J. (2021). *Spectral analysis of climate dynamics with operator-theoretic approaches*. Nature Communications **12**, 6570. [doi:10.1038/s41467-021-26357-x](https://doi.org/10.1038/s41467-021-26357-x).

## Requirements

macOS/Linux (or WSL), Bash, Conda, and curl. 
Python dependencies are specified in [environment.yml](environment.yml).

Precomputed MATLAB results for the ERSSTv5 Koopman analysis are included in `SRC/enso_ERSSTv5_2026/myFolder/`. MATLAB is not required when using these files.

## Run

Run from the repository root:

```bash
conda env create -f environment.yml
bash SRC/download_sst_data.sh
bash SRC/enso_ERSSTv5_2026/run_all.sh
```

The download script obtains [monthly ERSSTv5 SST](https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc) from NOAA and saves it to:

```text
SRC/enso_ERSSTv5_2026/CLIMATE_DATA/sst.mnmean.nc
```

Manual downloads should use the same path. NOAA input data are not bundled.

Figures are saved in `SRC/figure/`.

License: [BSD 3-Clause](LICENSE). Third-party components retain their own licenses.
