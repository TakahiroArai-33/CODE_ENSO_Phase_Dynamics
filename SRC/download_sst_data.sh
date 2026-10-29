#!/usr/bin/env bash
# Download the NOAA PSL monthly ERSSTv5 input without replacing an existing file.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="${script_dir}/enso_ERSSTv5_2026/CLIMATE_DATA"
destination="${data_dir}/sst.mnmean.nc"
url="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
mkdir -p "${data_dir}"
if [[ -s "${destination}" ]]; then
    printf 'Already present: %s\n' "${destination}"
    exit 0
fi
curl --fail --location --retry 5 --retry-delay 5 \
    --output "${destination}.part" "${url}"
mv "${destination}.part" "${destination}"
printf 'Saved: %s\n' "${destination}"
