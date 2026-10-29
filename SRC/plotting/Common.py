#!/usr/bin/env python3
# coding: utf-8
"""Read ERSSTv5_2026 Koopman outputs and SST fields for paper figures."""

from __future__ import annotations

from plot_config import FIGURE_OUTPUT_DIR

from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Literal

import numpy as np
import scipy.io


DatasetKey = Literal["ERSSTv5_2026"]


_REPO_ROOT = Path(__file__).resolve().parents[2]
_DATA_ROOT = _REPO_ROOT / "SRC"
_TWOPI = 2.0 * np.pi

# _DATASETS[key] の tuple は次の順番:
# (dataset_folder, annual_index, annual_conj_index, enso_index, enso_conj_index, phi_a_shift_rad, phi_e_shift_rad)
# mode index は Python の 0-based 列番号。論文・notebook 上の g 番号は index + 1。
# phi_*_shift_rad は common_fixed.py と同じ位相原点合わせの回転量。
_DATASETS: dict[str, tuple[str, int, int, int, int, float, float]] = {
    "ERSSTv5_2026": ("enso_ERSSTv5_2026", 1, 2, 8, 9, 6.280577758709, 5.912686606382),
}

# _SST_DATASETS[key] の tuple は次の順番:
# (SST NetCDFの相対path, longitude変数名, latitude変数名, SST変数名)
_SST_DATASETS: dict[str, tuple[tuple[str, ...], str, str, str]] = {
    "ERSSTv5_2026": (("CLIMATE_DATA", "sst.mnmean.nc"), "lon", "lat", "sst"),
}


def _matlab_datenum_to_datetime(value: float) -> datetime:
    """MATLAB datenum 1個を Python の ``datetime`` に変換する。"""
    ordinal = int(value)
    return datetime.fromordinal(ordinal) + timedelta(days=float(value) - ordinal) - timedelta(days=366)


def _matlab_datenum_array_to_datetime(values: np.ndarray) -> np.ndarray:
    """MATLAB datenum 配列を同じ長さの ``datetime`` 配列に変換する。"""
    return np.array([_matlab_datenum_to_datetime(v) for v in np.asarray(values).ravel()])


def _load_time_series_struct(path: Path, variable: str) -> dict[str, np.ndarray]:
    """``Z.mat`` や ``Z_left.mat`` の MATLAB struct を Python 用の向きで読む。"""
    mat_data = scipy.io.loadmat(path, variable_names=(variable,))
    data = mat_data[variable]
    return {
        "time": _matlab_datenum_array_to_datetime(data["time"][0, 0]),
        "idx": data["idx"][0, 0].T,
        "angle": data["angle"][0, 0].T,
    }


def _load_nino34(path: Path) -> dict[str, np.ndarray]:
    """``Nino34_data.mat`` を読み、時刻・複素 Nino3.4 指標・角度を返す。"""
    mat_data = scipy.io.loadmat(path, variable_names=("Nino34",))
    nino34 = mat_data["Nino34"]
    return {
        "time": _matlab_datenum_array_to_datetime(nino34["time"][0, 0]),
        "idx": np.real_if_close(nino34["idx"][0, 0].ravel()),
        "angle": nino34["angle"][0, 0].ravel(),
    }


def _load_lambda(path: Path) -> dict[str, np.ndarray]:
    """``Lambda.mat`` から Koopman 固有値と ``gamma = lambda * 12 / (2*pi)`` を読む。"""
    mat_data = scipy.io.loadmat(path, variable_names=("Lambda",))
    lambda_data = mat_data["Lambda"]
    return {
        "eigenvalue": lambda_data["value"][0, 0].ravel(),
        "gamma": lambda_data["gamma"][0, 0].ravel(),
    }


def _load_mu(myfolder: Path) -> np.ndarray | None:
    """``mu.mat`` が存在すれば重みベクトルを返し、なければ ``None`` を返す。"""
    path = myfolder / "mu.mat"
    if not path.exists():
        return None
    return scipy.io.loadmat(path, variable_names=("mu",))["mu"].ravel()


def _optional_z_left_path(myfolder: Path) -> Path | None:
    """``myFolder`` 内の left eigenfunction ファイル名を吸収して探す。"""
    for name in ("Z_left.mat", "Zleft.mat"):
        path = myfolder / name
        if path.exists():
            return path
    return None


def _netcdf_time_to_datetime(time_var) -> np.ndarray:
    """NetCDF の time 座標を Python の ``datetime`` 配列に変換する。"""
    import netCDF4 as nc4

    dates = nc4.num2date(
        time_var[:],
        units=time_var.units,
        calendar=getattr(time_var, "calendar", "standard"),
        only_use_cftime_datetimes=False,
        only_use_python_datetimes=True,
    )
    return np.array([datetime(d.year, d.month, d.day) for d in dates])


def _load_sst_aligned_to_z(key: DatasetKey, dataset_dir: Path, z_dates: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray, np.ma.MaskedArray]:
    """
    SST NetCDF を ``Z.time`` の年月に合わせて読み、ENSO用・annual用のSST場を返す。

    返り値は ``lon, lat, sst_monthly_anomaly_celsius, sst_zscore``。
    ``sst_monthly_anomaly_celsius`` は同じ暦月の平均を引いた monthly anomaly [C]で、
    ENSO composite に使う。
    ``sst_zscore`` は monthly anomaly ではない SST field を各格子点で全期間標準化した場で、
    annual composite に使う。
    """
    import netCDF4 as nc4

    path_parts, lon_name, lat_name, sst_name = _SST_DATASETS[key]
    path = dataset_dir.joinpath(*path_parts)
    with nc4.Dataset(path) as ds:
        dates = _netcdf_time_to_datetime(ds["time"])
        by_month = {(d.year, d.month): i for i, d in enumerate(dates)}
        try:
            #* SST場は Z.time と同じ年月だけを同じ順番で取り出す。
            use_time = np.array([by_month[(d.year, d.month)] for d in z_dates])
        except KeyError as exc:
            raise ValueError(f"SST NetCDF does not contain {exc.args[0]} needed by Z.time.") from exc

        lon = np.asarray(ds[lon_name][:], dtype=float)
        lat = np.asarray(ds[lat_name][:], dtype=float)
        sst = ds[sst_name][use_time, :, :]

    lon_order = np.argsort(np.mod(lon, 360.0))
    lat_order = np.argsort(lat)
    lon = np.mod(lon[lon_order], 360.0)
    lat = lat[lat_order]
    sst = sst[:, lat_order, :][:, :, lon_order]

    use_lat = (-75.0 <= lat) & (lat <= 75.0)
    lat = lat[use_lat]
    sst = np.ma.asarray(sst[:, use_lat, :], dtype=float)
    sst = np.ma.masked_invalid(sst)
    sst = np.ma.masked_where((sst < -100.0) | (100.0 < sst), sst)

    sst_monthly_anomaly_celsius = np.ma.empty(sst.shape, dtype=float)
    months = np.array([d.month for d in z_dates])
    for month in range(1, 13):
        use_month = months == month
        block = sst[use_month]
        sst_monthly_anomaly_celsius[use_month] = block - np.ma.mean(block, axis=0)

    sst_mean = np.ma.mean(sst, axis=0)
    sst_std = np.ma.masked_less_equal(np.ma.std(sst, axis=0, ddof=1), 0.0)
    sst_zscore = (sst - sst_mean) / sst_std

    return lon, lat, sst_monthly_anomaly_celsius, sst_zscore


def _phase_from_complex(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """複素時系列から unwrap した位相と振幅を返す。"""
    z = np.asarray(z)
    phase = np.unwrap(np.mod(np.angle(z), _TWOPI))
    amplitude = np.abs(z)
    return phase, amplitude


def _amplitude_filter_mask(amplitude: np.ndarray, amp: float) -> tuple[np.ndarray, float]:
    """既存スクリプトの AMP と同じ percentile rank で振幅 filter を作る。"""
    percentile_rank = float(amp) * 100.0 if float(amp) <= 1.0 else float(amp)
    if not 0.0 <= percentile_rank <= 100.0:
        raise ValueError("amp must be in [0, 1] or [0, 100].")
    threshold = np.percentile(np.asarray(amplitude, dtype=float), percentile_rank)
    return np.asarray(amplitude, dtype=float) >= threshold, percentile_rank


def _phase_mask(phase: np.ndarray, phase_index: int, S: int) -> np.ndarray:
    """位相 ``phase_index`` の wedge に入るサンプルを選ぶ。"""
    width = _TWOPI / S
    center = width * phase_index
    return np.abs(np.mod(phase - center + np.pi, _TWOPI) - np.pi) <= 0.5 * width



#* -------------------------


def get_phase_series(data: dict[str, object]) -> dict[str, np.ndarray]:
    """
    ``get_data()`` の戻り値から annual/ENSO の位相と振幅を取り出す。

    元の各プロジェクトの ``get_Rangle()`` と同じ定義。
    位相は ``np.angle`` の ``[-pi, pi]`` を一度 ``[0, 2*pi)`` に直し、
    その後 ``np.unwrap`` で時間方向に連続な位相系列へ変換する。
    振幅は各 ``Z_*`` 複素時系列の絶対値 ``np.abs``。

    任意の複素配列ではなく、必ず ``data = get_data(...)`` で得た辞書を渡す。
    """
    if "annual_enso" not in data:
        raise KeyError("data must be the result of get_data(), containing 'annual_enso'.")

    annual_enso = data["annual_enso"]
    if not isinstance(annual_enso, dict):
        raise TypeError("data['annual_enso'] must be a dict.")

    required = ("date", "Z_annual", "Z_enso", "Z_annual_conj", "Z_enso_conj")
    missing = [key for key in required if key not in annual_enso]
    if missing:
        raise KeyError(f"data['annual_enso'] is missing keys: {missing}")

    #* 原点を合わせた固有関数(Z)の複素時系列から位相と振幅を計算する。
    phase_annual, amp_annual = _phase_from_complex(annual_enso["Z_annual"])
    phase_enso, amp_enso = _phase_from_complex(annual_enso["Z_enso"])
    phase_annual_conj, amp_annual_conj = _phase_from_complex(annual_enso["Z_annual_conj"])
    phase_enso_conj, amp_enso_conj = _phase_from_complex(annual_enso["Z_enso_conj"])

    return {
        #* 時刻
        "date": annual_enso["date"],
        #* 位相
        "phase_annual": phase_annual,
        "phase_enso": phase_enso,
        "phase_annual_conj": phase_annual_conj,
        "phase_enso_conj": phase_enso_conj,
        #* 振幅
        "amp_annual": amp_annual,
        "amp_enso": amp_enso,
        "amp_annual_conj": amp_annual_conj,
        "amp_enso_conj": amp_enso_conj,
    }



def get_data(key: DatasetKey = "ERSSTv5_2026") -> dict[str, object]:
    """
    指定データセットの ``myFolder`` から時系列・固有値・補助データを読む。

    返り値は2階層に分ける。
    ``annual_enso`` は annual/ENSO の選択 mode とその共役 pair だけを集めた辞書。
    ``all`` は ``Z``, ``Z_left``, ``Nino34``, ``Lambda``, ``mu`` など、読み込んだ全情報の辞書。

    ``date`` は MATLAB datenum から変換した ``datetime`` 配列。
    ``annual_enso`` の ``Z_annual`` と ``Z_enso`` は ``common_fixed.py`` と同じ
    位相原点合わせを適用した系列で、``raw_Z_*`` は ``Z.mat``/``Z_left.mat`` に保存された raw 系列。
    ``gamma`` は MATLAB 側の ``Lambda.value * 12 / (2 * pi)`` で、虚部を cycles/year として読める。
    """
    if key not in _DATASETS:
        valid = ", ".join(_DATASETS)
        raise ValueError(f"Unknown key {key!r}. Expected one of: {valid}")

    folder, annual_index, annual_conj_index, enso_index, enso_conj_index, phi_a_shift, phi_e_shift = _DATASETS[key]
    dataset_dir = _DATA_ROOT / folder
    myfolder = dataset_dir / "myFolder"

    z = _load_time_series_struct(myfolder / "Z.mat", "Z")
    nino34 = _load_nino34(myfolder / "Nino34_data.mat")
    lambda_dict = _load_lambda(myfolder / "Lambda.mat")
    mu = _load_mu(myfolder)
    if z["idx"].shape[0] != nino34["idx"].size:
        raise ValueError("Z.mat and Nino34_data.mat lengths differ.")

    if max(annual_index, annual_conj_index, enso_index, enso_conj_index) >= z["idx"].shape[1]:
        raise IndexError(f"Mode indices are out of bounds for Z.idx shape {z['idx'].shape}.")
    if max(annual_index, annual_conj_index, enso_index, enso_conj_index) >= lambda_dict["eigenvalue"].size:
        raise IndexError(f"Mode indices are out of bounds for Lambda length {lambda_dict['eigenvalue'].size}.")

    raw_z_annual = z["idx"][:, annual_index]
    raw_z_annual_conj = z["idx"][:, annual_conj_index]
    raw_z_enso = z["idx"][:, enso_index]
    raw_z_enso_conj = z["idx"][:, enso_conj_index]

    annual_factor = np.exp(-1j * phi_a_shift)
    enso_factor = np.exp(-1j * phi_e_shift)
    z_annual = raw_z_annual * annual_factor
    z_annual_conj = raw_z_annual_conj * np.conj(annual_factor)
    z_enso = raw_z_enso * enso_factor
    z_enso_conj = raw_z_enso_conj * np.conj(enso_factor)

    z_left_path = _optional_z_left_path(myfolder)
    z_left = None
    raw_z_annual_left = raw_z_annual_conj_left = raw_z_enso_left = raw_z_enso_conj_left = None
    if z_left_path is not None:
        z_left_variable = "Z_left" if z_left_path.name == "Z_left.mat" else "Zleft"
        z_left = _load_time_series_struct(z_left_path, z_left_variable)
        if max(annual_index, annual_conj_index, enso_index, enso_conj_index) >= z_left["idx"].shape[1]:
            raise IndexError(f"Mode indices are out of bounds for Z_left.idx shape {z_left['idx'].shape}.")
        raw_z_annual_left = z_left["idx"][:, annual_index]
        raw_z_annual_conj_left = z_left["idx"][:, annual_conj_index]
        raw_z_enso_left = z_left["idx"][:, enso_index]
        raw_z_enso_conj_left = z_left["idx"][:, enso_conj_index]

    mode_info = {
        "indices": {
            "annual": annual_index,
            "annual_conj": annual_conj_index,
            "enso": enso_index,
            "enso_conj": enso_conj_index,
        },
        "numbers": {
            "annual": annual_index + 1,
            "annual_conj": annual_conj_index + 1,
            "enso": enso_index + 1,
            "enso_conj": enso_conj_index + 1,
        },
    }
    phase_origin = {
        "phi_a_shift_rad": phi_a_shift,
        "phi_e_shift_rad": phi_e_shift,
        "mode_indices": {"annual": annual_index, "enso": enso_index},
        "conj_mode_indices": {"annual": annual_conj_index, "enso": enso_conj_index},
        "mode_numbers": {"annual": annual_index + 1, "enso": enso_index + 1},
        "conj_mode_numbers": {"annual": annual_conj_index + 1, "enso": enso_conj_index + 1},
        "phase_rotation_applied_in": "SRC.plotting.Common.get_data",
        "source": f"{folder}/SRC_py/common_fixed.py",
    }
    annual_enso = {
        #* 時刻
        "date": z["time"],
        #* nino3.4 指標
        "nino3.4": nino34["idx"],
        #* 固有関数(nino3.4を用いて0を合わせたもの / 右固有関数)
        "Z_annual": z_annual,
        "Z_enso": z_enso,
        #* 選んだ固有関数の共役(nino3.4を用いて0を合わせたもの / 右固有関数)
        "Z_annual_conj": z_annual_conj,
        "Z_enso_conj": z_enso_conj,
        #* rawは回転していない右固有関数
        "raw_Z_annual": raw_z_annual,
        "raw_Z_enso": raw_z_enso,
        "raw_Z_annual_conj": raw_z_annual_conj,
        "raw_Z_enso_conj": raw_z_enso_conj,
        #* rawは回転していない左固有関数
        "raw_Z_annual_left": raw_z_annual_left,
        "raw_Z_enso_left": raw_z_enso_left,
        "raw_Z_annual_conj_left": raw_z_annual_conj_left,
        "raw_Z_enso_conj_left": raw_z_enso_conj_left,
        #* 固有値
        "lambda_annual": lambda_dict["eigenvalue"][annual_index],
        "lambda_enso": lambda_dict["eigenvalue"][enso_index],
        "lambda_annual_conj": lambda_dict["eigenvalue"][annual_conj_index],
        "lambda_enso_conj": lambda_dict["eigenvalue"][enso_conj_index],
        #* gamma: imag(lambda)をcycles/yearで読めるようにした値 (lambda*12/(2*pi))
        "gamma_annual": lambda_dict["gamma"][annual_index],
        "gamma_enso": lambda_dict["gamma"][enso_index],
        "gamma_annual_conj": lambda_dict["gamma"][annual_conj_index],
        "gamma_enso_conj": lambda_dict["gamma"][enso_conj_index],
        #* mode番号と位相原点
        "mode_info": mode_info,
        "phase_origin": phase_origin,
    }

    all_data = {
        "key": key,
        "dataset_dir": dataset_dir,
        "myfolder": myfolder,
        "date": z["time"],
        "Z": z,
        "Z_left": z_left,
        "nino34": nino34,
        "Lambda": lambda_dict,
        "mu": mu,
        "mode_info": mode_info,
        "phase_origin": phase_origin,
    }
    return {
        "annual_enso": annual_enso,
        "all": all_data,
    }





def get_SST_phase_composite(key1: DatasetKey = "ERSSTv5_2026",
                            key2: Literal["ENSO", "Annual"] = "ENSO",
                            S: int = 8,
                            amp: float = 0.3) -> Dict[str, object]:
    """
    指定データセットの SST phase composite を返す。

    ``key2="ENSO"`` では monthly anomaly field [°C] を使う。
    ``key2="Annual"`` では monthly anomaly ではない SST field を各格子点で標準化した z-score を使う。
    振幅 filter は既存の phase composite スクリプトと同じく ENSO mode の振幅で決める。

    Returns:
        Dict lon, lat, composite=((1, X_1), (2, X_2), ..., (S, X_S)) where S is the number of phases.
        X_sは2d MaskedArrayで、axes[0]が緯度、axes[1]が経度のSST compositeを表す。
        単位は ``key2="ENSO"`` では[°C]、``key2="Annual"`` では[std.]。
        phase_indexは各時刻のphase番号で、振幅filter外の時刻は0、filter内は1..S。
        NetCDFの欠損値・陸域maskは既存の各プロジェクトと同じく ``np.ma`` として保持する。
        amp=0.3は既存コードの AMP=30 と同じで、振幅が30 percentile以上のデータを採用。
    """
    if key1 not in _DATASETS:
        valid = ", ".join(_DATASETS)
        raise ValueError(f"Unknown key1 {key1!r}. Expected one of: {valid}")
    if key2 not in ("ENSO", "Annual"):
        raise ValueError("key2 must be 'ENSO' or 'Annual'.")
    if S <= 0:
        raise ValueError("S must be positive.")

    folder = _DATASETS[key1][0]
    dataset_dir = _DATA_ROOT / folder
    data_all = get_data(key1)
    data = data_all["annual_enso"]
    phase_data = get_phase_series(data_all)
    lon, lat, sst_monthly_anomaly_celsius, sst_zscore = _load_sst_aligned_to_z(key1, dataset_dir, data["date"])

    #* ENSO mode の振幅 filter は既存 phase composite と同じく ENSO/Annual 両方で使う
    valid, amp_percentile = _amplitude_filter_mask(phase_data["amp_enso"], amp)

    if key2 == "ENSO":
        #* ENSO phase composite は ENSO mode の位相で分け、monthly anomaly field 側を平均する
        phase = phase_data["phase_enso"]
        sst_data = sst_monthly_anomaly_celsius
    elif key2 == "Annual":
        #* Annual phase composite は Annual mode の位相で分け、標準化した SST field 側を平均する
        phase = phase_data["phase_annual"]
        sst_data = sst_zscore
    else:
        raise ValueError("key2 must be 'ENSO' or 'Annual'.")

    phase_mod = np.mod(phase, _TWOPI)
    composites = []
    counts = np.empty(S, dtype=int)
    phase_index = np.zeros(phase_mod.shape, dtype=int)
    for i in range(S):
        mask = valid & _phase_mask(phase_mod, i, S)
        phase_index[mask] = i + 1
        counts[i] = int(np.count_nonzero(mask))
        composite = np.ma.mean(sst_data[mask], axis=0)
        composites.append((i + 1, np.ma.asarray(composite, dtype=float)))

    return {
        #* 経度・緯度
        "lon": lon,
        "lat": lat,
        #* phaseごとのSST composite: ((1, X_1), ..., (S, X_S))
        "composite": tuple(composites),
        #* 各phaseのcomposite計算に使ったサンプル数
        "counts": counts,
        #* 各時刻が属するphase番号。振幅filter外は0、filter内は1..S。
        "phase_index": phase_index,
        #* 入力指定
        "key1": key1,
        "key2": key2,
        "S": S,
        "amp": amp,
        #* composite fieldの単位
        "field_unit": "degC" if key2 == "ENSO" else "std.",
        #* amp filterで使ったpercentile rank
        "amp_percentile": amp_percentile,
    }
