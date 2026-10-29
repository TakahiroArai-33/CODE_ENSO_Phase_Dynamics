#!/usr/bin/env python3
# coding: utf-8

"""Shared helpers for the corrected ENSO phase analysis.

修正点:
- MATLAB datenum は ``datetime.fromordinal(d) - timedelta(days=366)`` で読む。
- ``Z.idx`` は MATLAB の ``MyCodeOutputData.m`` で保存した生の Koopman eigenfunction。
- annual 位相は ``PHI_a_SHIFT`` で 1月の circular mean が 0 になるように合わせる。
- ENSO 位相は ``PHI_e_SHIFT`` で Nino3.4 ピーク平均位相が 0 になるように合わせる。
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import netCDF4 as nc4
import numpy as np
import scipy.io


SCRIPT_DIR = Path(__file__).resolve().parent
ENSO_DIR = SCRIPT_DIR.parent
TWOPI = 2.0 * np.pi
IDX_G2 = 1
# 1970-2025 の ERSSTv5 では ENSO 周波数の mode が g9 に移る。
IDX_ENSO = 8
NINO34_PEAK_PERCENTILE = 90.0

#* 確認用_位相回転量_g2_g9.ipynb で求めた既定の回転量。
#* get_data() では g * exp(-i * PHI_*_SHIFT) として適用する。
#* PHI_a_SHIFT: g2 の 1月 circular mean が phi_a=0 になる回転量。
#* PHI_e_SHIFT: g9 の Nino3.4 上位ピーク circular mean が phi_e=0 になる回転量。
PHI_a_SHIFT = 6.280577758709
PHI_e_SHIFT = 5.912686606382


def matlab_datenum_to_datetime(value: float) -> datetime:
    """
    MATLAB の datenum を Python の ``datetime`` に変換する。

    MATLAB datenum は MATLAB 独自の「1 日ごとの通し番号」で、Python の
    ``datetime`` が使う ordinal と原点が違う。そのため単純に
    ``datetime(1, 1, 1) + timedelta(days=value)`` とすると、例えば
    MATLAB datenum 720319 が本来の 1972-03-01 ではなく 1973-03-03 として
    読まれてしまう。

    ここでは整数部を ``datetime.fromordinal`` で日付に直し、最後に
    366 日を引く。これで MATLAB の ``datestr`` が表示する年月日と同じ
    Python 日付になる。小数部 ``fraction`` は、将来 datenum に時刻成分が
    含まれていた場合に、その時刻を残すために足している。
    """
    ordinal = int(value)
    fraction = float(value) - ordinal
    return datetime.fromordinal(ordinal) + timedelta(days=fraction) - timedelta(days=366)


def matlab_datenum_array_to_datetime(values) -> np.ndarray:
    """MATLAB datenum の配列を、同じ長さの Python 日付配列に変換する。"""
    return np.array([matlab_datenum_to_datetime(v) for v in np.asarray(values).ravel()])


def phase_from_complex(z: np.ndarray) -> np.ndarray:
    """
    複素数時系列から位相を取り出し、時間方向に連続になるよう unwrap する。

    ``np.angle`` は ``[-pi, pi]`` の範囲に戻るため、そのまま微分すると
    ``pi`` 付近で大きなジャンプが出る。ここでは一度 ``[0, 2*pi)`` に
    直したあと ``np.unwrap`` し、位相速度を計算しやすい連続位相にする。
    """
    return np.unwrap(np.mod(np.angle(z), TWOPI))


def get_Rangle(g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    複素時系列 ``g`` から、位相 ``phi`` と振幅 ``amp`` の時系列を取り出す。

    研究上は、この関数を annual mode と ENSO mode の両方に共通して使う。
    ``phi`` は時間方向に unwrap した連続位相、``amp`` は複素数の絶対値である。
    """
    g = np.asarray(g)
    return phase_from_complex(g), np.abs(g)


def circular_mean_phase(phase: np.ndarray) -> float:
    """
    位相の平均を円周上で計算し、結果を ``[0, 2*pi)`` に戻す。

    通常の平均では 0 と 2*pi 付近の値が離れて見えてしまう。ここでは
    ``exp(i*phase)`` の平均角を使うので、位相として自然な平均になる。
    """
    phase = np.mod(np.asarray(phase, dtype=float), TWOPI)
    return float(np.mod(np.angle(np.mean(np.exp(1j * phase))), TWOPI))


def zero_phase_error(phase: np.ndarray) -> float:
    """Return the smaller circular distance between the mean phase and zero."""
    phase0 = circular_mean_phase(phase)
    return min(phase0, TWOPI - phase0)


def align_annual_to_january(phase: np.ndarray, dates: np.ndarray) -> tuple[np.ndarray, float]:
    """
    annual 位相を回転し、1 月の平均位相が 0 になるように合わせる。

    annual 位相 ``phi_a`` は季節との対応を読みやすくするため、1 月を
    基準位相 0 としてそろえる。
    """
    january = np.array([d.month == 1 for d in dates])
    shift = circular_mean_phase(phase[january])
    aligned = phase - shift
    err = zero_phase_error(aligned[january])
    print(f"phi_a origin: January circular mean shift={shift:.6f} rad, zero error={err:.3e} rad")
    return aligned, shift


def positive_peak_mask(values: np.ndarray, percentile: float = NINO34_PEAK_PERCENTILE) -> tuple[np.ndarray, float]:
    """
    Nino3.4 などの時系列から、上位 percentile 以上の局所ピークを選ぶ。

    ピーク位相を確認するための補助関数。条件に合う局所ピークが一つも
    見つからない場合は、全期間で最大の点を一つだけピークとして返す。
    """
    values = np.asarray(values, dtype=float).reshape(-1)
    local_peak = np.r_[False, (values[1:-1] >= values[:-2]) & (values[1:-1] >= values[2:]), False]
    threshold = float(np.nanpercentile(values, percentile))
    peak = local_peak & (values >= threshold)
    if not np.any(peak):
        peak = np.zeros(values.size, dtype=bool)
        peak[np.nanargmax(values)] = True
    return peak, threshold


def load_z(enso_dir: Path = ENSO_DIR) -> dict[str, np.ndarray]:
    """
    MATLAB で保存された ``Z.mat`` を読み、Python で扱いやすい形に直す。

    ``Z.idx`` は MATLAB 側で保存された生の Koopman eigenfunction である。
    配列は ``time x mode`` の向きに転置して返すので、
    この 2026 用ディレクトリでは既定で ``idx[:, 8]``、
    つまり g9 を ENSO mode の時系列として使う。
    """
    mat_data = scipy.io.loadmat(enso_dir / "myFolder" / "Z.mat")
    Z = mat_data["Z"]
    return {
        "time": matlab_datenum_array_to_datetime(Z["time"][0, 0]),
        "idx": Z["idx"][0, 0].T,
        "angle": Z["angle"][0, 0].T,
    }


def load_lambda(enso_dir: Path = ENSO_DIR) -> dict[str, np.ndarray]:
    """Load MATLAB ``Lambda.mat`` eigenvalues and Koopman frequencies."""
    mat_data = scipy.io.loadmat(enso_dir / "myFolder" / "Lambda.mat")
    Lambda = mat_data["Lambda"]
    return {
        "eigenvalue": Lambda["value"][0, 0].ravel(),
        "gamma": Lambda["gamma"][0, 0].ravel(),
    }


def load_nino34_mat(enso_dir: Path = ENSO_DIR) -> dict[str, np.ndarray]:
    """Load MATLAB ``Nino34_data.mat`` with corrected datenum conversion."""
    mat_data = scipy.io.loadmat(enso_dir / "myFolder" / "Nino34_data.mat")
    Nino34 = mat_data["Nino34"]
    return {
        "time": matlab_datenum_array_to_datetime(Nino34["time"][0, 0]),
        "idx": Nino34["idx"][0, 0].ravel(),
        "angle": Nino34["angle"][0, 0].ravel(),
    }


def get_data(
    enso_dir: Path = ENSO_DIR,
    nino34_peak_percentile: float = NINO34_PEAK_PERCENTILE,
    phi_a_shift_rad: float | None = PHI_a_SHIFT,
    phi_e_shift_rad: float = PHI_e_SHIFT,
    annual_index: int = IDX_G2,
    enso_index: int = IDX_ENSO,
) -> dict[str, object]:
    """
    修正版で使う annual, ENSO, date 時系列を一箇所で取り出す。

    ここが複素時系列の入口である。デフォルトでは annual mode として
    ``Z.idx[:, 1]``、ENSO mode として ``Z.idx[:, 8]`` を使う。将来、
    別のデータや別の mode を使う場合は ``annual_index`` と ``enso_index`` を
    変える。

    ``annual`` は 1 月の平均位相が 0 になるように回転して返す。
    ``enso`` は Nino3.4 ピーク平均位相が 0 になるように
    ``PHI_e_SHIFT`` で回転して返す。

    返り値の主なキー:
    - ``annual``: annual の補正済み複素時系列
    - ``enso``: ENSO の補正済み複素時系列
    - ``date``: MATLAB datenum を修正変換した日付
    - ``raw_annual``, ``raw_enso``: 回転前の ``Z.idx`` 複素時系列
    - ``phase_origin``: どの回転を掛けたかを保存するメタデータ
    """
    z = load_z(enso_dir)
    nino34 = load_nino34_mat(enso_dir)

    if z["idx"].shape[0] != nino34["idx"].size:
        raise ValueError("Z.mat and Nino34_data.mat lengths differ.")

    date = z["time"]
    raw_annual = z["idx"][:, annual_index]
    raw_enso = z["idx"][:, enso_index]

    raw_phi_a, _ = get_Rangle(raw_annual)


    #* annual 位相は 1 月の circular mean が 0 になるように回転する。
    if phi_a_shift_rad is None:
        january = np.array([d.month == 1 for d in date])
        phi_a_shift_rad = circular_mean_phase(raw_phi_a[january])

    #* phi_a_shift_rad, phi_e_shift_radだけ回転させる。
    annual = raw_annual * np.exp(-1j * float(phi_a_shift_rad))
    enso = raw_enso * np.exp(-1j * float(phi_e_shift_rad))

    phi_a, _ = get_Rangle(annual)
    january = np.array([d.month == 1 for d in date])
    phi_a_zero_error = zero_phase_error(phi_a[january])
    print(
        "annual origin: January circular mean "
        f"shift={float(phi_a_shift_rad):.6f} rad, zero error={phi_a_zero_error:.3e} rad"
    )

    phi_e, _ = get_Rangle(enso)
    nino_peak_mask, nino_peak_threshold = positive_peak_mask(
        np.real(nino34["idx"]),
        percentile=nino34_peak_percentile,
    )
    peak_phase = np.mod(phi_e[nino_peak_mask], TWOPI)
    peak_mean_phase = circular_mean_phase(peak_phase)
    print(
        "ENSO origin: using raw Z.idx plus Python shift; "
        f"ENSO shift={float(phi_e_shift_rad):.6f} rad; "
        f"Nino34_data peak mean phase={peak_mean_phase:.6f} rad, "
        f"threshold={nino_peak_threshold:.3f}, n={nino_peak_mask.sum()}"
    )

    phase_origin = {
        "phi_a_shift_rad": float(phi_a_shift_rad),
        "phi_e_shift_rad": float(phi_e_shift_rad),
        "enso_peak_percentile": nino34_peak_percentile,
        "enso_peak_mask": nino_peak_mask,
        "enso_peak_threshold": nino_peak_threshold,
        "enso_peak_mean_phase_rad": peak_mean_phase,
        "mode_indices": {"annual": int(annual_index), "enso": int(enso_index)},
        "phase_convention": "matlab_datenum_fixed; get_data returns January-aligned annual mode and Nino3.4-peak-aligned ENSO mode by default",
    }
    return {
        "annual": annual,
        "enso": enso,
        "date": date,
        "raw_annual": raw_annual,
        "raw_enso": raw_enso,
        "z": z,
        "nino34": nino34,
        "phase_origin": phase_origin,
    }


def netcdf_time_to_datetime(time_var) -> np.ndarray:
    """Convert a NetCDF time coordinate to plain Python ``datetime`` objects."""
    dates = nc4.num2date(
        time_var[:],
        units=time_var.units,
        calendar=getattr(time_var, "calendar", "standard"),
        only_use_cftime_datetimes=False,
        only_use_python_datetimes=True,
    )
    return np.array([datetime(d.year, d.month, d.day) for d in dates])


def load_sst_aligned_to_z(z_dates: np.ndarray, enso_dir: Path = ENSO_DIR):
    """
    ERSSTv5 の SST 場を、修正済み ``Z.time`` と同じ年月だけ取り出す。

    位相 composite で日付がずれないように、``z_dates`` の各年月に一致する
    SST だけを使う。返す ``sst_anom`` は月別 climatology からの標準化偏差、
    ``sst_std`` は全期間平均・標準偏差で標準化した場である。
    """
    path = enso_dir / "CLIMATE_DATA" / "sst.mnmean.nc"
    with nc4.Dataset(path) as ds:
        lat = ds["lat"][:]
        lon = ds["lon"][:]
        dates = netcdf_time_to_datetime(ds["time"])
        sst = ds["sst"][:]

    order = np.argsort(np.mod(lon, 360.0))
    lon = np.mod(lon[order], 360.0)
    sst = sst[:, :, order]

    by_month = {(d.year, d.month): i for i, d in enumerate(dates)}
    use_time = np.array([by_month[(d.year, d.month)] for d in z_dates])
    use_lat = (-75.0 <= lat) & (lat <= 75.0)
    lat = lat[use_lat]
    sst = sst[use_time][:, use_lat, :]

    sst_anom = np.ma.empty_like(sst)
    months = np.array([d.month for d in z_dates])
    for month in range(1, 13):
        use_month = months == month
        block = sst[use_month]
        sst_anom[use_month] = (block - np.ma.mean(block, axis=0)) / np.ma.std(block, axis=0)

    sst_std = (sst - np.ma.mean(sst, axis=0)) / np.ma.std(sst, axis=0, ddof=1)
    return lon, lat, sst, sst_anom, sst_std


def prepare_phase_data(
    enso_dir: Path = ENSO_DIR,
    nino34_peak_percentile: float = NINO34_PEAK_PERCENTILE,
) -> dict[str, object]:
    """
    旧形式の ``phase_data`` 辞書を作るための互換用ラッパー。

    新しい修正版スクリプトでは、この関数を経由せず、基本的に
    ``get_data()`` で複素時系列を取り出し、``get_Rangle()`` で位相と振幅を
    作る。これは、古い notebook や一時的な確認コードが
    ``phase_data["angle"]`` などを期待している場合のために残している。

    ここで決めている約束は次の通り。
    ENSO 位相 ``phi_e`` は ``PHI_e_SHIFT`` で Nino3.4 ピーク平均位相が
    0 になるように回転する。annual 位相 ``phi_a`` は 1 月が位相 0 に
    なるように回転する。返り値には
    元の複素時系列、位相、振幅、Nino3.4、位相基準のメタデータをまとめる。

    返り値の中身:
    - ``z``: ``load_z()`` で読んだ MATLAB ``Z.mat`` の内容。``time``、
      ``idx``、``angle`` を含む。元データの確認や SST との日付合わせに使う。
    - ``nino34``: ``load_nino34_mat()`` で読んだ Nino3.4 時系列。ENSO
      ピーク位相の確認や図中マーカーに使う。
    - ``series``: 回転補正後の複素時系列。``annual`` は 1 月基準に回転済み、
      ``enso`` は Nino3.4 ピーク基準に回転済み、``date`` は対応する日付。
      複素平面や実部・虚部の時系列図ではこの系列を使う。
    - ``raw_angle``: 回転補正前の ``raw_annual`` / ``raw_enso`` から計算した
      unwrap 済み位相。補正前後の位相を比較したいときの確認用。
    - ``angle``: 回転補正後の ``series`` から計算した unwrap 済み位相。
      ``phi_a`` と ``phi_e`` として phase composite や Q 推定に使う。
    - ``amp``: 回転補正後の ``series`` の振幅。回転は絶対値を変えないので、
      実質的には元の振幅と同じ。amp filter の判定に使う。
    - ``phase_origin``: どの位相基準を使ったかを記録するメタデータ。
      ``phi_a_shift_rad``、``phi_e_shift_rad``、ENSO ピーク情報、mode index
      などを含み、推定結果ファイルにも保存する。
    """
    data = get_data(enso_dir=enso_dir)
    z = data["z"]
    nino34 = data["nino34"]
    phi_a, amp_a = get_Rangle(data["annual"])
    phi_e, amp_e = get_Rangle(data["enso"])

    series = {"enso": data["enso"], "annual": data["annual"], "date": data["date"]}
    raw_angle = {
        "enso": phase_from_complex(data["raw_enso"]),
        "annual": phase_from_complex(data["raw_annual"]),
        "date": data["date"],
    }

    angle = {"enso": phi_e, "annual": phi_a, "date": data["date"]}
    amp = {
        "enso": amp_e,
        "annual": amp_a,
        "date": data["date"],
    }
    return {
        "z": z,
        "nino34": nino34,
        "series": series,
        "raw_angle": raw_angle,
        "angle": angle,
        "amp": amp,
        "phase_origin": data["phase_origin"],
    }


prepare_phase_data_no_enso_shift = prepare_phase_data


def amplitude_filter_mask(amp: np.ndarray, percentile_rank: float) -> np.ndarray:
    """Select samples whose amplitude is at or above the requested percentile."""
    threshold = np.percentile(np.asarray(amp, dtype=float), float(percentile_rank))
    return np.asarray(amp, dtype=float) >= threshold


def amplitude_filter_segment_mask(amp: np.ndarray, percentile_rank: float) -> np.ndarray:
    """
    位相 composite 用に、振幅 filter 後の有効点を選ぶ。

    振幅が閾値以上でも、前後と連続していない孤立点は位相速度を作れない。
    そのため、連続区間の長さが 2 点以上ある部分だけを ``True`` にする。
    """
    indices = np.where(amplitude_filter_mask(amp, percentile_rank))[0]
    split_points = np.where(np.diff(indices) != 1)[0] + 1
    valid = np.zeros(np.asarray(amp).size, dtype=bool)
    for segment in np.split(indices, split_points):
        if segment.size > 1:
            valid[segment] = True
    return valid
