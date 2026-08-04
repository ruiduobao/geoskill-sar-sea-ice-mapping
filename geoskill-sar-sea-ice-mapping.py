#!/usr/bin/env python3
"""sar-sea-ice-mapping — SAR 海冰制图

从单时相 SAR σ⁰（线性功率）制图海冰类型与密集度。物理依据：

- **开放水面**：镜面反射，σ⁰ 极低、纹理均匀（低 GLCM 对比度）。
- **新冰（young ice）**：初生冰面，σ⁰ 中等、纹理较均匀。
- **多年冰（multi-year ice）**：经反复冻融，表面粗糙、含气泡少，σ⁰ 高且
  纹理强（高 GLCM 对比度）。

方法流程：

1. **GLCM 对比度纹理**（水平 + 垂直邻接差平方滑窗均值）。
2. **σ⁰ + 纹理联合分类**：用 Otsu 在 σ⁰ 上分出水面 / 冰；在冰内再按纹理
   强弱分新冰 / 多年冰。分位数阈值可按 ``--season`` 微调。
3. **密集度**：滑窗内冰像元占比（ice concentration ∈ [0,1]）。

输出海冰类型 GeoTIFF（0=水 1=新冰 2=多年冰）+ 密集度 GeoTIFF + 面积统计 JSON。

数据源：本地 SAR σ⁰ GeoTIFF（线性功率），或 ``--synthetic`` 生成海面背景 +
不同 σ⁰ / 纹理冰区的模拟场景用于离线验证。

隐私声明 / Privacy：
- 默认完全离线，``--synthetic`` 无网络。
- 所有处理本地完成，不上传用户数据。

Usage:
    python sar-sea-ice-mapping.py --input sigma0.tif --output-dir ./out
    python sar-sea-ice-mapping.py --bbox 120 75 122 77 --season winter --output-dir ./out

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

VERSION = "1.0.0"
SKILL_NAME = "sar-sea-ice-mapping"

try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover - fallback minimal definitions
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


# 类别编码
CLASS_WATER = 0
CLASS_YOUNG_ICE = 1
CLASS_MULTIYEAR_ICE = 2
CLASS_NAMES = {0: "open_water", 1: "young_ice", 2: "multiyear_ice"}


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 核心算法
# ---------------------------------------------------------------------------
def otsu_threshold(values: np.ndarray, n_bins: int = 256) -> float:
    """Otsu 最大类间方差阈值（取最大平台中点）。"""
    v = values[np.isfinite(values)]
    if v.size == 0:
        return 0.0
    vmin, vmax = float(v.min()), float(v.max())
    if vmax <= vmin:
        return float(vmin)
    hist, edges = np.histogram(v, bins=n_bins, range=(vmin, vmax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = hist.sum()
    if total == 0:
        return float(vmin)
    hist = hist.astype(np.float64)
    w_bg = np.cumsum(hist)
    w_fg = total - w_bg
    mean_cum = np.cumsum(hist * centers)
    mean_total = mean_cum[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        m_bg = mean_cum / w_bg
        m_fg = (mean_total - mean_cum) / w_fg
    between = w_bg * w_fg * (m_bg - m_fg) ** 2
    between = np.nan_to_num(between, nan=0.0, posinf=0.0, neginf=0.0)
    maxval = between.max()
    idxs = np.flatnonzero(between == maxval)
    idx = int(idxs[len(idxs) // 2]) if idxs.size else int(np.argmax(between))
    return float(centers[idx])


def multi_otsu(values: np.ndarray, n_bins: int = 128) -> Tuple[float, float]:
    """3 类 Otsu：返回两个阈值 ``(t_low, t_high)``，把样本分成低 / 中 / 高三类。

    最大化三类间方差。对双峰（退化）场景返回 ``(thr, thr)``（两阈值相等），
    调用方据此回退到 2 类逻辑。
    """
    v = values[np.isfinite(values)]
    if v.size == 0:
        return 0.0, 0.0
    vmin, vmax = float(v.min()), float(v.max())
    if vmax <= vmin:
        return vmin, vmin
    hist, edges = np.histogram(v, bins=n_bins, range=(vmin, vmax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    P = hist.astype(np.float64)
    total = P.sum()
    if total == 0:
        return vmin, vmin
    P = P / total
    omega = np.cumsum(P)
    mu = np.cumsum(P * centers)
    mu_T = float(mu[-1])
    best = -1.0
    t1 = t2 = -1
    for i in range(0, n_bins - 2):
        w0 = omega[i]
        if w0 <= 0:
            continue
        mu0 = mu[i] / w0
        for j in range(i + 1, n_bins - 1):
            w1 = omega[j] - omega[i]
            if w1 <= 0:
                continue
            w2 = 1.0 - omega[j]
            if w2 <= 0:
                continue
            mu1 = (mu[j] - mu[i]) / w1
            mu2 = (mu_T - mu[j]) / w2
            s = w0 * (mu0 - mu_T) ** 2 + w1 * (mu1 - mu_T) ** 2 + w2 * (mu2 - mu_T) ** 2
            if s > best:
                best = s
                t1, t2 = i, j
    if t1 < 0:
        thr = otsu_threshold(v, n_bins)
        return thr, thr
    return float(centers[t1]), float(centers[t2])


def glcm_contrast(gray: np.ndarray, levels: int = 32, window: int = 7) -> np.ndarray:
    """局部 GLCM 对比度纹理（水平 + 垂直邻接差平方的滑窗均值）。"""
    g = np.nan_to_num(np.asarray(gray, dtype=np.float32), nan=0.0)
    gmin, gmax = float(g.min()), float(g.max())
    if gmax <= gmin:
        return np.zeros_like(g, dtype=np.float32)
    levels = max(int(levels), 2)
    q = ((g - gmin) / (gmax - gmin) * (levels - 1)).round().astype(np.float32)
    dx2 = (q[:, 1:] - q[:, :-1]) ** 2
    dy2 = (q[1:, :] - q[:-1, :]) ** 2
    dx2 = np.pad(dx2, ((0, 0), (0, 1)), mode="edge")
    dy2 = np.pad(dy2, ((0, 1), (0, 0)), mode="edge")
    diff_sq = 0.5 * (dx2 + dy2)
    from scipy.ndimage import uniform_filter
    return uniform_filter(diff_sq, size=max(int(window), 1), mode="reflect").astype(np.float32)


def classify_ice(
    sigma0: np.ndarray,
    season: str = "winter",
    tex_levels: int = 32,
    tex_window: int = 7,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """σ⁰ + 纹理联合海冰类型分类。

    1. 在 **dB 刻度**上 3 类 Otsu 分 水面 / 新冰 / 多年冰（水面镜面反射 σ⁰
       极低；对数刻度 + 三阈值避免新冰被并入水面）。双峰场景自动回退 2 类。
    2. 纹理精炼：多年冰 = 高 σ⁰ 且高 GLCM 对比度；其余冰判为新冰。

    返回 ``(class_map_uint8, params)``，类别码见 CLASS_*。
    """
    s = np.nan_to_num(np.asarray(sigma0, dtype=np.float32), nan=0.0)
    s = np.clip(s, 0.0, None)
    if s.size == 0:
        raise ValidationError("input is empty")

    db = (10.0 * np.log10(np.clip(s, 1e-6, None))).astype(np.float32)
    t_low, t_high = multi_otsu(db)
    bimodal = t_high <= t_low
    if bimodal:
        t_low = otsu_threshold(db)
        t_high = t_low

    ice_mask = db > t_low
    texture = glcm_contrast(s, levels=tex_levels, window=tex_window)

    class_map = np.full(s.shape, CLASS_WATER, dtype=np.uint8)
    tex_thr = None
    if ice_mask.any():
        tex_thr = otsu_threshold(texture[ice_mask])
        multiyear = ice_mask & (db > t_high) & (texture > tex_thr)
        young = ice_mask & ~multiyear
        class_map[young] = CLASS_YOUNG_ICE
        class_map[multiyear] = CLASS_MULTIYEAR_ICE

    # 季节先验：夏季融冰减弱多年冰纹理 → 提高纹理门限（更保守地判多年冰）
    season_factor = 1.0 if season == "winter" else 1.25
    if tex_thr is not None and season_factor != 1.0:
        adj_thr = tex_thr * season_factor
        class_map2 = np.full(s.shape, CLASS_WATER, dtype=np.uint8)
        multiyear2 = ice_mask & (db > t_high) & (texture > adj_thr)
        young2 = ice_mask & ~multiyear2
        class_map2[young2] = CLASS_YOUNG_ICE
        class_map2[multiyear2] = CLASS_MULTIYEAR_ICE
        class_map = class_map2

    params = {
        "season": season,
        "season_factor": season_factor,
        "water_sigma0_threshold": float(10.0 ** (t_low / 10.0)),
        "water_sigma0_threshold_db": float(t_low),
        "multiyear_sigma0_threshold_db": float(t_high),
        "texture_threshold": None if tex_thr is None else float(tex_thr * season_factor),
        "texture_used": True,
        "bimodal_fallback": bool(bimodal),
    }
    return class_map, params


def ice_concentration(
    class_map: np.ndarray, window: int = 9
) -> np.ndarray:
    """冰密集度 = 滑窗内冰像元（非水）占比，∈[0,1]。"""
    ice = (class_map != CLASS_WATER).astype(np.float32)
    from scipy.ndimage import uniform_filter
    conc = uniform_filter(ice, size=max(int(window), 1), mode="reflect")
    return np.clip(conc, 0.0, 1.0).astype(np.float32)


def pixel_area_km2(bbox: List[float], height: int, width: int) -> float:
    lat_mid = 0.5 * (bbox[1] + bbox[3])
    km_per_deg_lon = 111.0 * float(np.cos(np.deg2rad(lat_mid)))
    px_w = (bbox[2] - bbox[0]) / max(width, 1) * km_per_deg_lon
    px_h = (bbox[3] - bbox[1]) / max(height, 1) * 111.0
    return float(abs(px_w * px_h))


def ice_statistics(
    class_map: np.ndarray, concentration: np.ndarray,
    bbox: List[float], params: Dict[str, Any],
) -> Dict[str, Any]:
    """逐类面积统计 + 平均密集度。"""
    h, w = class_map.shape
    px = pixel_area_km2(bbox, h, w)
    total = class_map.size
    per_class: Dict[str, Any] = {}
    for code, name in CLASS_NAMES.items():
        n = int((class_map == code).sum())
        per_class[name] = {
            "pixels": n,
            "fraction": float(n / total),
            "area_km2": float(n * px),
        }
    ice_frac = float((class_map != CLASS_WATER).mean())
    stats = {
        "total_pixels": total,
        "total_area_km2": float(total * px),
        "pixel_area_km2": px,
        "ice_fraction": ice_frac,
        "ice_area_km2": float(ice_frac * total * px),
        "mean_concentration": float(np.mean(concentration)),
        "per_class": per_class,
    }
    stats.update(params)
    return stats


# ---------------------------------------------------------------------------
# 合成数据：海面 + 新冰 + 多年冰
# ---------------------------------------------------------------------------
def generate_synthetic(
    bbox: List[float],
    season: str = "winter",
    width: int = 64,
    height: int = 64,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """生成 σ⁰ 场景：开放水面（低 σ⁰ 低纹理）+ 新冰（中 σ⁰）+ 多年冰（高 σ⁰ 高纹理）。

    返回 ``(sigma0, truth_class_map, info)``。
    """
    rng = np.random.default_rng(seed)
    sigma0 = np.full((height, width), 0.004, dtype=np.float32)  # 水面 ~-24 dB
    # 乘性斑点
    sigma0 = sigma0 * np.exp(rng.normal(0, 0.12, (height, width))).astype(np.float32)

    truth = np.zeros((height, width), dtype=np.uint8)  # 0=水

    # 新冰区：中 σ⁰（~-15 dB），低纹理（均匀）
    yi = (slice(6, 26), slice(6, 26))
    sigma0[yi] = 0.03 + rng.normal(0, 0.003, (20, 20)).astype(np.float32)
    truth[yi] = CLASS_YOUNG_ICE

    # 多年冰区：高 σ⁰（~-6 dB）+ 高纹理（随机团块）
    mi = (slice(34, 58), slice(34, 58))
    rr, cc = 24, 24
    base = 0.25
    # 高纹理：叠加随机高低斑块
    texture = rng.uniform(0.0, 0.25, (rr, cc)).astype(np.float32)
    checker = ((np.mgrid[0:rr, 0:cc].sum(axis=0)) % 2).astype(np.float32)
    sigma0[mi] = base + 0.12 * checker + texture
    truth[mi] = CLASS_MULTIYEAR_ICE

    sigma0 = np.clip(sigma0, 1e-5, None).astype(np.float32)
    info = {
        "bbox": bbox, "width": width, "height": height, "seed": seed,
        "season": season,
        "truth_water_fraction": float((truth == CLASS_WATER).mean()),
        "truth_young_fraction": float((truth == CLASS_YOUNG_ICE).mean()),
        "truth_multiyear_fraction": float((truth == CLASS_MULTIYEAR_ICE).mean()),
    }
    return sigma0, truth, info


# ---------------------------------------------------------------------------
# GeoTIFF I/O
# ---------------------------------------------------------------------------
def write_geotiff(
    path: str, cube: np.ndarray, bbox: List[float],
    nodata: float = -9999.0, dtype: str = "float32",
) -> None:
    import rasterio
    from rasterio.transform import from_bounds
    if cube.ndim == 2:
        cube = cube[np.newaxis, ...]
    nb, h, w = cube.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": dtype, "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(cube[b].astype(dtype), b + 1)


def read_geotiff(path: str) -> Tuple[np.ndarray, List[float]]:
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
    return cube, bbox


def read_geotiff_full(path: str) -> Tuple[np.ndarray, List[float], Optional[float]]:
    """扩展版 read：同时返回 nodata 值（若无则为 None）。"""
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
        nodata = src.nodata
        if nodata is not None:
            nodata = float(nodata)
    return cube, bbox, nodata


# ---------------------------------------------------------------------------
# 输入校验
# ---------------------------------------------------------------------------
def validate_bbox(bbox: Optional[List[float]]) -> List[float]:
    """校验地理 bbox 合法性，失败抛 ValidationError（exit 6）。"""
    if bbox is None or len(bbox) != 4:
        raise ValidationError("bbox must be 4 floats [W S E N]")
    try:
        w, s, e, n = [float(x) for x in bbox]
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"bbox entries must be numeric: {exc}")
    if not (-90.0 <= s <= 90.0 and -90.0 <= n <= 90.0):
        raise ValidationError(f"latitude out of [-90,90]: S={s}, N={n}")
    if not (-180.0 <= w <= 180.0 and -180.0 <= e <= 180.0):
        raise ValidationError(f"longitude out of [-180,180]: W={w}, E={e}")
    if s >= n:
        raise ValidationError(
            f"S >= N (S={s}, N={n}); bbox inverted (S must be < N)"
        )
    if w >= e:
        raise ValidationError(
            f"W >= E (W={w}, E={e}); cross-180° bbox not supported. "
            f"Split into two non-antipodal bboxes."
        )
    if (e - w) < 0.001 or (n - s) < 0.001:
        raise ValidationError(
            f"bbox too small ({(e-w):.6f}°×{(n-s):.6f}°); min span is 0.001°"
        )
    return [w, s, e, n]


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(
    output_dir, args, outputs, qa, started_at, exit_code, bbox,
    input_nodata: Optional[float] = None,
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME, skill_version=VERSION, command=cmd,
        started_at=started_at, finished_at=_utc_now(), exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "season": getattr(args, "season", None),
            "synthetic": bool(getattr(args, "synthetic", False)),
            "input_nodata": input_nodata,
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir
    bbox = list(args.bbox) if args.bbox else None
    truth = None
    synth_info: Optional[Dict[str, Any]] = None
    input_nodata: Optional[float] = None
    n_valid_pixels: Optional[int] = None

    # 校验 CLI 参数（前置）
    if args.window < 1:
        raise ValidationError(
            f"--window must be >= 1 (got {args.window})"
        )

    if args.input and not args.synthetic:
        cube, file_bbox, src_nodata = read_geotiff_full(args.input)
        input_nodata = src_nodata
        if bbox is not None:
            bbox = validate_bbox(bbox)
        else:
            bbox = validate_bbox(file_bbox)
        # NoData 处理
        if src_nodata is not None:
            n_total = int(cube[0].size)
            n_nd = int(np.count_nonzero(cube[0] == src_nodata))
            n_valid_pixels = n_total - n_nd
            if n_valid_pixels == 0:
                raise ValidationError(
                    f"input raster has no valid pixels "
                    f"(all {n_nd}/{n_total} are NoData={src_nodata})",
                    path=args.input, nodata=src_nodata,
                )
            cube = np.where(cube == src_nodata, np.nan, cube).astype(np.float32)
        else:
            n_valid_pixels = int(cube[0].size)
        sigma0 = cube[0]
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <raster>")
        bbox = validate_bbox(bbox)
        sigma0, truth, synth_info = generate_synthetic(bbox, season=args.season)
        n_valid_pixels = int(sigma0.size)
        source_note = "synthetic"

    if sigma0.size == 0:
        raise ValidationError("input raster is empty")
    if not np.any(np.isfinite(sigma0)):
        raise ValidationError("input raster has no finite values")

    # 校验通过后再建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 2) 分类 + 密集度
    class_map, params = classify_ice(sigma0, season=args.season)
    concentration = ice_concentration(class_map, window=args.window)
    stats = ice_statistics(class_map, concentration, bbox, params)

    # 3) 写出
    cls_tif = os.path.join(output_dir, "ice_type.tif")
    conc_tif = os.path.join(output_dir, "ice_concentration.tif")
    write_geotiff(cls_tif, class_map.astype("uint8"), bbox, nodata=255, dtype="uint8")
    write_geotiff(conc_tif, concentration, bbox)

    stats_path = os.path.join(output_dir, "ice_statistics.json")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note,
        "season": args.season,
        "ice_fraction": stats["ice_fraction"],
        "mean_concentration": stats["mean_concentration"],
        "per_class_fraction": {k: v["fraction"] for k, v in stats["per_class"].items()},
        "n_valid_pixels": int(n_valid_pixels) if n_valid_pixels is not None else None,
        "input_nodata": input_nodata,
    }
    if synth_info is not None:
        qa["synthetic_truth"] = {
            "water": synth_info["truth_water_fraction"],
            "young": synth_info["truth_young_fraction"],
            "multiyear": synth_info["truth_multiyear_fraction"],
        }

    outputs = [
        {"path": cls_tif, "kind": "raster", "crs_epsg": 4326, "bbox_wgs84": bbox, "band_count": 1},
        {"path": conc_tif, "kind": "raster", "crs_epsg": 4326, "bbox_wgs84": bbox, "band_count": 1},
        {"path": stats_path, "kind": "json"},
    ]
    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox,
                              input_nodata=input_nodata)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}  season: {args.season}")
        print(f"[{SKILL_NAME}] ice fraction: {stats['ice_fraction']:.4f}  "
              f"mean concentration: {stats['mean_concentration']:.4f}")
        for name, d in stats["per_class"].items():
            print(f"[{SKILL_NAME}]   {name}: {d['fraction']:.4f} ({d['area_km2']:.2f} km²)")
        print(f"[{SKILL_NAME}] output: {cls_tif}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="SAR sea ice type and concentration mapping from backscatter + GLCM texture.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input SAR σ⁰ GeoTIFF (linear power)")
    p.add_argument("--season", default="winter", choices=["winter", "summer"],
                   help="season prior for multi-year ice (default: winter)")
    p.add_argument("--window", type=int, default=9,
                   help="concentration sliding-window size in pixels (default: 9)")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a synthetic sea-ice SAR scene (offline)")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
