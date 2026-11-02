# SAR Sea Ice Mapping (geoskill-sar-sea-ice-mapping)

> Maps sea ice type and concentration from SAR σ⁰: three-class Otsu in dB scale separates open water / young ice / multi-year ice, GLCM texture refines multi-year ice, sliding windows compute ice concentration, outputting a sea ice type GeoTIFF, a concentration GeoTIFF and area statistics JSON. SAR sea ice type and concentration mapping via multi-Otsu backscatter thresholds and GLCM texture.

---

## 1. Overview

Maps sea ice type and concentration from a single-temporal SAR σ⁰ (linear power). Physical basis: - **Open water**: specular reflection, very low σ⁰ (~−24 dB), uniform texture. - **Young ice**: newly formed ice surface, moderate σ⁰, fairly uniform texture. - **Multi-year ice**: repeated freeze-thaw, rough surface, high σ⁰ and strong texture. Method flow: 1. **Three-class Otsu in dB scale**: apply 3-class Otsu on `10·log10(σ⁰)` to separate open water / young ice / multi-year ice (log scale + three thresholds prevent young ice from being merged into water; automatically falls back to 2 classes for bimodal scenes). 2. **GLCM texture refinement**: multi-year ice = high σ⁰ and high GLCM contrast, remaining ice is classified as young ice. `--season summer` raises the texture threshold (summer melt weakens texture, so multi-year ice is judged more conservatively). 3. **Concentration**: fraction of ice pixels within the `--window` sliding window (ice concentration ∈ [0,1]).

## 3. Quick Start

```bash
pip install -r requirements.txt
python geoskill-sar-sea-ice-mapping.py --bbox 116 39 117 40 --synthetic --output-dir ./out
```

## 4. CLI Parameters

Run `python <skill>.py --help` for the full list. Common parameters:

| Parameter | Type | Description |
|---|---|---|
| `--bbox` | `float[4]` | WGS84 bounding box `min_lon min_lat max_lon max_lat` |
| `--input` | `path` | Local input file (GeoJSON/GeoTIFF/etc.) |
| `--output-dir` | `path` | Output directory (default `./output`) |
| `--synthetic` | `flag` | Use synthetic data instead of real input |
| `--quiet` | `flag` | Suppress non-essential stdout |

## 5. Input / Output

| File | Format | Description |
|---|---|---|
| `ice_type.tif` | GeoTIFF (uint8) | Sea ice type (0=open water 1=young ice 2=multi-year ice), EPSG:4326 |
| `ice_concentration.tif` | GeoTIFF (float32) | Ice concentration [0,1] |
| `ice_statistics.json` | JSON | Per-class pixels / share / area, mean concentration, thresholds |
| `output-manifest.json` | JSON | Run manifest |

## 6. Technical Principle

(see SKILL.md for details)

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

> 中文完整版：[references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)
