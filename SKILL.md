---
name: geoskill-sar-sea-ice-mapping
description: '从SAR σ⁰制图海冰类型与密集度：dB刻度三类Otsu分水面/新冰/多年冰，GLCM纹理精炼多年冰，滑窗计算冰密集度，输出海冰类型GeoTIFF、密集度GeoTIFF与面积统计JSON。SAR sea ice type and concentration mapping via multi-Otsu backscatter thresholds and GLCM texture.'
---

# SAR 海冰制图 | SAR Sea Ice Mapping

Maps sea ice type and concentration from single-temporal SAR σ⁰ (linear power). Physical basis:

- **Open water**: specular reflection, very low σ⁰ (~−24 dB) and uniform texture.
- **Young ice**: newly formed ice surface, moderate σ⁰ and fairly uniform texture.
- **Multi-year ice**: repeatedly frozen and thawed, rough surface, high σ⁰ and strong texture.

Processing workflow:

1. **Three-class Otsu on a dB scale**: 3-class Otsu on `10·log10(σ⁰)` separates open water / young ice / multi-year ice (log scale + three thresholds prevent young ice from being merged into open water; bimodal scenes automatically fall back to 2 classes).
2. **GLCM texture refinement**: multi-year ice = high σ⁰ and high GLCM contrast; the remaining ice is classified as young ice. `--season summer` raises the texture threshold (summer melt weakens texture, so multi-year ice is judged more conservatively).
3. **Concentration**: fraction of ice pixels within the `--window` sliding window (ice concentration ∈ [0,1]).

## Dependencies / 依赖

```bash
pip install numpy rasterio scipy
```

## Usage / 使用方法

### Basic usage (bbox only, synthetic data auto-generated)

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120.0 75.0 122.0 77.0 --output-dir ./out
```

### Example 1: synthetic data (offline, winter)

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120 75 122 77 --synthetic --season winter --output-dir ./syn
```

### Example 2: real SAR σ⁰ imagery

```bash
python geoskill-sar-sea-ice-mapping.py --input sigma0_linear.tif --output-dir ./real
```

### Example 3: summer scenario

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120 75 122 77 --season summer --output-dir ./summer --quiet
```

### Example 4: custom concentration window

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120 75 122 77 --window 15 --output-dir ./w15 --quiet
```

## Output / 输出

| File | Format | Description |
|---|---|---|
| `ice_type.tif` | GeoTIFF (uint8) | Sea ice type (0=open water 1=young ice 2=multi-year ice), EPSG:4326 |
| `ice_concentration.tif` | GeoTIFF (float32) | Ice concentration [0,1] |
| `ice_statistics.json` | JSON | Per-class pixels / share / area, mean concentration, thresholds |
| `output-manifest.json` | JSON | Run manifest |

## Data Source / 数据源 / Source

- **Real mode**: local SAR σ⁰ GeoTIFF (linear power).
- **Synthetic mode**: locally generated sea-surface background + young / multi-year ice areas (with different σ⁰ and texture).

## Privacy / 隐私声明 / Privacy

- Fully offline by default; `--synthetic` mode makes no network calls.
- All processing is done locally; no user data is uploaded.

## License / License

MIT

---

> 中文完整版：需要时再读 [references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)，不要默认全量加载。
