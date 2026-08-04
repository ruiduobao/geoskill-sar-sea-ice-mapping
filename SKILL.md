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

<!-- ===== 中文原文 (Chinese Original) ===== -->

---
name: geoskill-sar-sea-ice-mapping
description: '从SAR σ⁰制图海冰类型与密集度：dB刻度三类Otsu分水面/新冰/多年冰，GLCM纹理精炼多年冰，滑窗计算冰密集度，输出海冰类型GeoTIFF、密集度GeoTIFF与面积统计JSON。SAR sea ice type and concentration mapping via multi-Otsu backscatter thresholds and GLCM texture.'
---

# SAR 海冰制图 | SAR Sea Ice Mapping

从单时相 SAR σ⁰（线性功率）制图海冰类型与密集度。物理依据：

- **开放水面**：镜面反射，σ⁰ 极低（~−24 dB）、纹理均匀。
- **新冰（young ice）**：初生冰面，σ⁰ 中等、纹理较均匀。
- **多年冰（multi-year ice）**：反复冻融、表面粗糙，σ⁰ 高且纹理强。

方法流程：

1. **dB 刻度三类 Otsu**：在 `10·log10(σ⁰)` 上用 3 类 Otsu 分水面 / 新冰 /
   多年冰（对数刻度 + 三阈值避免新冰被并入水面；双峰场景自动回退 2 类）。
2. **GLCM 纹理精炼**：多年冰 = 高 σ⁰ 且高 GLCM 对比度，其余冰判为新冰。
   `--season summer` 会提高纹理门限（夏季融冰减弱纹理，更保守判多年冰）。
3. **密集度**：`--window` 滑窗内冰像元占比（ice concentration ∈ [0,1]）。

## 依赖

```bash
pip install numpy rasterio scipy
```

## 使用方法

### 基本用法（仅给 bbox，自动合成）

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120.0 75.0 122.0 77.0 --output-dir ./out
```

### 示例 1：合成数据（离线，冬季）

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120 75 122 77 --synthetic --season winter --output-dir ./syn
```

### 示例 2：真实 SAR σ⁰ 影像

```bash
python geoskill-sar-sea-ice-mapping.py --input sigma0_linear.tif --output-dir ./real
```

### 示例 3：夏季场景

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120 75 122 77 --season summer --output-dir ./summer --quiet
```

### 示例 4：自定义密集度窗口

```bash
python geoskill-sar-sea-ice-mapping.py --bbox 120 75 122 77 --window 15 --output-dir ./w15 --quiet
```

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `ice_type.tif` | GeoTIFF (uint8) | 海冰类型（0=水面 1=新冰 2=多年冰），EPSG:4326 |
| `ice_concentration.tif` | GeoTIFF (float32) | 冰密集度 [0,1] |
| `ice_statistics.json` | JSON | 逐类像元 / 占比 / 面积、平均密集度、阈值 |
| `output-manifest.json` | JSON | 运行清单 |

## 数据源 / Source

- **真实模式**：本地 SAR σ⁰ GeoTIFF（线性功率）。
- **合成模式**：本地生成海面背景 + 新冰 / 多年冰区（不同 σ⁰ 与纹理）。

## 隐私声明 / Privacy

- 默认完全离线运行，`--synthetic` 无任何网络。
- 所有处理本地完成，不上传用户数据。

## License

MIT
