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
