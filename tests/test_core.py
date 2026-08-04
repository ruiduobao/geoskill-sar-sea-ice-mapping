"""Core algorithm tests for sar-sea-ice-mapping."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod as si


class TestOtsu:
    def test_bimodal(self):
        rng = np.random.default_rng(0)
        data = np.concatenate([rng.normal(5, 1, 3000), rng.normal(50, 1, 3000)]).astype(np.float32)
        thr = si.otsu_threshold(data)
        assert 5 < thr < 50

    def test_constant(self):
        assert si.otsu_threshold(np.full(50, 3.0, dtype=np.float32)) == pytest.approx(3.0)

    def test_empty(self):
        assert si.otsu_threshold(np.full(5, np.nan, dtype=np.float32)) == 0.0


class TestMultiOtsu:
    def test_three_modes(self):
        rng = np.random.default_rng(0)
        data = np.concatenate([
            rng.normal(-24, 0.5, 3000),
            rng.normal(-15, 0.5, 1000),
            rng.normal(-5, 0.5, 1000),
        ]).astype(np.float32)
        t_low, t_high = si.multi_otsu(data)
        assert -24 < t_low < -15
        assert -15 < t_high < -5

    def test_two_modes_degenerate(self):
        # 仅两个离散取值（直方图只有 2 个占用 bin）→ 无法切出非空中类 → 退化
        data = np.concatenate([np.full(2000, 5.0), np.full(2000, 50.0)]).astype(np.float32)
        t_low, t_high = si.multi_otsu(data)
        assert t_low == t_high


class TestGlcm:
    def test_uniform_zero(self):
        c = si.glcm_contrast(np.full((32, 32), 100.0, dtype=np.float32))
        assert c.max() == pytest.approx(0.0)

    def test_textured_higher(self):
        yy, xx = np.mgrid[0:32, 0:32]
        check = (((yy + xx) % 2) * 200).astype(np.float32)
        flat = np.full((32, 32), 100.0, dtype=np.float32)
        assert si.glcm_contrast(check).mean() > si.glcm_contrast(flat).mean() + 1.0


class TestClassify:
    def test_separates_water_ice(self):
        sigma0 = np.full((64, 64), 0.004, dtype=np.float32)  # 水
        sigma0[20:40, 20:40] = 0.3  # 冰（高 σ⁰）
        cm, params = si.classify_ice(sigma0, season="winter")
        # 冰区不是水
        assert (cm[20:40, 20:40] != si.CLASS_WATER).mean() > 0.9
        # 水区大部分是水
        water_bg = np.ones((64, 64), dtype=bool)
        water_bg[20:40, 20:40] = False
        assert (cm[water_bg] == si.CLASS_WATER).mean() > 0.9

    def test_multiyear_vs_young_texture(self):
        """高纹理冰判为多年冰，低纹理冰判为新冰。"""
        sigma0 = np.full((64, 64), 0.004, dtype=np.float32)
        # 新冰：均匀中 σ⁰
        sigma0[5:25, 5:25] = 0.05
        # 多年冰：高 σ⁰ + 强纹理（棋盘）
        yy, xx = np.mgrid[0:20, 0:20]
        checker = ((yy + xx) % 2).astype(np.float32)
        sigma0[38:58, 38:58] = 0.25 + 0.2 * checker
        cm, params = si.classify_ice(sigma0, season="winter")
        assert (cm[38:58, 38:58] == si.CLASS_MULTIYEAR_ICE).mean() > 0.5
        assert (cm[5:25, 5:25] == si.CLASS_YOUNG_ICE).mean() > 0.5

    def test_class_codes_valid(self):
        sigma0, truth, info = si.generate_synthetic([120, 75, 122, 77])
        cm, _ = si.classify_ice(sigma0)
        assert set(np.unique(cm)).issubset({0, 1, 2})


class TestConcentration:
    def test_all_ice(self):
        cm = np.full((16, 16), si.CLASS_YOUNG_ICE, dtype=np.uint8)
        conc = si.ice_concentration(cm)
        assert conc.mean() == pytest.approx(1.0, abs=0.05)

    def test_all_water(self):
        cm = np.zeros((16, 16), dtype=np.uint8)
        conc = si.ice_concentration(cm)
        assert conc.max() == pytest.approx(0.0)

    def test_half(self):
        cm = np.zeros((40, 40), dtype=np.uint8)
        cm[:, :20] = si.CLASS_MULTIYEAR_ICE
        conc = si.ice_concentration(cm, window=5)
        # 中心区接近 0.5（过渡带有平滑）
        assert 0.35 < conc[20, 20] < 0.65


class TestSynthetic:
    def test_shapes(self):
        sigma0, truth, info = si.generate_synthetic([120, 75, 122, 77])
        assert sigma0.shape == (64, 64)
        assert truth.shape == (64, 64)
        assert info["truth_young_fraction"] > 0
        assert info["truth_multiyear_fraction"] > 0

    def test_detection_matches_truth(self):
        """分类结果应与注入真值有较好一致性。"""
        sigma0, truth, info = si.generate_synthetic([120, 75, 122, 77], seed=7)
        cm, _ = si.classify_ice(sigma0, season="winter")
        # 冰 / 水分离的准确度
        ice_truth = truth != si.CLASS_WATER
        ice_pred = cm != si.CLASS_WATER
        iou = np.logical_and(ice_truth, ice_pred).sum() / max(np.logical_or(ice_truth, ice_pred).sum(), 1)
        assert iou > 0.7, f"ice/water IoU too low: {iou:.3f}"

    def test_positive_sigma0(self):
        sigma0, _, _ = si.generate_synthetic([120, 75, 122, 77])
        assert sigma0.min() > 0


class TestStatistics:
    def test_per_class_fractions_sum_one(self):
        cm = np.zeros((10, 10), dtype=np.uint8)
        cm[0:5, :] = si.CLASS_YOUNG_ICE
        cm[5:10, 0:5] = si.CLASS_MULTIYEAR_ICE
        conc = si.ice_concentration(cm)
        stats = si.ice_statistics(cm, conc, [120, 75, 122, 77], {"season": "winter"})
        total_frac = sum(d["fraction"] for d in stats["per_class"].values())
        assert total_frac == pytest.approx(1.0)
        assert stats["ice_fraction"] == pytest.approx(0.75)
        assert stats["per_class"]["multiyear_ice"]["pixels"] == 25


class TestIO:
    def test_roundtrip(self, tmp_path):
        arr = np.random.randint(0, 3, (16, 16)).astype(np.uint8)
        path = str(tmp_path / "c.tif")
        si.write_geotiff(path, arr, [120.0, 75.0, 122.0, 77.0], nodata=255, dtype="uint8")
        back, bbox = si.read_geotiff(path)
        np.testing.assert_allclose(back[0], arr)

    def test_missing_raises(self):
        with pytest.raises(si.UsageError):
            si.read_geotiff("/nonexistent/x.tif")
