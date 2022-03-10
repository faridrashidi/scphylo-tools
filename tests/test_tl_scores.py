import numpy as np
import pytest

import scphylo as scp


class TestScores:
    def test_gs(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        gs = scp.tl.gs(grnd, sol)
        assert np.abs(gs - 0.9895) < 0.0001

    def test_ad(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        ad = scp.tl.ad(grnd, sol)
        assert np.abs(ad - 0.9778) < 0.0001

    def test_dl(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        dl = scp.tl.dl(grnd, sol)
        assert np.abs(dl - 0.9880) < 0.0001

    def test_cc(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        scp.tl.cc(grnd, sol)
        assert True

    def test_tpted(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        tpted = scp.tl.tpted(grnd, sol)
        assert np.abs(tpted - 0.8811) < 0.0001

    def test_caset(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        caset = scp.tl.caset(grnd, sol)
        assert np.abs(caset - 0.7847) < 0.0001

    def test_disc(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        disc = scp.tl.disc(grnd, sol)
        assert np.abs(disc - 0.7762) < 0.0001

    def test_mp3(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        mp3 = scp.tl.mp3(grnd, sol)
        assert np.abs(mp3 - 0.6582) < 0.002

    def test_rf(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        rf = scp.tl.rf(grnd, sol)
        assert np.abs(rf - 0.4864) < 0.0001

    @pytest.mark.skip(
        reason="Using MLTD in two tests is taking so long in test_scores!"
    )
    def test_mltd(self, test_cf_data_1, test_cf_data_2):
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        mltd = scp.tl.mltd(grnd, sol)
        assert np.abs(mltd["normalized_similarity"] - 0.7800) < 0.0001
