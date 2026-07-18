"""Verify tree and genotype reconstruction scores."""

import numpy as np
import pandas as pd

import scphylo as scp


class TestScores:
    """Compare score implementations with known fixture values."""

    def test_gs(self, test_cf_data_1, test_cf_data_2):
        """Verify the genotype-similarity score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        gs = scp.tl.gs(grnd, sol)
        assert np.abs(gs - 0.9895) < 0.0001

    def test_ad(self, test_cf_data_1, test_cf_data_2):
        """Verify the ancestor-descendant score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        ad = scp.tl.ad(grnd, sol)
        assert np.abs(ad - 0.9778) < 0.0001

    def test_dl(self, test_cf_data_1, test_cf_data_2):
        """Verify the different-lineages score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        dl = scp.tl.dl(grnd, sol)
        assert np.abs(dl - 0.9880) < 0.0001

    def test_cc(self, test_cf_data_1, test_cf_data_2):
        """Verify that the co-clustering score can be computed."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        scp.tl.cc(grnd, sol)
        assert True

    def test_tpted(self, test_cf_data_1, test_cf_data_2):
        """Verify the tree-path tree-edit-distance score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        tpted = scp.tl.tpted(grnd, sol)
        assert np.abs(tpted - 0.8811) < 0.0001

    def test_tpted_ignores_sibling_order(self):
        """Ignore sibling order introduced by differing clone frequencies."""
        columns = ["a", "b", "c"]
        ground = pd.DataFrame(
            [
                [1, 1, 0],
                [1, 0, 1],
                [1, 0, 1],
                [0, 0, 0],
            ],
            columns=columns,
        ).rename(index=str)
        reordered = pd.DataFrame(
            [
                [1, 1, 0],
                [1, 1, 0],
                [1, 0, 1],
                [0, 0, 0],
            ],
            columns=columns,
        ).rename(index=str)
        chain = pd.DataFrame(
            [
                [1, 0, 0],
                [1, 0, 0],
                [1, 0, 1],
                [0, 1, 0],
                [0, 1, 0],
                [0, 0, 0],
            ],
            columns=columns,
        ).rename(index=str)

        assert scp.tl.tpted(ground, reordered) == 1.0
        assert scp.tl.tpted(reordered, ground) == 1.0
        assert scp.tl.tpted(ground, chain) == 0.75
        assert scp.tl.tpted(chain, ground) == 0.75

    def test_caset(self, test_cf_data_1, test_cf_data_2):
        """Verify the CASet score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        caset = scp.tl.caset(grnd, sol)
        assert np.abs(caset - 0.7847) < 0.0001

    def test_disc(self, test_cf_data_1, test_cf_data_2):
        """Verify the DISC score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        disc = scp.tl.disc(grnd, sol)
        assert np.abs(disc - 0.8226) < 0.0001

    def test_disc_mutations_in_same_node(self):
        """Treat matching empty distinctly inherited sets as identical."""
        same_node = pd.DataFrame(
            [[1, 1], [0, 0]],
            index=["tumor", "normal"],
            columns=["a", "b"],
        )
        chain = pd.DataFrame(
            [[1, 1], [1, 0], [0, 0]],
            index=["both", "a_only", "normal"],
            columns=["a", "b"],
        )

        assert scp.tl.disc(same_node, same_node) == 1.0
        assert scp.tl.disc(same_node, chain) == 0.5
        assert scp.tl.disc(chain, same_node) == 0.5

    def test_mp3(self, test_cf_data_1, test_cf_data_2):
        """Verify the MP3 similarity score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        mp3 = scp.tl.mp3(grnd, sol)
        assert np.abs(mp3 - 0.6582) < 0.002

    def test_rf(self, test_cf_data_1, test_cf_data_2):
        """Verify the Robinson-Foulds similarity score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        rf = scp.tl.rf(grnd, sol)
        assert np.abs(rf - 0.4864) < 0.0001

    def test_mltd(self, test_cf_data_1, test_cf_data_2):
        """Verify the multi-labeled tree distance similarity score."""
        grnd = scp.io.read(test_cf_data_1)
        sol = scp.io.read(test_cf_data_2)
        mltd = scp.tl.mltd(grnd, sol)
        assert np.abs(mltd["normalized_similarity"] - 0.7800) < 0.0001
