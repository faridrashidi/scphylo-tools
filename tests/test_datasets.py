"""Verify simulation, noise generation, and bundled dataset loaders."""

import numpy as np
import pytest

import scphylo as scp

from ._helpers import skip_rpy2


class TestDatasets:
    """Exercise dataset generation and loading helpers."""

    @skip_rpy2("oncoNEM")
    def test_simulate_1(self):
        """Verify that noisy OncoNEM simulation creates conflicts."""
        df_in = scp.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0.001, beta=0.4, missing=0.2
        )
        is_cf = scp.ul.is_conflict_free_gusfield(df_in)
        assert not is_cf

    @skip_rpy2("oncoNEM")
    def test_simulate_2(self):
        """Verify that adding noise to a clean simulation creates conflicts."""
        df_ground = scp.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0, beta=0, missing=0
        )
        df_noisy = scp.datasets.add_noise(df_ground, alpha=0.001, beta=0.4, missing=0.2)
        assert not scp.ul.is_conflict_free_gusfield(df_noisy)

    def test_add_noise(self):
        """Verify that added observation noise creates genotype conflicts."""
        df = scp.datasets.test()
        df[df == 3] = 0
        df_noisy = scp.datasets.add_noise(df, alpha=0.001, beta=0.1, missing=0.2)
        assert not scp.ul.is_conflict_free_gusfield(df_noisy)

    def test_add_readcount(self):
        """Verify that simulated read-count layers preserve matrix shape."""
        df = scp.datasets.test()
        adata = scp.datasets.add_readcount(df, mean_coverage=60, seed=5)
        assert adata.layers["total"].shape == df.shape
        assert adata.layers["mutant"].shape == df.shape

    @pytest.mark.parametrize(
        "func, n, m", [(scp.datasets.example, 83, 452), (scp.datasets.test, 20, 20)]
    )
    def test_load_datasets_grid(self, func, n, m):
        """Verify the expected shape of each parameterized dataset."""
        adata = func()
        assert adata.shape == (n, m)

    def test_colorectal2_variants(self):
        """Verify the full genotype and read-count variants stay aligned."""
        genotype = scp.datasets.colorectal2()
        readcount = scp.datasets.colorectal2(readcount=True)

        assert genotype.shape == readcount.shape == (182, 36)
        assert genotype.obs_names.equals(readcount.obs_names)
        assert genotype.var_names.equals(readcount.var_names)
        np.testing.assert_array_equal(genotype.X, readcount.X)
        values, counts = np.unique(genotype.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 5174, 1: 866, 3: 512}
        assert (genotype.X == 1).any(axis=1).sum() == 86

        assert genotype.obs["phiscs_fig7"].sum() == 78
        assert genotype.var["phiscs_fig7"].sum() == 25
        assert genotype.obs["copy_number_profile"].value_counts().to_dict() == {
            0: 117,
            1: 32,
            2: 10,
            3: 23,
        }

        phiscs_fig7 = genotype.uns["phiscs_fig7"]
        assert len(phiscs_fig7["obs_names"]) == 78
        assert len(phiscs_fig7["var_names"]) == 25
        assert phiscs_fig7["solution_fig7a"].shape == (78, 25)
        assert phiscs_fig7["solution_fig7b"].shape == (78, 25)

        assert {"mutant", "total"} <= set(readcount.layers)
        mutant = readcount.layers["mutant"]
        total = readcount.layers["total"]
        assert np.isfinite(mutant).all() and np.isfinite(total).all()
        assert np.all(mutant >= 0)
        assert np.all(mutant <= total)
        assert np.equal(mutant, np.floor(mutant)).all()
        assert np.equal(total, np.floor(total)).all()
        assert mutant.sum() == 208955
        assert total.sum() == 1994705
        assert np.all(total.sum(axis=1) > 0)
        assert readcount["MA_94", "TP53_chr17_7577548"].layers["mutant"].item() == 248
        assert readcount["MA_94", "TP53_chr17_7577548"].layers["total"].item() == 250

        bulk_columns = [column for column in genotype.var if "Count_CO8" in column]
        assert len(bulk_columns) == 4
        assert (
            genotype.var.loc[genotype.var["phiscs_fig7"], bulk_columns]
            .notna()
            .all()
            .all()
        )
        assert (
            genotype.var.loc[~genotype.var["phiscs_fig7"], bulk_columns]
            .isna()
            .all()
            .all()
        )

    def test_load_datasets(self):
        """Verify the expected shapes of all bundled datasets."""
        adata = scp.datasets.example()
        assert adata.shape == (83, 452)
        adata = scp.datasets.test()
        assert adata.shape == (20, 20)

        adata = scp.datasets.acute_lymphocytic_leukemia1()
        assert adata.shape == (111, 20)
        adata = scp.datasets.acute_lymphocytic_leukemia2()
        assert adata.shape == (102, 16)
        adata = scp.datasets.acute_lymphocytic_leukemia3()
        assert adata.shape == (150, 49)
        adata = scp.datasets.acute_lymphocytic_leukemia4()
        assert adata.shape == (143, 78)
        adata = scp.datasets.acute_lymphocytic_leukemia5()
        assert adata.shape == (96, 105)
        adata = scp.datasets.acute_lymphocytic_leukemia6()
        assert adata.shape == (146, 10)
        adata = scp.datasets.colorectal1()
        assert adata.shape == (178, 16)
        adata = scp.datasets.colorectal2()
        assert adata.shape == (182, 36)
        # adata = scp.datasets.colorectal3()
        adata = scp.datasets.erbc()
        assert adata.shape == (47, 40)
        adata = scp.datasets.high_grade_serous_ovarian_cancer_3celllines()
        assert adata.shape == (891, 14068)
        adata = scp.datasets.melanoma20()
        assert adata.shape == (20, 2367)
        adata = scp.datasets.muscle_invasive_bladder()
        assert adata.shape == (44, 443)
        adata = scp.datasets.myeloproliferative_neoplasms18()
        assert adata.shape == (58, 18)
        adata = scp.datasets.myeloproliferative_neoplasms78()
        assert adata.shape == (58, 78)
        adata = scp.datasets.myeloproliferative_neoplasms712()
        assert adata.shape == (58, 712)
        adata = scp.datasets.oligodendroglioma_idh_mutated_tumor()
        assert adata.shape == (579, 77)
        adata = scp.datasets.renal_cell_carcinoma()
        assert adata.shape == (17, 35)
        adata = scp.datasets.tnbc()
        assert adata.shape == (16, 20)
        # adata = scp.datasets.high_grade_serous_ovarian_cancer1()
        # adata = scp.datasets.high_grade_serous_ovarian_cancer2()
        # adata = scp.datasets.high_grade_serous_ovarian_cancer3()
        # adata = scp.datasets.acute_lymphocytic_leukemia_many()
