import pytest

import scphylo as scp

from ._helpers import skip_rpy2


class TestDatasets:
    @skip_rpy2("oncoNEM")
    def test_simulate_1(self):
        df_in = scp.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0.001, beta=0.4, missing=0.2
        )
        is_cf = scp.ul.is_conflict_free_gusfield(df_in)
        assert not is_cf

    @skip_rpy2("oncoNEM")
    def test_simulate_2(self):
        df_ground = scp.datasets.simulate(
            n_cells=100, n_muts=100, n_clones=5, alpha=0, beta=0, missing=0
        )
        df_noisy = scp.datasets.add_noise(df_ground, alpha=0.001, beta=0.4, missing=0.2)
        assert not scp.ul.is_conflict_free_gusfield(df_noisy)

    @pytest.mark.parametrize(
        "func, n, m", [(scp.datasets.example, 83, 452), (scp.datasets.test, 20, 20)]
    )
    def test_load_datasets_grid(self, func, n, m):
        adata = func()
        assert adata.shape == (n, m)

    def test_load_datasets(self):
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
        assert adata.shape == (78, 25)
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
