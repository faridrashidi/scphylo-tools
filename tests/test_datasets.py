"""Verify simulation, noise generation, and bundled dataset loaders."""

import hashlib

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

    def test_colorectal3_complete_matrix_and_subset(self):
        """Verify the full SiFit input and historical 48-cell view."""
        adata = scp.datasets.colorectal3()

        assert adata.shape == (61, 77)
        assert adata.obs_names.tolist() == (
            [f"c{i}" for i in range(1, 36)]
            + [f"ap{i}" for i in range(1, 14)]
            + [f"n{i}" for i in range(1, 14)]
        )
        assert adata.var_names.tolist() == [f"mutation_{i}" for i in range(1, 78)]
        np.testing.assert_array_equal(adata.var["sifit_index"], np.arange(1, 78))
        assert adata.obs_names.is_unique
        assert adata.var_names.is_unique

        values, counts = np.unique(adata.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 2867, 1: 1388, 3: 442}
        matrix_hash = hashlib.sha256(
            np.asarray(adata.X, dtype=np.int8).tobytes(order="C")
        ).hexdigest()
        assert (
            matrix_hash
            == "b57bff2721e598facc04e0a628846a343dc94a2cb486458f6c08ee91859701cd"
        )

        assert adata.obs["source_tissue"].value_counts(sort=False).to_dict() == {
            "cancer_biopsy": 35,
            "adenomatous_polyp": 13,
            "matched_normal_colorectal": 13,
        }
        non_normal = adata.obs["non_normal_tissue"].to_numpy()
        assert non_normal.sum() == 48
        assert adata.obs_names[non_normal].tolist() == (
            [f"c{i}" for i in range(1, 36)] + [f"ap{i}" for i in range(1, 14)]
        )

        historical = np.asarray(adata[non_normal].X, dtype=np.int8)
        values, counts = np.unique(historical, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 1996, 1: 1387, 3: 313}
        assert (
            hashlib.sha256(historical.tobytes(order="C")).hexdigest()
            == "bda1fb67cea6789754c3478ca79dd274442df23b8cbe61c82a631b4c3b1d119a"
        )
        derivative = adata.uns["non_normal_tissue_48"]
        assert (
            derivative["obs_names"].astype(str).tolist()
            == adata.obs_names[non_normal].tolist()
        )
        assert "SiFit Figure 5 itself used all 61 cells" in derivative["relationship"]

        four_gamete_violations = 0
        for first in range(adata.n_vars):
            for second in range(first + 1, adata.n_vars):
                pair = np.asarray(adata.X[:, [first, second]], dtype=np.int8)
                pair = pair[(pair != 3).all(axis=1)]
                observed = {tuple(row) for row in pair}
                four_gamete_violations += len(observed) == 4
        assert four_gamete_violations == 1847

        provenance = adata.uns["provenance"]
        assert (
            provenance["wu_fig3a"]["figure_sha256"]
            == "661cf7b640139c7c0b506c42011434730dc5e8ca7e8038e91ace651c88af8690"
        )
        assert (
            provenance["sifit_fig5"]["supplement_pdf_sha256"]
            == "019015b66bf82757fc1c1bf0d23e5b1bb4c625c49ccfbc7dd75205d53ab2a865"
        )
        assert (
            provenance["sifit_fig5"]["embedded_image_sha256"]
            == "b0214a5f8e02d1d5fd3f6ed7b3f32c9bb1742a4cad1190e65344fc95f2da0219"
        )
        assert provenance["sifit_fig5"]["missing_entries"] == 442
        assert provenance["sifit_fig5"]["four_gamete_violations"] == 1847
        assert (
            provenance["sifit_archive"]["runner_sha256"]
            == "33dcaa43b18eb5e4edf1a88cc3135df277fdfdc01cc14b76213541709e6188d4"
        )

    def test_acute_lymphocytic_leukemia2_variants(self):
        """Verify the full ALL2 input and historical PhISCS derivative."""
        adata = scp.datasets.acute_lymphocytic_leukemia2()

        assert adata.shape == (115, 16)
        assert adata.obs_names.tolist() == [f"cell{i}" for i in range(115)]
        values, counts = np.unique(adata.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 1139, 1: 701}
        np.testing.assert_array_equal(
            adata.X.sum(axis=0),
            [18, 48, 27, 20, 40, 50, 27, 11, 41, 91, 89, 48, 99, 18, 38, 36],
        )

        excluded = [
            "cell2",
            "cell8",
            "cell18",
            "cell24",
            "cell30",
            "cell36",
            "cell40",
            "cell43",
            "cell46",
            "cell100",
            "cell102",
            "cell105",
            "cell107",
        ]
        assert adata.obs["phiscs_fig9"].sum() == 102
        assert adata.var["phiscs_fig9"].all()
        assert adata.obs_names[~adata.obs["phiscs_fig9"]].tolist() == excluded
        np.testing.assert_array_equal(
            adata[excluded].X.sum(axis=1),
            [9, 7, 10, 16, 6, 6, 5, 9, 8, 10, 7, 15, 7],
        )

        phiscs_fig9 = adata.uns["phiscs_fig9"]
        obs_names = phiscs_fig9["obs_names"].astype(str).tolist()
        var_names = phiscs_fig9["var_names"].astype(str).tolist()
        assert obs_names == adata.obs_names[adata.obs["phiscs_fig9"]].tolist()
        assert var_names == adata.var_names[adata.var["phiscs_fig9"]].tolist()

        historical = adata[obs_names, var_names]
        values, counts = np.unique(historical.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 1046, 1: 586}

        solution = phiscs_fig9["solution_fig9"]
        assert solution.shape == (102, 16)
        assert np.isnan(solution).sum() == 306
        eliminated = np.isnan(solution).all(axis=0)
        assert adata.var_names[eliminated].tolist() == ["CMTM8", "RIMS2", "RRP8"]
        solution_df = historical.to_df().loc[:, ~eliminated]
        solution_df.iloc[:, :] = solution[:, ~eliminated]
        assert scp.ul.is_conflict_free_gusfield(solution_df)
        assert scp.ul.count_flips(historical.X, solution)[:2] == (62, 4)

        assert "solution_fig9" not in adata.layers
        assert "params_fig9" not in adata.uns
        assert phiscs_fig9["params_fig9"] == {
            "alpha": 0.001,
            "beta": 0.181749,
            "delta": 0.2,
            "kmax": 3,
            "w": 0,
        }
        assert (
            adata.uns["provenance"]["genotype"]["infscite_sha256"]
            == "088d754af50dcc1f8f6cc2b08aa3c878eca311461328e61c81959ac56da94e21"
        )

    def test_acute_lymphocytic_leukemia3_provenance(self):
        """Verify the published ALL3 matrix and distinct raw-cell cohort."""
        adata = scp.datasets.acute_lymphocytic_leukemia3()

        assert adata.shape == (150, 49)
        assert adata.obs_names.tolist() == [f"cell{i}" for i in range(150)]
        values, counts = np.unique(adata.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 4399, 1: 2951}
        matrix_hash = hashlib.sha256(
            np.asarray(adata.X, dtype="<i8").tobytes(order="C")
        ).hexdigest()
        assert (
            matrix_hash
            == "5ebfcfbe915970bb754e09d37345f657672afcc861174c4713600040cda7dd79"
        )
        assert adata.obs["gawad_cluster"].value_counts(sort=False).to_dict() == {
            1: 31,
            2: 46,
            3: 9,
            4: 34,
            5: 30,
        }

        provenance = adata.uns["provenance"]
        assert (
            provenance["genotype"]["source_sha256"]
            == "ad25e8216801c79b0043b4c8eefd2e7a371d4c63e6feefcb350a168b6dd12723"
        )
        reanalysis = provenance["later_255_cell_reanalysis"]
        assert reanalysis["accession"] == "SRP044380"
        assert (
            reanalysis["sciphi_sample_map_sha256"]
            == "d41672fea422d9adbe74d75dca8c701532036e184ded1eeafb15ffb606b0bc91"
        )
        assert reanalysis["scistree_called_snv_sites"] == 406
        assert "not a 255 x 49 extension" in reanalysis["relationship"]
        assert adata.uns["params_scite"] == {
            "alpha": 0.000001,
            "beta": 0.2521591,
            "gama": 0,
        }

    def test_high_grade_serous_ovarian_cancer1_derivatives(self):
        """Verify the full HGSOC1 assay and the infSCITE Figure S22 view."""
        adata = scp.datasets.high_grade_serous_ovarian_cancer1()

        assert adata.shape == (588, 43)
        values, counts = np.unique(adata.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 10774,
            1: 4747,
            3: 9763,
        }
        assert (
            hashlib.sha256(
                np.asarray(adata.X, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "a290d79c2bd1f3885a22643d8f3811bdb04baf65e20d514d4dd93e9ccaf5af31"
        )

        assert adata.obs_names.is_unique
        assert adata.var_names.is_unique
        assert adata.obs["sample_id"].value_counts(sort=False).to_dict() == {
            "omentum_site_1": 126,
            "omentum_site_2": 126,
            "right_ovary_site_1": 84,
            "right_ovary_site_2": 252,
        }
        assert adata.obs["variation"].value_counts(sort=False).to_dict() == {
            "barcode": 5,
            "standard": 583,
        }
        obs_payload = (
            "\n".join(
                f"{str(index)}\t{str(row.sample_id)}\t{int(row.well_id)}\t"
                f"{int(row.library)}\t{str(row.variation)}\t"
                f"{int(row.total_depth)}\t"
                f"{int(bool(row.mcpherson_clonal_analysis))}"
                for index, row in zip(
                    adata.obs_names, adata.obs.itertuples(), strict=True
                )
            )
            + "\n"
        )
        assert (
            hashlib.sha256(obs_payload.encode()).hexdigest()
            == "a7866741c951240ce776ee8163a28289d371b8745cf7b6138d75b51f7d16c982"
        )

        assert {
            "allele_state",
            "alt_p_value",
            "mutant",
            "ref_p_value",
            "total",
        } <= set(adata.layers)
        allele_state = adata.layers["allele_state"]
        values, counts = np.unique(allele_state, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 10774,
            1: 3642,
            2: 1105,
            3: 9763,
        }
        assert (
            hashlib.sha256(
                np.asarray(allele_state, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "9f3c7fff6a97e87fc27aa2a8de016a307f5a0593b2dd0792b20ce94125787c41"
        )
        np.testing.assert_array_equal(
            adata.X,
            np.where(allele_state == 3, 3, (allele_state > 0).astype(np.int8)),
        )

        mutant = adata.layers["mutant"]
        total = adata.layers["total"]
        called = total >= 50
        ref_present = adata.layers["ref_p_value"] < 1e-6
        alt_present = adata.layers["alt_p_value"] < 1e-6
        expected_state = np.full(adata.shape, 3, dtype=np.int8)
        expected_state[called & ref_present & ~alt_present] = 0
        expected_state[called & ref_present & alt_present] = 1
        expected_state[called & ~ref_present & alt_present] = 2
        np.testing.assert_array_equal(allele_state, expected_state)
        assert (
            hashlib.sha256(
                np.asarray(adata.layers["ref_p_value"], dtype="<f8").tobytes(order="C")
            ).hexdigest()
            == "43aeefc5dec6b480de0045b4e34c754c07678adedc5753525213d19e01f065e2"
        )
        assert (
            hashlib.sha256(
                np.asarray(adata.layers["alt_p_value"], dtype="<f8").tobytes(order="C")
            ).hexdigest()
            == "c3b8eeb3966d6a89444d5ba1d1aa75dcd5a95c7c255943fe20708d24ec7e702f"
        )
        assert mutant.sum() == 4748593
        assert total.sum() == 32181831
        assert np.all(mutant >= 0)
        assert np.all(mutant <= total)
        assert (
            hashlib.sha256(
                np.asarray(mutant, dtype="<i8").tobytes(order="C")
            ).hexdigest()
            == "05d77be325d05d2c4abddc98c480568f036009e6f7bb03f3f4879a01e6850d34"
        )
        assert (
            hashlib.sha256(
                np.asarray(total, dtype="<i8").tobytes(order="C")
            ).hexdigest()
            == "bfd72f4a6b661bd4cad6522467d0171a9fb2081e005fe303471407ca849bbfa0"
        )

        mcpherson = adata.obs["mcpherson_clonal_analysis"].to_numpy()
        np.testing.assert_array_equal(mcpherson, called.any(axis=1))
        assert mcpherson.sum() == 491
        assert np.all(adata.X[~mcpherson] == 3)

        assert adata.var["category"].value_counts(sort=False).to_dict() == {
            "normal_marker": 6,
            "snv": 37,
        }
        infscite = adata.var["infscite_fig_s22"].to_numpy()
        np.testing.assert_array_equal(infscite, adata.var["category"] == "snv")
        assert infscite.sum() == 37
        infscite_matrix = adata[:, infscite].X
        values, counts = np.unique(infscite_matrix, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 9552,
            1: 3993,
            3: 8211,
        }
        assert (
            hashlib.sha256(
                np.asarray(infscite_matrix, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "663cc02b40002530c8293ff34e768765f74a6a709cf352a5e8e6142e0153a0a6"
        )

        normal_markers = {
            (str(row.chrom), int(row.position), str(row.ref), str(row.alt))
            for row in adata.var.loc[~infscite].itertuples()
        }
        assert normal_markers == {
            ("10", 95791613, "T", "C"),
            ("10", 99220707, "C", "T"),
            ("10", 104491702, "G", "A"),
            ("10", 120445794, "G", "A"),
            ("10", 131321853, "T", "C"),
            ("17", 1732300, "T", "C"),
        }
        var_payload = (
            "\n".join(
                f"{str(index)}\t{str(row.chrom)}\t{int(row.position)}\t"
                f"{str(row.ref)}\t{str(row.alt)}\t{str(row.gene_name)}\t"
                f"{str(row.effect)}\t{str(row.category)}\t"
                f"{int(bool(row.infscite_fig_s22))}"
                for index, row in zip(
                    adata.var_names, adata.var.itertuples(), strict=True
                )
            )
            + "\n"
        )
        assert (
            hashlib.sha256(var_payload.encode()).hexdigest()
            == "2a8bf9037d09a73da1cc7ecbe3f101bc71599c8057e46b77843f2335f8f29754"
        )

        provenance = adata.uns["provenance"]
        assert (
            provenance["source"]["workbook_sha256"]
            == "166bf71ae2c8f686fa11e97784de1619ef1b6339becc65a647c2a66fc78bdbd6"
        )
        assert (
            provenance["source"]["methods_pdf_sha256"]
            == "7915b4e3772878a1ddb568abd574e6b868ac7fd130cab497de7be3e9deec6844"
        )
        figure_s22 = provenance["infscite_fig_s22"]
        assert (
            figure_s22["supplement_pdf_sha256"]
            == "c1a07a492a40d117884d43fd6c90077709ce16498babebbfe2aed0e1907a7de5"
        )
        assert figure_s22["false_positive_rate"] == 0.015
        assert figure_s22["false_signal_rate"] == 0.465
        assert figure_s22["recurrent_gene"] == "AC004538.3"
        full_panel = provenance["infscite_full_panel"]
        assert full_panel["n_loci"] == 43
        assert full_panel["false_positive_rate"] == 0.020
        assert full_panel["false_signal_rate"] == 0.436

    def test_high_grade_serous_ovarian_cancer2_derivatives(self):
        """Verify the full HGSOC2 assay and the infSCITE Figure S23 view."""
        adata = scp.datasets.high_grade_serous_ovarian_cancer2()

        assert adata.shape == (672, 84)
        values, counts = np.unique(adata.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 18283,
            1: 10299,
            3: 27866,
        }
        assert (
            hashlib.sha256(
                np.asarray(adata.X, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "419839ab9e60320ddcb35b0a2ea5ea6a717151d650204fbceb874e3a03d17fe1"
        )

        assert adata.obs_names.is_unique
        assert adata.var_names.is_unique
        assert adata.obs["sample_id"].value_counts(sort=False).to_dict() == {
            "adnexa_site_1": 168,
            "omentum_site_1": 168,
            "right_ovary_site_1": 168,
            "right_ovary_site_2": 168,
        }
        obs_payload = (
            "\n".join(
                f"{str(index)}\t{str(row.sample_id)}\t{int(row.well_id)}\t"
                f"{int(row.library)}\t{int(row.total_depth)}\t"
                f"{int(bool(row.mcpherson_clonal_analysis))}"
                for index, row in zip(
                    adata.obs_names, adata.obs.itertuples(), strict=True
                )
            )
            + "\n"
        )
        assert (
            hashlib.sha256(obs_payload.encode()).hexdigest()
            == "e0987bd38e75d099ef819e2e8c2aac61614b447f6276adbab7744530f81d064b"
        )

        assert {
            "allele_state",
            "alt_p_value",
            "mutant",
            "ref_p_value",
            "total",
        } <= set(adata.layers)
        allele_state = adata.layers["allele_state"]
        values, counts = np.unique(allele_state, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 18283,
            1: 8533,
            2: 1766,
            3: 27866,
        }
        assert (
            hashlib.sha256(
                np.asarray(allele_state, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "5f186d0526743635244259b4a49d5e4e5450ea816a1fb0fe3fa6984d6215d011"
        )
        np.testing.assert_array_equal(
            adata.X,
            np.where(allele_state == 3, 3, (allele_state > 0).astype(np.int8)),
        )

        mutant = adata.layers["mutant"]
        total = adata.layers["total"]
        called = total >= 50
        ref_present = adata.layers["ref_p_value"] < 1e-6
        alt_present = adata.layers["alt_p_value"] < 1e-6
        expected_state = np.full(adata.shape, 3, dtype=np.int8)
        expected_state[called & ref_present & ~alt_present] = 0
        expected_state[called & ref_present & alt_present] = 1
        expected_state[called & ~ref_present & alt_present] = 2
        np.testing.assert_array_equal(allele_state, expected_state)
        assert (
            hashlib.sha256(
                np.asarray(adata.layers["ref_p_value"], dtype="<f8").tobytes(order="C")
            ).hexdigest()
            == "366620ae2634a83bc893d05c418da4f763cfa7c91493698b6cb4db144a153356"
        )
        assert (
            hashlib.sha256(
                np.asarray(adata.layers["alt_p_value"], dtype="<f8").tobytes(order="C")
            ).hexdigest()
            == "ff559a21f225794615ac4406abb1807961d92ba30b9a9d9aca1d2f96170f50b8"
        )
        assert mutant.sum() == 3963366
        assert total.sum() == 20062693
        assert np.all(mutant >= 0)
        assert np.all(mutant <= total)
        assert (
            hashlib.sha256(
                np.asarray(mutant, dtype="<i8").tobytes(order="C")
            ).hexdigest()
            == "1f01cc77f2fbcb844b6e723fdc81b1f1eeb279adf6e36ad9f0672a0643595701"
        )
        assert (
            hashlib.sha256(
                np.asarray(total, dtype="<i8").tobytes(order="C")
            ).hexdigest()
            == "a18f57295cf9485dc7e55436c76e53f35bd874b4813a32ba45c17299ec621549"
        )

        mcpherson = adata.obs["mcpherson_clonal_analysis"].to_numpy()
        np.testing.assert_array_equal(mcpherson, called.any(axis=1))
        assert mcpherson.sum() == 636
        assert np.all(adata.X[~mcpherson] == 3)

        assert adata.var["source_category"].value_counts(sort=False).to_dict() == {
            "snv": 84
        }
        assert adata.var["category"].value_counts(sort=False).to_dict() == {
            "normal_marker": 24,
            "snv": 60,
        }
        infscite = adata.var["infscite_fig_s23"].to_numpy()
        np.testing.assert_array_equal(infscite, adata.var["category"] == "snv")
        assert infscite.sum() == 60
        infscite_matrix = adata[:, infscite].X
        values, counts = np.unique(infscite_matrix, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 13547,
            1: 6782,
            3: 19991,
        }
        assert (
            hashlib.sha256(
                np.asarray(infscite_matrix, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "0467b9e97691b6f6d826996689de56e810f35a15cde4d086797c02c8b79dc16f"
        )

        normal_markers = {
            (str(row.chrom), int(row.position), str(row.ref), str(row.alt))
            for row in adata.var.loc[~infscite].itertuples()
        }
        assert normal_markers == {
            ("5", 55662214, "C", "A"),
            ("5", 56681043, "G", "A"),
            ("9", 121883165, "A", "G"),
            ("9", 128901785, "C", "A"),
            ("9", 134175658, "A", "G"),
            ("11", 17818550, "G", "A"),
            ("11", 41693309, "C", "T"),
            ("12", 14348017, "A", "G"),
            ("12", 19105851, "G", "A"),
            ("12", 32738983, "A", "G"),
            ("12", 43149770, "A", "T"),
            ("12", 86807499, "T", "A"),
            ("12", 120614093, "C", "T"),
            ("14", 78572568, "A", "G"),
            ("15", 46405908, "T", "G"),
            ("15", 46820930, "A", "G"),
            ("15", 48112038, "G", "A"),
            ("15", 52306097, "T", "C"),
            ("15", 92460694, "A", "G"),
            ("17", 10015197, "G", "A"),
            ("17", 33989694, "C", "G"),
            ("21", 15729517, "T", "C"),
            ("21", 15849874, "G", "A"),
            ("21", 16227138, "A", "T"),
        }
        var_payload = (
            "\n".join(
                f"{str(index)}\t{str(row.chrom)}\t{int(row.position)}\t"
                f"{str(row.ref)}\t{str(row.alt)}\t{str(row.gene_name)}\t"
                f"{str(row.effect)}\t{str(row.source_category)}\t"
                f"{str(row.category)}\t{int(bool(row.infscite_fig_s23))}"
                for index, row in zip(
                    adata.var_names, adata.var.itertuples(), strict=True
                )
            )
            + "\n"
        )
        assert (
            hashlib.sha256(var_payload.encode()).hexdigest()
            == "0757a674c7bc6247494c0b8faa3084d8772a94e3582769d3e9866d519ffe2864"
        )

        provenance = adata.uns["provenance"]
        source = provenance["source"]
        assert (
            source["workbook_sha256"]
            == "166bf71ae2c8f686fa11e97784de1619ef1b6339becc65a647c2a66fc78bdbd6"
        )
        assert (
            source["methods_pdf_sha256"]
            == "7915b4e3772878a1ddb568abd574e6b868ac7fd130cab497de7be3e9deec6844"
        )
        assert "source_category preserves" in source["category_rule"]
        figure_s23 = provenance["infscite_fig_s23"]
        assert (
            figure_s23["supplement_pdf_sha256"]
            == "c1a07a492a40d117884d43fd6c90077709ce16498babebbfe2aed0e1907a7de5"
        )
        assert figure_s23["n_loci"] == 60
        assert figure_s23["false_positive_rate"] == 0.0098
        assert figure_s23["false_signal_rate"] == 0.438
        assert figure_s23["recurrent_chrom"] == "5"
        assert figure_s23["recurrent_position"] == 50248591
        assert figure_s23["recurrent_substitution"] == "T>C"
        assert (
            figure_s23["included_locus_ids_sha256"]
            == "0b36d2ef7af8c0f8314dba91db2c8b3c6eecda108c01b17e20a3abf3934abd56"
        )
        full_panel = provenance["infscite_full_panel"]
        assert full_panel["n_loci"] == 84
        assert full_panel["false_positive_rate"] == 0.016
        assert full_panel["false_signal_rate"] == 0.442
        assert (
            full_panel["excluded_locus_ids_sha256"]
            == "14654890081ff41f04d78fd02caed5ccc50a62473e5a836c6d125ffdd8f95cd0"
        )

    def test_high_grade_serous_ovarian_cancer3_derivatives(self):
        """Verify the full HGSOC3 assay and its published subset dimensions."""
        adata = scp.datasets.high_grade_serous_ovarian_cancer3()

        assert adata.shape == (420, 43)
        values, counts = np.unique(adata.X, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 8099, 1: 4390, 3: 5571}
        matrix_hash = hashlib.sha256(
            np.asarray(adata.X, dtype=np.int8).tobytes(order="C")
        ).hexdigest()
        assert (
            matrix_hash
            == "7e59ef1d8c004f09b6cff293639a2af2dbc523fc400ea413963538ec1341d792"
        )

        assert adata.obs_names.is_unique
        assert adata.var_names.is_unique
        assert adata.obs["sample_id"].value_counts(sort=False).to_dict() == {
            "left_ovary_site_1": 84,
            "left_ovary_site_2": 84,
            "omentum_site_1": 84,
            "omentum_site_2": 84,
            "right_ovary_site_1": 84,
        }
        obs_payload = (
            "\n".join(
                f"{str(index)}\t{str(row.sample_id)}\t{int(row.well_id)}\t"
                f"{int(row.library)}\t{int(row.total_depth)}\t"
                f"{int(bool(row.mcpherson_clonal_analysis))}\t"
                f"{int(bool(row.sciphi_scvilp))}"
                for index, row in zip(
                    adata.obs_names, adata.obs.itertuples(), strict=True
                )
            )
            + "\n"
        )
        assert (
            hashlib.sha256(obs_payload.encode()).hexdigest()
            == "87e4298d0970005bb711fe898ed698cabd844b1d3d1bac1360ae7c1f18366a79"
        )

        assert {
            "allele_state",
            "alt_p_value",
            "mutant",
            "ref_p_value",
            "total",
        } <= set(adata.layers)
        allele_state = adata.layers["allele_state"]
        values, counts = np.unique(allele_state, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {
            0: 8099,
            1: 3804,
            2: 586,
            3: 5571,
        }
        state_hash = hashlib.sha256(
            np.asarray(allele_state, dtype=np.int8).tobytes(order="C")
        ).hexdigest()
        assert (
            state_hash
            == "fe368c36858bee47973b5f5ea1f360789d93a91af5d73cc995acca4953b5aa7b"
        )
        np.testing.assert_array_equal(
            adata.X,
            np.where(allele_state == 3, 3, (allele_state > 0).astype(np.int8)),
        )

        mutant = adata.layers["mutant"]
        total = adata.layers["total"]
        called = total >= 50
        ref_present = adata.layers["ref_p_value"] < 1e-6
        alt_present = adata.layers["alt_p_value"] < 1e-6
        expected_state = np.full(adata.shape, 3, dtype=np.int8)
        expected_state[called & ref_present & ~alt_present] = 0
        expected_state[called & ref_present & alt_present] = 1
        expected_state[called & ~ref_present & alt_present] = 2
        np.testing.assert_array_equal(allele_state, expected_state)
        assert (
            hashlib.sha256(
                np.asarray(adata.layers["ref_p_value"], dtype="<f8").tobytes(order="C")
            ).hexdigest()
            == "8ac765021982d29527c77212ee23d3b91cf237bfec5e3377f1fff8af2c42bcd9"
        )
        assert (
            hashlib.sha256(
                np.asarray(adata.layers["alt_p_value"], dtype="<f8").tobytes(order="C")
            ).hexdigest()
            == "19c8bf9c999334b3819129fcbe6e0ed242ae6956682bf7eabf578a469ddbd708"
        )
        assert mutant.sum() == 11979881
        assert total.sum() == 50114355
        assert np.all(mutant >= 0)
        assert np.all(mutant <= total)
        assert (
            hashlib.sha256(
                np.asarray(mutant, dtype="<i8").tobytes(order="C")
            ).hexdigest()
            == "666c8d56112dc5e1ebae49176d4f9483bdb782e16068833240dc3610eaac15fc"
        )
        assert (
            hashlib.sha256(
                np.asarray(total, dtype="<i8").tobytes(order="C")
            ).hexdigest()
            == "cb847f65ac6919d2a88253f2d3e73ea0d9de812c8dbb220899ac75a989101d47"
        )

        mcpherson = adata.obs["mcpherson_clonal_analysis"].to_numpy()
        np.testing.assert_array_equal(mcpherson, (total >= 50).any(axis=1))
        assert mcpherson.sum() == 402
        assert np.all(adata.X[~mcpherson] == 3)

        sciphi_scvilp = adata.obs["sciphi_scvilp"].to_numpy()
        cell_depth = total.sum(axis=1)
        np.testing.assert_array_equal(sciphi_scvilp, cell_depth > 10000)
        assert sciphi_scvilp.sum() == 370
        assert cell_depth[~sciphi_scvilp].max() == 9993
        assert cell_depth[sciphi_scvilp].min() == 10649
        assert adata.obs.loc[sciphi_scvilp, "sample_id"].value_counts(
            sort=False
        ).to_dict() == {
            "left_ovary_site_1": 76,
            "left_ovary_site_2": 73,
            "omentum_site_1": 67,
            "omentum_site_2": 83,
            "right_ovary_site_1": 71,
        }

        assert adata.var["category"].value_counts(sort=False).to_dict() == {
            "normal_marker": 6,
            "snv": 37,
        }
        infscite = adata.var["infscite_fig_s24"].to_numpy()
        np.testing.assert_array_equal(infscite, adata.var["category"] == "snv")
        assert infscite.sum() == 37
        infscite_matrix = adata[:, infscite].X
        values, counts = np.unique(infscite_matrix, return_counts=True)
        assert dict(zip(values, counts, strict=True)) == {0: 6354, 1: 4212, 3: 4974}
        assert (
            hashlib.sha256(
                np.asarray(infscite_matrix, dtype=np.int8).tobytes(order="C")
            ).hexdigest()
            == "ca719fa00fd47844c5ae50e9d655de4d3f708f16d7afaac96d0893d0ae885ea9"
        )

        normal_markers = {
            (str(row.chrom), int(row.position), str(row.ref), str(row.alt))
            for row in adata.var.loc[~infscite].itertuples()
        }
        assert normal_markers == {
            ("17", 41245466, "G", "A"),
            ("17", 57754591, "A", "G"),
            ("17", 78333840, "T", "C"),
            ("17", 79511135, "A", "G"),
            ("22", 50869692, "G", "C"),
            ("X", 100395663, "G", "T"),
        }
        var_payload = (
            "\n".join(
                f"{str(index)}\t{str(row.chrom)}\t{int(row.position)}\t"
                f"{str(row.ref)}\t{str(row.alt)}\t{str(row.gene_name)}\t"
                f"{str(row.effect)}\t{str(row.category)}\t"
                f"{int(bool(row.infscite_fig_s24))}"
                for index, row in zip(
                    adata.var_names, adata.var.itertuples(), strict=True
                )
            )
            + "\n"
        )
        assert (
            hashlib.sha256(var_payload.encode()).hexdigest()
            == "30ec96fe2c1a41f944835723db290054b73ec184548789d970d7683eacaf54ed"
        )

        provenance = adata.uns["provenance"]
        assert (
            provenance["source"]["workbook_sha256"]
            == "166bf71ae2c8f686fa11e97784de1619ef1b6339becc65a647c2a66fc78bdbd6"
        )
        assert (
            provenance["source"]["methods_pdf_sha256"]
            == "7915b4e3772878a1ddb568abd574e6b868ac7fd130cab497de7be3e9deec6844"
        )
        assert (
            provenance["infscite_fig_s24"]["supplement_pdf_sha256"]
            == "c1a07a492a40d117884d43fd6c90077709ce16498babebbfe2aed0e1907a7de5"
        )
        assert (
            provenance["sciphi_scvilp"]["scvilp_counts_sha256"]
            == "540cdef971fd0dc1508338960f9d5a5a08aeeb6747e57138be28a2cb100468b8"
        )
        assert (
            provenance["sciphi_scvilp"]["sciphi_supplement_sha256"]
            == "3a750b7fcaedd6587084a54c653221d1a6ecb1c682dcd953d98b439b7d832847"
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
        assert adata.shape == (115, 16)
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
        adata = scp.datasets.colorectal3()
        assert adata.shape == (61, 77)
        adata = scp.datasets.erbc()
        assert adata.shape == (47, 40)
        adata = scp.datasets.high_grade_serous_ovarian_cancer_3celllines()
        assert adata.shape == (891, 14068)
        adata = scp.datasets.high_grade_serous_ovarian_cancer1()
        assert adata.shape == (588, 43)
        adata = scp.datasets.high_grade_serous_ovarian_cancer2()
        assert adata.shape == (672, 84)
        adata = scp.datasets.high_grade_serous_ovarian_cancer3()
        assert adata.shape == (420, 43)
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
        # adata = scp.datasets.acute_lymphocytic_leukemia_many()
