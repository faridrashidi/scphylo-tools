import pytest
from click.testing import CliRunner

from scphylo.commands.scphylo import cli

from ._helpers import skip_graphviz


class TestCommandsSolver:
    def setup_method(self):
        self.runner = CliRunner()

    def test_scistree(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "scistree",
                test_data,
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_huntress_both(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "huntress",
                test_data,
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_scite(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "scite",
                test_data,
                "0.0000001",
                "0.1",
                "-r 3",
                "-l 1000",
            ],
        )
        assert result.exit_code == 0

    def test_scite_experiment(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "scite",
                test_data,
                "0.0000001",
                "0.1",
                "-e",
                "-t 1",
                "-s 1",
            ],
        )
        assert result.exit_code == 0

    def test_phiscsb(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "phiscsb",
                test_data,
                "0.0000001",
                "0.1",
            ],
        )
        assert result.exit_code == 0

    def test_consensus(
        self, test_dir, test_consensus_biorxiv_fig3b, test_consensus_biorxiv_figs18a
    ):
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "consensus",
                test_consensus_biorxiv_fig3b,
                test_consensus_biorxiv_figs18a,
                f"{test_dir}/consensus.CFMatrix",
            ],
        )
        assert result.exit_code == 0

    def test_cf2newick(self, test_dir):
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "cf2newick",
                f"{test_dir}/test.phiscsb.CFMatrix",
            ],
        )
        assert result.exit_code == 0

    def test_search(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "search",
                test_data,
                "-p 2",
            ],
        )
        assert result.exit_code == 0

    def test_score(self, test_cf_data_1, test_cf_data_2):
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "score",
                test_cf_data_1,
                test_cf_data_2,
            ],
        )
        assert result.exit_code == 0

    @skip_graphviz
    def test_cf2tree(self, test_dir):
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "cf2tree",
                f"{test_dir}/test.phiscsb.CFMatrix",
            ],
        )
        assert result.exit_code == 0

    @pytest.mark.skip(reason="Error Don't know!")
    def test_partf(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "partf",
                test_data,
                "0.0001",
                "0.1",
                "--n_threads 2",
                "--n_samples 100",
            ],
        )
        assert result.exit_code == 0

    @pytest.mark.skip(reason="PyTest issue with multithreading!")
    def test_booster(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "booster",
                test_data,
                "0.0000001",
                "0.1",
                "--solver scite",
                "--n_samples 100",
                "--sample_size 15",
                "--n_jobs 2",
                "--n_iterations 10000",
            ],
        )
        assert result.exit_code == 0

    @pytest.mark.skip(reason="pyBnB issue with CLI!")
    def test_bnb(self, test_data):
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "bnb",
                test_data,
                "-b simulated",
            ],
        )
        assert result.exit_code == 0
