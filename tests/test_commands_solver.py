"""Exercise solver and utility commands through the public CLI."""

from click.testing import CliRunner

from scphylo.commands import cli


class TestCommandsSolver:
    """Verify command-line solver and utility workflows."""

    def setup_method(self):
        """Create an isolated Click command runner for each test."""
        self.runner = CliRunner()

    def test_scistree(self, test_data):
        """Verify that the ScisTree command completes successfully."""
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
        """Verify that the HUNTRESS command handles both error types."""
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
        """Verify that the SCITE command completes successfully."""
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
        """Verify that SCITE supports experiment-mode output."""
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
        """Verify that the binary PhISCS command completes successfully."""
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
        """Verify that the consensus utility combines two matrix files."""
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
        """Verify conversion from a conflict-free matrix to Newick."""
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
        """Verify that the search utility accepts a process count."""
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
        """Verify that the score utility compares two matrix files."""
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

    def test_cf2tree(self, test_dir):
        """Verify conversion from a conflict-free matrix to a tree."""
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "cf2tree",
                f"{test_dir}/test.phiscsb.CFMatrix",
            ],
        )
        assert result.exit_code == 0

    def test_partf(self, test_data):
        """Verify that the partition-function utility completes successfully."""
        result = self.runner.invoke(
            cli,
            [
                "utils",
                "partf",
                test_data,
                "0.0001",
                "0.1",
                "--n_threads",
                "2",
                "--n_samples",
                "100",
            ],
        )
        assert result.exit_code == 0

    def test_booster(self, test_data):
        """Verify that the booster command runs with the SCITE backend."""
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "booster",
                test_data,
                "0.0000001",
                "0.1",
                "--solver",
                "scite",
                "--n_samples",
                "100",
                "--sample_size",
                "15",
                "--n_jobs",
                "2",
                "--n_iterations",
                "10000",
            ],
        )
        assert result.exit_code == 0

    def test_bnb(self, test_data):
        """Verify that the branch-and-bound command accepts simulated bounds."""
        result = self.runner.invoke(
            cli,
            [
                "solver",
                "bnb",
                test_data,
                "-b",
                "simulated",
            ],
        )
        assert result.exit_code == 0
