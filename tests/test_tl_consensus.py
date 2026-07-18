"""Verify consensus-tree algorithms against published fixtures."""

import networkx as nx

import scphylo as scp


class TestConsensus:
    """Exercise the primary consensus-tree implementation."""

    def test_consensus_1(
        self, test_consensus_biorxiv_fig3b, test_consensus_biorxiv_figs18a
    ):
        """Verify the supplementary figure 18 consensus tree size."""
        # result in biorxiv.figs18b
        sc1 = scp.io.read(test_consensus_biorxiv_fig3b)
        sc2 = scp.io.read(test_consensus_biorxiv_figs18a)
        final_tree = scp.tl.consensus(sc1, sc2)
        assert len(final_tree.nodes) == 19

    def test_consensus_2(
        self, test_consensus_recomb_fig1a, test_consensus_recomb_fig1b
    ):
        """Verify the RECOMB figure 1 consensus tree size."""
        # result in recomb.fig1c
        sc1 = scp.io.read(test_consensus_recomb_fig1a)
        sc2 = scp.io.read(test_consensus_recomb_fig1b)
        final_tree = scp.tl.consensus(sc1, sc2)
        assert len(final_tree.nodes) == 21

    def test_consensus_3(
        self, test_consensus_biorxiv_fig4b, test_consensus_biorxiv_fig4c
    ):
        """Verify the BioRxiv figure 4 consensus tree size."""
        # result in biorxiv.fig4d
        sc1 = scp.io.read(test_consensus_biorxiv_fig4b)
        sc2 = scp.io.read(test_consensus_biorxiv_fig4c)
        final_tree = scp.tl.consensus(sc1, sc2)
        assert len(final_tree.nodes) == 34

    def test_consensus_4(
        self, test_consensus_biorxiv_fig3b, test_consensus_biorxiv_fig3c
    ):
        """Verify the BioRxiv figure 3 consensus tree size."""
        # result in biorxiv.fig3f
        sc1 = scp.io.read(test_consensus_biorxiv_fig3b)
        sc2 = scp.io.read(test_consensus_biorxiv_fig3c)
        final_tree = scp.tl.consensus(sc1, sc2)
        assert len(final_tree.nodes) == 19


class TestConsensusDay:
    """Compare the Day consensus implementation with the primary algorithm."""

    def test_consensus_day_1(
        self, test_consensus_biorxiv_fig3b, test_consensus_biorxiv_figs18a
    ):
        """Verify isomorphism for the supplementary figure 18 inputs."""
        # result in biorxiv.figs18b
        sc1 = scp.io.read(test_consensus_biorxiv_fig3b)
        sc2 = scp.io.read(test_consensus_biorxiv_figs18a)
        tris_tree = scp.tl.consensus(sc1, sc2)
        day_tree = scp.tl.consensus_day(sc1, sc2)
        assert nx.is_isomorphic(tris_tree, day_tree)

    def test_consensus_day_2(
        self, test_consensus_recomb_fig1a, test_consensus_recomb_fig1b
    ):
        """Verify isomorphism for the RECOMB figure 1 inputs."""
        # result in recomb.fig1c
        sc1 = scp.io.read(test_consensus_recomb_fig1a)
        sc2 = scp.io.read(test_consensus_recomb_fig1b)
        tris_tree = scp.tl.consensus(sc1, sc2)
        day_tree = scp.tl.consensus_day(sc1, sc2)
        assert nx.is_isomorphic(tris_tree, day_tree)

    def test_consensus_day_3(
        self, test_consensus_biorxiv_fig4b, test_consensus_biorxiv_fig4c
    ):
        """Verify isomorphism for the BioRxiv figure 4 inputs."""
        # result in biorxiv.fig4d
        sc1 = scp.io.read(test_consensus_biorxiv_fig4b)
        sc2 = scp.io.read(test_consensus_biorxiv_fig4c)
        tris_tree = scp.tl.consensus(sc1, sc2)
        day_tree = scp.tl.consensus_day(sc1, sc2)
        assert nx.is_isomorphic(tris_tree, day_tree)

    def test_consensus_day_4(
        self, test_consensus_biorxiv_fig3b, test_consensus_biorxiv_fig3c
    ):
        """Verify isomorphism for the BioRxiv figure 3 inputs."""
        # result in biorxiv.fig3f
        sc1 = scp.io.read(test_consensus_biorxiv_fig3b)
        sc2 = scp.io.read(test_consensus_biorxiv_fig3c)
        tris_tree = scp.tl.consensus(sc1, sc2)
        day_tree = scp.tl.consensus_day(sc1, sc2)
        assert nx.is_isomorphic(tris_tree, day_tree)
