"""Cover plotting branches with rendering backends replaced by in-memory fakes."""

from contextlib import nullcontext
from types import SimpleNamespace
from unittest.mock import Mock

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd

from scphylo.pl import _data, _helper, _trees


def test_heatmap_alternate_inputs_and_distance_histogram(monkeypatch):
    """Cover row colors, alternate layers, raw frames, and distance plotting."""
    cluster = Mock()
    histogram = Mock()
    monkeypatch.setattr(_data.sns, "clustermap", cluster)
    monkeypatch.setattr(_data.plt, "hist", histogram)
    adata = ad.AnnData(
        np.zeros((2, 2)),
        obs=pd.DataFrame({"group": ["a", "b"]}, index=["c1", "c2"]),
    )
    adata.obsm["custom"] = np.ones((2, 2))
    adata.obsp["distance"] = np.array([[0.0, 2.0], [2.0, 0.0]])

    _data.heatmap(adata, color_attrs=["group"], layer="custom", rvb=["white", "red"])
    _data.heatmap(pd.DataFrame([[0, 1]]), color_attrs=None, rvb=["white", "black"])
    _data.plot_dist(adata, "distance")

    assert cluster.call_count == 2
    np.testing.assert_array_equal(histogram.call_args.args[0], [2.0])


def test_plot_helpers_build_chromosome_and_structured_mutation_metadata():
    """Cover chromosome command generation and parsed mutation annotations."""
    assert "chr1" in _helper._add_chromplot()
    assert "chrX" in _helper._add_chromplot_helper("X", 10)
    matrix = pd.DataFrame(
        [[1], [0]],
        index=["cell1", "cell2"],
        columns=["ENSG_GENE.chr1.123.A.T"],
    )
    tree = _helper.scp.ul.to_tree(matrix)

    newick, info, mutations = _helper._newick_info2_mutation_list(tree)

    assert newick.endswith(";")
    assert "root" in info["nmuts_label"].values
    assert mutations.loc[0, "Gene"] == "GENE"


def test_networkx_tree_rendering_attribute_modes(monkeypatch):
    """Cover ID and custom node/edge attribute rendering without Graphviz."""
    rendered = []

    class FakePydot:
        def create_png(self):
            rendered.append(True)
            return b"png"

    monkeypatch.setattr(
        _trees.nx.drawing.nx_pydot, "to_pydot", Mock(return_value=FakePydot())
    )
    monkeypatch.setattr(_trees, "Image", Mock(return_value=object()))
    monkeypatch.setattr(_trees, "display", Mock())
    tree = nx.DiGraph([(0, 1), (1, 2)])
    nx.set_node_attributes(tree, {0: "root", 1: "inner", 2: "leaf"}, "name")
    nx.set_edge_attributes(tree, {(0, 1): "a", (1, 2): "b"}, "mutation")

    _trees.networkx_tree(tree, n_attr="id", e_attr="mutation")
    _trees.networkx_tree(tree, n_attr="name", e_attr="mutation")

    assert len(rendered) == 2


def test_newick_tree_uses_mocked_r_rendering(monkeypatch):
    """Exercise both leaf-color and tip-label paths without requiring ggtree."""
    import rpy2.robjects as ro
    import rpy2.robjects.lib.grdevices as grdevices
    import rpy2.robjects.packages as packages

    class FakeR:
        def __call__(self, command):
            return command

        def show(self, _plot):
            return None

    class FakeImage:
        def getvalue(self):
            return b"png"

    monkeypatch.setattr(
        _trees.scp.ul, "import_rpy2", Mock(return_value=(object(), False))
    )
    monkeypatch.setattr(packages, "importr", Mock(return_value=object()))
    monkeypatch.setattr(ro, "r", FakeR())
    monkeypatch.setattr(
        grdevices,
        "render_to_bytesio",
        lambda *_args, **_kwargs: nullcontext(FakeImage()),
    )
    monkeypatch.setattr(_trees, "Image", Mock(return_value=object()))
    monkeypatch.setattr(_trees, "display", Mock(return_value="shown"))
    adata = SimpleNamespace(
        obs=pd.DataFrame({"group_color": ["red"]}, index=["cell1"]),
        uns={"tree": "(cell1);"},
    )

    assert (
        _trees.newick_tree(
            adata, leaf_attr="group_color", show_tiplab=True, width=10, height=10
        )
        == "shown"
    )
    assert _trees.newick_tree(adata, leaf_attr=None, width=10, height=10) == "shown"
