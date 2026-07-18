import networkx as nx


def to_png(tree, filepath, dpi=150):
    """Render a phylogenetic tree as a PNG image."""
    mygraph = nx.drawing.nx_agraph.to_agraph(tree)
    mygraph.graph_attr["dpi"] = dpi
    mygraph.layout(prog="dot")
    mygraph.draw(filepath)
