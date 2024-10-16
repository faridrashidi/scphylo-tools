import networkx as nx
import numpy as np

import scphylo as scp


def collapse(tree):
    tc2 = tree.copy()
    root = scp.ul.root_id(tree)
    nnodes = 1
    for _ in range(len(tree.nodes)):
        d_in = tc2.in_degree(tc2)
        d_out = tc2.out_degree(tc2)
        for node in tc2.nodes():
            if d_out[node] == 1 and d_in[node] == 1:
                parent = list(x for x in tc2.predecessors(node))[0]
                child = list(x for x in tc2.successors(node))[0]
                if d_out[parent] < 2 and d_in[parent] == 1:
                    new_node = root + nnodes
                    nnodes += 1
                    new_label = (
                        f"{tc2.nodes[parent]['label']}"
                        f"{tree.graph['splitter_cell']}"
                        f"{tc2.nodes[node]['label']}"
                    )
                    new_edge = (
                        f"{tc2[parent][node]['label']}"
                        f"{tree.graph['splitter_mut']}"
                        f"{tc2[node][child]['label']}"
                    )

                    tc2 = nx.contracted_nodes(tc2, parent, node, self_loops=False)
                    mapping = {parent: new_node}
                    tc2 = nx.relabel_nodes(tc2, mapping)
                    tc2[new_node][child]["label"] = new_edge
                    tc2.nodes[new_node]["label"] = new_label
                    break
    d_in = tc2.in_degree(tc2)
    d_out = tc2.out_degree(tc2)
    nodes = []
    for node in tc2.nodes():
        if d_out[node] == 0:
            nodes.append(node)
    for node in nodes:
        parent = list(x for x in tc2.predecessors(node))[0]
        if d_out[parent] == 1 and d_in[parent] == 1:
            grandparent = list(x for x in tc2.predecessors(parent))[0]

            new_node = root + nnodes
            nnodes += 1
            new_label = (
                f"{tc2.nodes[parent]['label']}"
                f"{tree.graph['splitter_cell']}"
                f"{tc2.nodes[node]['label']}"
            )
            new_edge = (
                f"{tc2[grandparent][parent]['label']}"
                f"{tree.graph['splitter_mut']}"
                f"{tc2[parent][node]['label']}"
            )

            tc2 = nx.contracted_nodes(tc2, parent, node, self_loops=False)
            mapping = {parent: new_node}
            tc2 = nx.relabel_nodes(tc2, mapping)
            tc2[grandparent][new_node]["label"] = new_edge
            tc2.nodes[new_node]["label"] = new_label
    return tc2


def sample_from_tree(tree, ratio, axis="cell"):
    sampled = []
    if axis == "cell":
        for n in tree.nodes:
            if tree.in_degree(n) == 1 and "––" not in tree.nodes[n]["label"]:
                array = tree.nodes[n]["label"].split(tree.graph["splitter_cell"])
                array = np.random.choice(
                    array, replace=False, size=int(ratio * len(array))
                )
                sampled += list(array)
        return sampled
    elif axis == "mut":
        for _, _, lable in tree.edges(data=True):
            array = lable["label"].split(tree.graph["splitter_mut"])
            array = np.random.choice(array, replace=False, size=int(ratio * len(array)))
            sampled += list(array)
        return sampled
    else:
        return None
