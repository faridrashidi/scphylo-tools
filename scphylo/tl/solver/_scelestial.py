import math
import os
import sys
import time

import networkx as nx
import numpy as np
import pandas as pd

import scphylo as scp


def scelestial(df_input):
    """Solving using Scelestial.

    Fast and accurate single-cell lineage tree inference based on a Steiner tree
    approximation algorithm :cite:`Scelestial`.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Input genotype matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1), absence (0) and missing
        entires (3).

    Returns
    -------
    :class:`pandas.DataFrame`
        A conflict-free matrix in which rows are cells and columns are mutations.
        Values inside this matrix show the presence (1) and absence (0).
    """
    executable = scp.ul.executable("scelestial", "Scelestial")

    scp.logg.info("running Scelestial")

    tmpdir = scp.ul.tmpdirsys(suffix=".scelestial")

    np.savetxt(
        f"{tmpdir.name}/scelestial.SC.T", df_input.values.T, delimiter="\t", fmt="%1.0f"
    )
    with open(f"{tmpdir.name}/scelestial.cellNames", "w") as fout:
        fout.write("\n".join(df_input.index))
    with open(f"{tmpdir.name}/scelestial.mutNames", "w") as fout:
        fout.write("\n".join(df_input.columns))

    _convert_input(
        f"{tmpdir.name}/scelestial.SC.T", f"{tmpdir.name}/scelestial.input", "/dev/null"
    )

    cmd = (
        f"{executable} < {tmpdir.name}/scelestial.input > "
        f"{tmpdir.name}/scelestial.tree_clone"
    )

    s_time = time.time()
    os.system(cmd)
    e_time = time.time()
    running_time = e_time - s_time

    _steiner_to_seq(
        f"{tmpdir.name}/scelestial.tree_clone",
        f"{tmpdir.name}/scelestial.impute",
    )

    _stein_to_clone_tree(
        f"{tmpdir.name}/scelestial.input",
        f"{tmpdir.name}/scelestial.tree_clone",
        f"{tmpdir.name}/scelestial.tree",
        f"{tmpdir.name}/scelestial.clone",
    )

    _clone_tree_to_mu_tree_imput(
        f"{tmpdir.name}/scelestial.tree",
        f"{tmpdir.name}/scelestial.clone",
        f"{tmpdir.name}/scelestial.SC.T",
        f"{tmpdir.name}/scelestial.mutNames",
        f"{tmpdir.name}/scelestial.cellNames",
        f"{tmpdir.name}/scelestial.output",
    )

    tree = nx.DiGraph(nx.nx_pydot.read_dot(f"{tmpdir.name}/scelestial.output"))
    tree = nx.relabel_nodes(tree, lambda x: int(x))
    tree.graph["splitter_mut"] = "\n"
    tree.graph["splitter_cell"] = "\n"
    tree.add_node(-1, label="root")
    tree.add_edge(-1, scp.ul.root_id(tree))
    for i, e in enumerate(tree.edges()):
        tree.edges[e]["label"] = f"mut{i}"
        tree.edges[e]["weight"] = 1
    for n in tree.nodes():
        if tree.nodes[n]["label"] == '""':
            tree.nodes[n]["label"] = "––"
    df = scp.ul.to_cfmatrix(tree)
    df2 = pd.DataFrame(
        np.identity(df.shape[0]),
        index=df.index,
        columns=[f"mut{x}" for x in range(df.shape[0], 2 * df.shape[0])],
        dtype=int,
    )
    df_output = df.merge(df2, right_index=True, left_index=True)

    tmpdir.cleanup()

    scp.ul.stat(df_output, df_output, 0, 0, running_time)

    return df_output


def _convert_input(sciteFile_path, imputeFile_path, bpFile_path):
    sciteFile = open(sciteFile_path)
    imputeFile = open(imputeFile_path, "w+")
    bpFile = open(bpFile_path, "w+")

    seq = []
    for line in sciteFile:
        x = line.strip().split()
        for i, v in enumerate(x):
            while i >= len(seq):
                seq.append([])
            seq[i].append(int(v))

    loc = len(seq[0])
    for i in range(loc):
        print(f"{i + 1}", end=" ", file=imputeFile)
        for j in range(len(seq)):
            s = seq[j][i]
            if s == 0:
                ss = "A/A"
            elif s == 1:
                ss = "C/C"
            elif s == 2:
                ss = "A/C"
            elif s == 3:
                ss = "./."
            else:
                print(f"OH! {s}")
                raise ValueError(f"OH! {s}")
            print(f"{ss}", end=" ", file=imputeFile)
        print(file=imputeFile)

    for i in range(loc):
        print(f"V{i + 1}", end="", file=bpFile)
        if i != loc - 1:
            print(",", end="", file=bpFile)
    print(file=bpFile)

    for s in seq:
        first = True
        for v in s:
            if v == 3:
                v = 2
            elif v == 2:
                v = 1
            if not first:
                print(",", end="", file=bpFile)
            first = False
            print(f"{v}", end="", file=bpFile)
        print(file=bpFile)

    sciteFile.close()
    imputeFile.close()
    bpFile.close()


def _stein_to_clone_tree(
    seqSciteFile_path, steinerFile_path, treeFile_path, cloneFile_path
):
    seqSciteFile = open(seqSciteFile_path)
    steinerFile = open(steinerFile_path)
    treeFile = open(treeFile_path, "w+")
    cloneFile = open(cloneFile_path, "w+")
    type(seqSciteFile)

    n = int(steinerFile.readline().strip())
    treeNodes = []
    cells = []
    maxAc = -1
    for _ in range(n):
        line = steinerFile.readline().strip()
        x = line.split()
        treeNodes.append(x[0])
        if x[1] == "1":
            cells.append(x[0])
        ac = sum(1 for s in x[2] if s == "A")
        if ac > maxAc:
            maxAc = ac
            maxAcLen = len(x[2])
            treeRootSeqIdx = int(x[0])

    seqIdxToExclude = -1
    if maxAcLen == maxAc and 5 < len(sys.argv) and sys.argv[5] == "-exclude-root":
        seqIdxToExclude = treeRootSeqIdx

    cells = set(cells)

    for t in treeNodes:
        print(t, end=" ", file=cloneFile)
        if t in cells and t != seqIdxToExclude and int(t) != seqIdxToExclude:
            print(int(t) + 1, end="", file=cloneFile)
        if int(t) == treeRootSeqIdx:
            treeRootCloneIdx = int(t)
        print(file=cloneFile)

    print(" ".join(treeNodes), file=treeFile)

    edges = {}

    m = int(steinerFile.readline().strip())
    for _ in range(m):
        line = steinerFile.readline().strip()
        x = line.split()
        if x[0] not in edges:
            edges[x[0]] = []
        if x[1] not in edges:
            edges[x[1]] = []
        edges[x[0]].append((x[1], x[2]))
        edges[x[1]].append((x[0], x[2]))

    mark = {}
    global dfsNumCounter
    dfsNumCounter = 0
    global dfsNum
    dfsNum = {}

    def dfs(v):
        global dfsNumCounter, dfsNum
        dfsNum[v] = dfsNumCounter
        dfsNumCounter += 1
        mark[v] = True
        for (u, _) in edges[v]:
            if u not in mark:
                dfs(u)

    dfs(str(treeRootCloneIdx))

    for v, nei in edges.items():
        for (u, w) in nei:
            if v < u:
                x, y = u, v
                if dfsNum[v] > dfsNum[u]:
                    x, y = v, u
                print(f"{x}->{y} {w}", file=treeFile)

    seqSciteFile.close()
    steinerFile.close()
    treeFile.close()
    cloneFile.close()


def _clone_tree_to_mu_tree_imput(
    treeFileName,
    cloneFileName,
    seqFileName,
    mutationInfoFileName,
    cellNamesFileName,
    outputFileName,
    margeClones=False,
    markMutations=False,
    compress=False,
    markMutationsSeparated=False,
):
    from graphviz import Digraph
    from scipy import stats

    def loadTree(treeFileName):
        treeFile = open(treeFileName)
        line = treeFile.readline()
        vertices = line.strip().split()
        edges = []
        treeParent = {}
        treeChildren = {}
        for line in treeFile:
            (e, w) = line.strip().split(" ")
            (v, u) = e.split("->")
            edges.append((v, u, float(w)))

            treeParent[v] = u
            if u not in treeChildren:
                treeChildren[u] = []
            treeChildren[u].append((v, float(w)))

        treeRoot = list(treeParent.keys())[0]
        while treeRoot in treeParent:
            treeRoot = treeParent[treeRoot]

        treeFile.close()
        return vertices, edges, treeParent, treeChildren, treeRoot

    def loadClones(cloneFileName):
        cloneFile = open(cloneFileName)
        treeNodeCells = {}

        for line in cloneFile:
            x = line.strip().split()
            treeNodeCells[x[0]] = x[1:]

        cloneFile.close()
        return treeNodeCells

    def compressedTree(treeChildren, treeNodeCells, treeRoot):
        compressedTreeChildren = {}

        def dfs(v):
            ret = (None, None)
            children = []
            if v in treeChildren:
                for u, cw in treeChildren[v]:
                    c, w = dfs(u)
                    if c is not None:
                        w += cw
                        children.append((c, w))
                        ret = (c, w)

            if (
                (v in treeNodeCells and len(treeNodeCells[v]) > 0)
                or len(children) > 1
                or v == treeRoot
            ):
                ret = (v, 0)

            if ret == (v, 0):
                compressedTreeChildren[v] = []
                for c, cw in children:
                    compressedTreeChildren[v].append((c, cw))

            return ret

        dfs(treeRoot)
        return compressedTreeChildren

    def loadSequenceFile(seqFileName):
        sequences = []
        seqFile = open(seqFileName)
        for line in seqFile:
            for num, val in enumerate(line.strip().split()):
                while num >= len(sequences):
                    sequences.append([])
                sequences[num].append(int(val))
        seqFile.close()
        return sequences

    def loadMutationInfoFile(mutationInfoFileName):
        mutationInfoFile = open(mutationInfoFileName)
        mutationInfo = []
        for line in mutationInfoFile:
            mutationInfo.append({id: line.strip()})
        mutationInfoFile.close()
        return mutationInfo

    def writeGraph(
        outputFileName, treeNodeCells, nodes, edges, treeNodeDescColor, treeEdgeLabel
    ):
        dot = Digraph(format="pdf")
        dot.graph_attr["rankdir"] = "LR"
        for treeNode in nodes:
            cells = []
            if treeNode in treeNodeCells:
                cells = treeNodeCells[treeNode]
            prop = treeNodeDescColor(treeNode, cells)
            if len(prop) == 2:
                desc, col, fontcol, fillcol = prop[0], prop[1], "black", "none"
            else:
                desc, col, fontcol, fillcol = prop[0], prop[1], prop[2], prop[3]
            dot.node(
                treeNode,
                desc,
                color=col,
                fillcolor=fillcol,
                style="filled",
                fontcolor=fontcol,
                gradientangle="0",
                penwidth="4",
                shape="circle",
                margin="0",
            )

        for v, u, w in edges:
            tup = treeEdgeLabel(v, u, w)
            if isinstance(tup, tuple):
                label, edgecol = tup[0], tup[1]
            else:
                label, edgecol = tup, "black"
            dot.edge(u, v, weight=str(w), label=label, color=edgecol)

        dot.render(outputFileName)

    def loadFileRows(fileName):
        f = open(fileName)
        r = []
        for line in f:
            r.append(line.strip())
        f.close()
        return r

    def loadCellNames(cellNamesFileName):
        return loadFileRows(cellNamesFileName)

    BLUE = "#69c5f0"
    BLACK = "black"

    _, edges, _, treeChildren, treeRoot = loadTree(treeFileName)
    treeNodeCells = loadClones(cloneFileName)

    sequences = loadSequenceFile(seqFileName)
    mutationInfo = loadMutationInfoFile(mutationInfoFileName)
    cellNames = loadCellNames(cellNamesFileName)

    if margeClones:

        def sameColon(v, u, w):
            m = len(sequences[0])
            pv = stats.binom_test(w, m, 0.2, alternative="less")
            print(f"same colon test: {w}/{m} = {pv}")
            return pv < 0.01

        mergeParent = {v: v for v in treeChildren.keys()}
        toBeMerged = {}
        for v, childrenDistList in treeChildren.items():
            for c, w in childrenDistList:
                if sameColon(v, c, w):
                    mergeParent[c] = v
                    toBeMerged[c] = True

        newTreeChildren = {v: [] for v in treeChildren.keys()}
        newTreeNodeCells = {v: [] for v in treeNodeCells.keys()}

        def dfs(v, firstKeptNode, wToHere):
            newTreeNodeCells[firstKeptNode] += treeNodeCells[v]
            if v in treeChildren:
                for u, w in treeChildren[v]:
                    if u not in toBeMerged:
                        newTreeChildren[firstKeptNode].append((u, wToHere + w))
                        dfs(u, u, 0)
                    else:
                        dfs(u, firstKeptNode, wToHere + w)

        dfs(treeRoot, treeRoot, 0)

        for v, _ in treeChildren.items():
            while mergeParent[mergeParent[v]] != mergeParent[v]:
                mergeParent[v] = mergeParent[mergeParent[v]]

        treeNodeCells = newTreeNodeCells
        treeChildren = newTreeChildren

    treeNodeMutations = {}

    def fillTreeNodeMutations(mutIndex=None):
        def dfs(v):
            myCellsStar = []
            if v in treeNodeCells:
                myCellsStar += [int(u) - 1 for u in treeNodeCells[v]]
            if v in treeChildren:
                for u, _ in treeChildren[v]:
                    myCellsStar = dfs(u) + myCellsStar

            nodeMutations[v] = []
            for i, _ in enumerate(mutationInfo):
                subTreeNormal, subTreeMutated = 0, 0
                for c in myCellsStar:
                    if sequences[c][i] == 0:
                        subTreeNormal += 1
                    if sequences[c][i] == 1:
                        subTreeMutated += 1

                if subTreeMutated > 0:
                    ignoremut = False
                    if v in treeChildren:
                        for u, _ in treeChildren[v]:
                            smi = [
                                submut for ii, submut, _ in nodeMutations[u] if i == ii
                            ]
                            if len(smi) == 1 and smi[0] == subTreeMutated:
                                ignoremut = True
                    if not ignoremut:
                        nodeMutations[v].append(
                            (
                                i,
                                subTreeMutated,
                                str(subTreeNormal) + "," + str(len(myCellsStar)),
                            )
                        )
            return myCellsStar

        nodeMutations = {}
        dfs(treeRoot)
        for v, muts in nodeMutations.items():
            treeNodeMutations[v] = [
                mutationInfo[i][id] + "/" + str(mut) + "," + str(desc)
                for i, mut, desc in muts
                if mutIndex is None or i == mutIndex
            ]

    if markMutations:
        fillTreeNodeMutations()

    nodes = list(treeNodeCells.keys())
    if compress:
        compressedTreeChildren = compressedTree(treeChildren, treeNodeCells, treeRoot)
        edges = []
        for par, cwList in compressedTreeChildren.items():
            for c, w in cwList:
                edges.append((c, par, w))
        nodes = list(compressedTreeChildren.keys())

    def treeNodeDescColor(treeNode, cells):
        desc = ", ".join([cellNames[int(c) - 1] for c in cells])
        if treeNodeMutations is not None and treeNode in treeNodeMutations:
            mutations = sorted(set(treeNodeMutations[treeNode]))
            seqmut = math.sqrt(len(mutations) / 7)
            if len(mutations) > 0:
                desc += " ["
                lastLineLen = 0
                for _, mut in enumerate(mutations):
                    desc += mut + " "
                    lastLineLen += 1
                    if lastLineLen >= seqmut:
                        desc += "\n"
                        lastLineLen = 0
                desc += "]"
        return desc, BLUE, BLACK, BLUE

    def treeEdgeLabel(v, u, w):
        return "", "#A2A2A2"

    if markMutationsSeparated:
        for i, _ in enumerate(mutationInfo):
            fillTreeNodeMutations(i)
            writeGraph(
                outputFileName + "-" + str(i),
                treeNodeCells,
                nodes,
                edges,
                treeNodeDescColor,
                treeEdgeLabel,
            )
    else:
        writeGraph(
            outputFileName,
            treeNodeCells,
            nodes,
            edges,
            treeNodeDescColor,
            treeEdgeLabel,
        )


def _steiner_to_seq(steinerFile_path, imputeFile_path):
    steinerFile = open(steinerFile_path)
    imputeFile = open(imputeFile_path, "w+")
    n = int(steinerFile.readline().strip())
    treeNodes = []
    cells = []
    maxAc = -1
    seq = []
    for _ in range(n):
        line = steinerFile.readline().strip()
        x = line.split()
        treeNodes.append(x[0])
        if x[1] == "1":
            cells.append(x[0])
        ac = sum(1 for s in x[2] if s == "A")
        if ac >= maxAc:
            maxAc = ac

        if x[1] == "1":
            if sum(1 for s in x[3] if s not in {"A", "C"}) > 0:
                print(
                    f"imputed sequence contains invalid char {x[3]}",
                    file=sys.stderr,
                )
                raise Exception("Invalid imputation")
            seq.append(x[3])

    if len(seq) == 0:
        exit()
    for j in range(len(seq[0])):
        print(
            " ".join(["0" if s[j] == "A" else "1" for s in seq]),
            file=imputeFile,
        )
    steinerFile.close()
    imputeFile.close()
