"""Cover lightweight core helpers with small, deterministic fixtures."""

import builtins
import math
import os
import runpy
import subprocess
import sys
import threading
import time
import types
from pathlib import Path

import anndata as ad
import joblib
import networkx as nx
import numpy as np
import pandas as pd
import pytest

import scphylo as scp
import scphylo.io._genotype as genotype_io
import scphylo.logging._logging as logging_impl
import scphylo.ul._external as external
import scphylo.ul._servers as servers
import scphylo.ul._trees as tree_utils
import scphylo.ul._utils as general_utils


def test_io_roundtrips_vcf_and_validation(tmp_path, monkeypatch):
    """Exercise compact matrix roundtrips, VCF output, and input validation."""
    matrix = pd.DataFrame(
        [[0, 1, 3, 0], [1, 0, 0, 1]],
        index=["cell_a", "cell_b"],
        columns=["mut_a", "mut_b", "mut_c", "mut_d"],
    )
    tsv_path = tmp_path / "matrix.tsv"
    csv_path = tmp_path / "matrix.csv"
    scp.io.write(matrix.copy(), tsv_path)
    matrix.to_csv(csv_path)
    pd.testing.assert_frame_equal(scp.io.read(tsv_path), matrix, check_names=False)
    pd.testing.assert_frame_equal(scp.io.read(csv_path), matrix)

    adata = ad.AnnData(matrix.astype(int))
    adata.layers["total"] = np.array([[0, 5, 5, 10], [5, 5, 0, 10]])
    adata.layers["mutant"] = np.array([[0, 5, 0, 5], [0, 2, 0, 5]])
    output_stem = str(tmp_path / "counts")
    scp.io.write(adata, output_stem)
    assert scp.io.read(output_stem + ".h5ad.gz").shape == adata.shape

    monkeypatch.setattr(genotype_io.np, "math", math, raising=False)
    vcf_path = tmp_path / "matrix.vcf"
    scp.io.to_vcf(adata, vcf_path)
    lines = vcf_path.read_text().splitlines()
    assert lines[0] == "##fileformat=VCFv4.3"
    assert len([line for line in lines if line.startswith("chr1\t")]) == 4
    assert "0/0:0,0,0" in lines[-4]
    assert "0/1:0,0,255" in lines[-3]
    assert "./.:255,0,0" in lines[-2]

    duplicate = pd.DataFrame([[0, 1]], columns=["mut", "mut"])
    with monkeypatch.context() as patcher:
        patcher.setattr(genotype_io.pd, "read_table", lambda *args, **kwargs: duplicate)
        with pytest.raises(RuntimeError):
            scp.io.read(tmp_path / "duplicate.tsv")
    with monkeypatch.context() as patcher:
        patcher.setattr(
            genotype_io.pd,
            "read_csv",
            lambda *args, **kwargs: pd.DataFrame([[1]], columns=["mut"]),
        )
        with pytest.raises(RuntimeError):
            scp.io.read(tmp_path / "ones.csv")
    with pytest.raises(RuntimeError):
        scp.io.read(tmp_path / "matrix.unsupported")
    with pytest.raises(RuntimeError):
        scp.io.write(object(), tmp_path / "object.tsv")


def test_tree_png_renderer(tmp_path):
    """Render a two-node graph without using the plotting workflows."""
    tree = nx.DiGraph([(0, 1)])
    output = tmp_path / "tree.png"
    scp.io.to_png(tree, output, dpi=72)
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_dataset_and_noise_edge_paths(monkeypatch):
    """Cover the few loader and simulator paths absent from dataset tests."""
    assert scp.datasets.example(is_expression=True).n_obs > 0
    assert scp.datasets.acute_myeloid_leukemia1().shape == (1430, 15)
    assert scp.datasets.acute_myeloid_leukemia2().shape == (1066, 21)

    matrix = pd.DataFrame(
        [[0, 1], [1, 0], [0, 1], [1, 1]],
        index=["a", "b", "c", "d"],
        columns=["x", "y"],
    )
    flipped = scp.datasets.add_noise(matrix.iloc[:1], alpha=1, beta=1, missing=0)
    np.testing.assert_array_equal(flipped.to_numpy(), [[1, 0]])
    missing = scp.datasets.add_noise(matrix.iloc[:1], alpha=0, beta=0, missing=1)
    np.testing.assert_array_equal(missing.to_numpy(), [[3, 3]])

    monkeypatch.setattr(
        np.random, "choice", lambda values, **kwargs: np.array([values[0]])
    )
    doublets = scp.datasets.add_doublets(
        matrix, matrix.copy(), alpha=0, beta=0, missing=0, doublet=0.25
    )
    assert doublets.shape == matrix.shape

    with pytest.raises(RuntimeError):
        scp.datasets.add_noise(pd.DataFrame([[2]]), alpha=0, beta=0, missing=0)
    random_state = np.random.get_state()
    try:
        with pytest.raises(RuntimeError):
            scp.datasets.add_readcount(pd.DataFrame([[2]]), mean_coverage=2, seed=0)
    finally:
        np.random.set_state(random_state)


def test_logging_destinations_and_verbosity(tmp_path, monkeypatch, capsys):
    """Exercise file logging, colored output, and all verbosity representations."""
    logfile = tmp_path / "scphylo.log"
    monkeypatch.setattr(scp.settings, "verbosity", "debug")
    monkeypatch.setattr(scp.settings, "logfile", str(logfile))
    scp.logg.msg("default verbosity", v=None)
    scp.logg.hint("hint")
    scp.logg.info("timed", time=True)
    assert "default verbosity" in logfile.read_text()
    assert "--> hint" in logfile.read_text()

    monkeypatch.setattr(scp.settings, "logfile", "")
    scp.logg.msg("colored", v=2, color="red", end="!")
    scp.logg.msg(v=2)
    assert "colored" in capsys.readouterr().out

    monkeypatch.setattr(scp.settings, "verbosity", 0)
    scp.logg.info("suppressed")
    assert capsys.readouterr().out == ""
    assert len(logging_impl._get_date_string()) == 16


def test_preprocessing_groups_annotations_and_tree_sampling():
    """Cover grouping, SnpEff filtering, and mutation sampling on tiny objects."""
    adata = ad.AnnData(
        np.zeros((4, 3), dtype=int),
        obs=pd.DataFrame(
            {"group": ["a", "a", "b", "b"]},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["m1", "m2", "m3"]),
    )
    adata.layers["genotype"] = np.array([[0, 2, 0], [0, 2, 0], [1, 3, 0], [1, 2, 0]])
    adata.layers["mutant"] = np.arange(12).reshape(4, 3)
    separated = scp.pp.mut_seperated_by_cell_group(adata, "group")
    assert separated.loc["m1", "b"] == 1
    grouped = scp.pp.group_obs_apply_func(adata, "group", layer="mutant")
    np.testing.assert_array_equal(grouped["a"], [3, 5, 7])

    selected = adata[:3].copy()
    scp.pp.remove_nonsense_cells(selected)
    assert selected.obs_names.tolist() == ["c3"]

    annotated = ad.AnnData(
        np.zeros((2, 3)),
        obs=pd.DataFrame(index=["A_Tumor", "NB"]),
        var=pd.DataFrame(
            {
                "Transcript_BioType": ["protein_coding", "other", "protein_coding"],
                "Feature_Type": ["transcript", "transcript", "transcript"],
                "ALT": ["A", "C", "A,C"],
            },
            index=["v1", "v2", "v3"],
        ),
    )
    annotated.layers["mutant"] = np.array([[2, 3, 4], [0, 0, 0]])
    annotated.layers["total"] = np.array([[4, 6, 8], [1, 1, 1]])
    scp.pp.filter_snpeff(annotated, exome=True)
    assert annotated.shape == (1, 1)
    assert annotated.var.loc["v1", "VAF"] == 0.5
    assert annotated.var.loc["v1", "SAMPLE"] == "A_Tumor"

    matrix = pd.DataFrame(
        [[1, 1], [1, 0], [0, 0]],
        index=["c1", "c2", "normal"],
        columns=["m1", "m2"],
    )
    tree = scp.ul.to_tree(matrix)
    assert set(scp.pp.sample_from_tree(tree, 1, axis="mut")) == {"m1", "m2"}
    assert scp.pp.sample_from_tree(tree, 1, axis="unsupported") is None


def test_general_matrix_file_and_timeout_helpers(tmp_path, monkeypatch):
    """Exercise fast matrix statistics, parsers, temporary paths, and decorators."""
    input_values = np.array([[0, 1, 3, 3]])
    output_values = np.array([[1, 0, 0, 1]])
    assert scp.ul.count_flips(input_values, output_values) == (1, 1, 1, 1)
    assert scp.ul.infer_rates(input_values, output_values) == (1, 1, 0.5)

    valid = pd.DataFrame([[1, 1], [0, 1]])
    conflict = pd.DataFrame([[1, 1], [0, 1], [1, 0]])
    invalid = pd.DataFrame([[0, 3]])
    assert scp.ul.is_conflict_free(valid)
    assert not scp.ul.is_conflict_free(conflict)
    assert not scp.ul.is_conflict_free(invalid)
    assert scp.ul.calc_nll_matrix(valid, valid, 0, 0.1) is None
    assert np.isfinite(scp.ul.calc_nll_matrix(valid, valid, 0.01, 0.1))

    log_path = tmp_path / "run.tool.log"
    log_path.write_text(
        "output -- time: 1.5s (0:00:01.5)\n"
        "output -- CF: True\n"
        "rates -- FN: 0.125\n"
        "rates -- FP: 0.25\n"
    )
    assert scp.ul.parse_log_file(log_path) == {
        "tool": "tool",
        "running_time": 1.5,
        "is_cf": True,
        "fn_rate": 0.125,
        "fp_rate": 0.25,
    }
    score_path = tmp_path / "scores.demo.txt"
    score_path.write_text("first=1.5\nsecond=2\n")
    assert scp.ul.parse_score_file(score_path) == {
        "tool": "demo",
        "first": 1.5,
        "second": 2.0,
    }
    assert scp.ul.parse_params_file("plain.txt") == {}
    assert scp.ul.split_mut("ENSG_GENE.chr1.123.A.T") == (
        "ENSG",
        "GENE",
        "chr1",
        "123",
        "A",
        "T",
    )
    assert scp.ul.split_mut("invalid") == (None, None, None, None, None, None)

    directory = scp.ul.tmpdir(dirname=tmp_path)
    assert Path(directory).is_dir()
    scp.ul.cleanup(directory)
    temporary_file = scp.ul.tmpfile(dirname=tmp_path)
    assert Path(temporary_file).is_file()
    scp.ul.remove(temporary_file)
    scp.ul.remove(temporary_file)
    with scp.ul.tmpdirsys() as system_directory:
        assert Path(system_directory).is_dir()
    nested = tmp_path / "one" / "two"
    assert scp.ul.mkdir(nested) == nested
    assert scp.ul.mkdir(nested) == nested
    assert Path(scp.ul.get_file("scphylo.datasets/test/test.tsv")).is_file()

    @scp.ul.with_timeout(1)
    def immediate():
        return "done"

    event = threading.Event()

    @scp.ul.with_timeout(0.001)
    def delayed():
        return event.wait(0.05)

    assert immediate() == "done"
    assert delayed() is None

    clock = iter([0.0, 61.0])
    fake_time = types.SimpleNamespace(
        time=lambda: next(clock),
        gmtime=time.gmtime,
        strftime=time.strftime,
    )
    monkeypatch.setattr(general_utils, "time", fake_time)

    @scp.ul.timeit
    def timed():
        return 7

    assert timed() == 7


def test_logging_statistics_and_joblib_progress(monkeypatch):
    """Cover aggregate logging and Joblib progress callback restoration."""
    monkeypatch.setattr(scp.settings, "verbosity", 0)
    source = pd.DataFrame([[0, 1, 3, 3]], columns=list("abcd"))
    result = pd.DataFrame([[1, 0, 0, 1]], columns=list("abcd"))
    scp.ul.stat(source, result, alpha=0.01, beta=0.1, running_time=1.25)
    scp.ul.log_output(pd.DataFrame([[1, 1], [0, 1], [1, 0]]), 0)

    class Progress:
        def __init__(self):
            self.updates = 0
            self.closed = False

        def update(self, n):
            self.updates += n

        def close(self):
            self.closed = True

    progress = Progress()
    original = joblib.parallel.BatchCompletionCallBack
    with scp.ul.tqdm_joblib(progress):
        values = joblib.Parallel(n_jobs=2, prefer="threads")(
            joblib.delayed(abs)(value) for value in [-1, -2]
        )
    assert values == [1, 2]
    assert progress.updates == 2
    assert progress.closed
    assert joblib.parallel.BatchCompletionCallBack is original


def test_external_resolution_and_failure_translation(tmp_path, monkeypatch):
    """Cover path normalization and every translated subprocess failure."""
    assert external._directory_values(None) == []
    assert external._directory_values([tmp_path]) == [tmp_path]
    path_list = os.pathsep.join([str(tmp_path), str(tmp_path / "other")])
    assert external._directory_values(path_list) == [tmp_path, tmp_path / "other"]

    executable = tmp_path / "tool"
    executable.write_text("#!/bin/sh\n")
    executable.chmod(0o755)
    assert scp.ul.resolve_executable(executable, "Tool") == str(executable.resolve())
    artifact = tmp_path / "tool.jar"
    artifact.write_text("jar")
    assert scp.ul.resolve_external_file(artifact, "Tool") == str(artifact.resolve())
    assert external._log_tail(artifact, limit=2) == "ar"
    assert external._log_tail(tmp_path / "missing.log") == ""
    assert external._log_tail(None) == ""

    non_executable = tmp_path / "not-executable"
    non_executable.write_text("text")
    with pytest.raises(scp.ul.ExternalToolNotFoundError, match="tool_path"):
        scp.ul.resolve_executable(
            non_executable,
            "Tool",
            path_option="tool_path",
        )
    with pytest.raises(scp.ul.ExternalToolNotFoundError, match="explicit path"):
        scp.ul.resolve_external_file(tmp_path / "missing.jar", "Tool")
    with pytest.raises(ValueError, match="at least one"):
        scp.ul.run_external([], "Tool")

    failures = [
        (FileNotFoundError(), scp.ul.ExternalToolNotFoundError, "was not found"),
        (PermissionError(), scp.ul.ExternalToolExecutionError, "not executable"),
        (OSError("bad format"), scp.ul.ExternalToolExecutionError, "Cannot start"),
        (
            subprocess.TimeoutExpired(["tool"], 0.25),
            scp.ul.ExternalToolExecutionError,
            "timed out",
        ),
        (
            subprocess.CalledProcessError(9, ["tool"], stderr="stderr detail"),
            scp.ul.ExternalToolExecutionError,
            "exit code 9",
        ),
    ]
    for error, expected, match in failures:
        with monkeypatch.context() as patcher:
            patcher.setattr(
                external.subprocess,
                "run",
                lambda *args, error=error, **kwargs: (_ for _ in ()).throw(error),
            )
            with pytest.raises(expected, match=match):
                scp.ul.run_external(["tool"], "Tool")


def test_optional_import_helpers_without_optional_dependencies(monkeypatch):
    """Cover successful and missing optional imports using in-memory modules."""
    monkeypatch.setattr(scp.settings, "verbosity", 0)
    helpers = [
        ("gurobipy", scp.ul.import_gurobi),
        ("mpi4py", scp.ul.import_mpi4py),
        ("pygraphviz", scp.ul.import_graphviz),
        ("graph_tool", scp.ul.import_graph_tool),
    ]
    for module_name, helper in helpers:
        fake_module = types.ModuleType(module_name)
        with monkeypatch.context() as patcher:
            patcher.setitem(sys.modules, module_name, fake_module)
            assert helper() == (fake_module, False)
        real_import = builtins.__import__

        def blocked_import(
            name,
            *args,
            module_name=module_name,
            real_import=real_import,
            **kwargs,
        ):
            if name == module_name or name.startswith(module_name + "."):
                raise ImportError
            return real_import(name, *args, **kwargs)

        with monkeypatch.context() as patcher:
            patcher.setattr(builtins, "__import__", blocked_import)
            assert helper() == (None, True)

    fake_r = types.ModuleType("rpy2.robjects")
    fake_r.r = object()
    fake_callbacks = types.ModuleType("rpy2.rinterface_lib.callbacks")
    fake_callbacks.logger = types.SimpleNamespace(setLevel=lambda level: None)
    fake_packages = types.ModuleType("rpy2.robjects.packages")

    class PackageNotInstalledError(Exception):
        pass

    fake_packages.PackageNotInstalledError = PackageNotInstalledError
    sentinel = object()
    fake_packages.importr = lambda name: sentinel
    fake_modules = {
        "rpy2": types.ModuleType("rpy2"),
        "rpy2.rinterface_lib": types.ModuleType("rpy2.rinterface_lib"),
        "rpy2.rinterface_lib.callbacks": fake_callbacks,
        "rpy2.robjects": fake_r,
        "rpy2.robjects.packages": fake_packages,
    }
    with monkeypatch.context() as patcher:
        for name, module in fake_modules.items():
            patcher.setitem(sys.modules, name, module)
        assert scp.ul.import_rpy2("demo") == (sentinel, False)
        fake_packages.importr = lambda name: (_ for _ in ()).throw(
            PackageNotInstalledError
        )
        assert scp.ul.import_rpy2("demo", "install demo") == (None, True)

    real_import = builtins.__import__

    def block_rpy2(name, *args, **kwargs):
        if name == "rpy2" or name.startswith("rpy2."):
            raise ImportError
        return real_import(name, *args, **kwargs)

    with monkeypatch.context() as patcher:
        patcher.setattr(builtins, "__import__", block_rpy2)
        assert scp.ul.import_rpy2() == (None, True)


def test_sample_discovery_and_command_file_generation(tmp_path, monkeypatch):
    """Cover paired-read discovery and local and Biowulf command generation."""
    paired = tmp_path / "paired"
    paired.mkdir()
    for name in [
        "one_1.fastq.gz",
        "one_2.fastq.gz",
        "two_1.fastq.gz",
        "two_2.fastq.gz",
    ]:
        (paired / name).touch()
    assert scp.ul.is_paired(paired) == [
        {"sample": "one", "is_paired": True},
        {"sample": "two", "is_paired": True},
    ]

    single = tmp_path / "single"
    single.mkdir()
    (single / "sample.fastq.gz").touch()
    assert scp.ul.is_paired(single) == [{"sample": "sample", "is_paired": False}]
    for name in ["normal.markdup_bqsr.bam", "tumor.markdup_bqsr.bam"]:
        (single / name).touch()
    assert set(scp.ul.get_samples_df(single)["sample"]) == {"normal", "tumor"}
    assert scp.ul.get_samples_df(single, normal="normal")["sample"].tolist() == [
        "tumor"
    ]

    frame = pd.DataFrame({"cmd": ["first", "second"]})
    monkeypatch.setattr(servers.os, "uname", lambda: ("", "local", "", "", ""))
    local = servers.write_cmds_get_main(
        frame,
        "local",
        "1:00:00",
        2,
        None,
        1,
        "test@example.com",
        str(tmp_path / "local"),
    )
    assert "parallel" in local

    monkeypatch.setattr(servers.os, "uname", lambda: ("", "biowulf-node", "", "", ""))
    swarm = servers.write_cmds_get_main(
        frame,
        "sra",
        "1:00:00",
        2,
        "module/1",
        4,
        "test@example.com",
        str(tmp_path / "swarm"),
        afterok="123",
    )
    assert "swarm" in swarm and "--gres=lscratch:20" in swarm
    sbatch_dependency = servers.write_cmds_get_main(
        frame.iloc[:1],
        "one",
        "1:00:00",
        2,
        None,
        1,
        "test@example.com",
        str(tmp_path / "sbatch-dependency"),
        afterok="123",
    )
    assert "--dependency=afterok:123" in sbatch_dependency
    sbatch = servers.write_cmds_get_main(
        frame.iloc[:1],
        "plain",
        "1:00:00",
        2,
        None,
        1,
        "test@example.com",
        str(tmp_path / "sbatch"),
    )
    assert sbatch.endswith("plain.sh")
    assert servers.cmd(["one", "two"]) == "one two && "
    assert servers.cmd(["one", "two"], islast=True) == "one two"


def test_private_tree_conversions_and_empty_queries(capsys):
    """Cover private serializations, label splitting, and empty tree queries."""
    matrix = pd.DataFrame(
        [[1, 1], [1, 0], [0, 0]],
        index=["c1", "c2", "normal"],
        columns=["m1", "m2"],
    )
    tree = scp.ul.to_tree(matrix)
    newick = tree_utils._to_newick(tree)
    assert newick.endswith(";")
    assert "c1" in newick and "c2" in newick

    string_tree = nx.DiGraph([(2, 0), (2, 1)])
    nx.set_node_attributes(string_tree, {2: "root", 0: "a", 1: "b"}, "label")
    assert tree_utils._to_apted(string_tree) == "{root{a}{b}}"

    def mutation_tree(reverse=False):
        result = nx.DiGraph([(10, 1), (1, 2), (10, 3)])
        labels = {
            10: ["root"],
            1: ["b", "a"] if reverse else ["a", "b"],
            2: ["e"],
            3: ["d", "c"] if reverse else ["c", "d"],
        }
        nx.set_node_attributes(result, labels, "label")
        return result

    split_one, split_two = tree_utils._split_labels(
        mutation_tree(), mutation_tree(reverse=True)
    )
    for split_tree in [split_one, split_two]:
        root = scp.ul.root_id(split_tree)
        assert all(
            isinstance(split_tree.nodes[node]["label"], str)
            for node in split_tree
            if node != root
        )

    duplicate = nx.DiGraph([(0, 1), (1, 2)])
    duplicate.graph["splitter_mut"] = "\n"
    duplicate.graph["splitter_cell"] = "\n"
    nx.set_node_attributes(duplicate, {0: "root", 1: "––", 2: "cell"}, "label")
    nx.set_edge_attributes(duplicate, {(0, 1): "m", (1, 2): "m"}, "label")
    assert scp.ul.to_cfmatrix(duplicate) == ["m", "m"]
    assert capsys.readouterr().out.strip() == "m"

    partition_tree = nx.DiGraph([("a", "b"), ("a", "c")])
    partition_tree.graph["splitter_cell"] = ", "
    partition_tree.graph["data"] = pd.DataFrame(index=["a", "b", "c", "d"])
    inside, outside = scp.ul.partition_cells(partition_tree, ["a"])
    assert set(inside) == {"a", "b", "c"}
    assert outside.tolist() == ["d"]

    tree.graph["mutation_list"] = pd.DataFrame({"Node": ["[1]"]}, index=["m1"])
    assert scp.ul.cells_rooted_at(tree, "[missing]").size == 0
    assert scp.ul.muts_rooted_at(tree, "[missing]").size == 0
    assert scp.ul.is_leaf(
        tree, next(node for node in tree if tree.in_degree(node) == 1)
    )

    cycle = nx.DiGraph([(0, 1), (1, 0)])
    assert scp.ul.root_id(cycle) is None
    with pytest.raises(RuntimeError, match=""):
        scp.ul.to_tree(pd.DataFrame([[1, 1], [0, 1], [1, 0]]))


def test_star_rsem_and_bamreadcount_scripts(tmp_path, monkeypatch):
    """Run three data-collection modules on tiny temporary inputs."""
    star_dir = tmp_path / "star"
    star_dir.mkdir()
    star_lines = [
        "Number of input reads | 100",
        "Uniquely mapped reads number | 80",
        "Uniquely mapped reads % | 80%",
        "Number of splices: Total | 10",
        "Number of splices: GC/AG | 1",
        "Insertion average length | 1.5",
        "Deletion average length | 2.5",
        "% of reads unmapped: too short | 3%",
        "Average mapped length | 75",
        "Deletion rate per base | 0.1%",
        "Mismatch rate per base, % | 0.2%",
        "Average input read length | 100",
        "Number of splices: AT/AC | 2",
        "Number of splices: Annotated (sjdb) | 7",
        "Number of splices: GT/AG | 6",
        "Number of reads mapped to too many loci | 2",
        "Number of reads unmapped: too many mismatches | 3",
        "% of reads unmapped: too many mismatches | 3%",
        "Number of reads unmapped: other | 4",
        "Insertion rate per base | 0.3%",
        "% of reads unmapped: other | 4%",
        "% of reads mapped to multiple loci | 5%",
        "Number of reads mapped to multiple loci | 5",
        "Number of splices: Non-canonical | 1",
        "Number of reads unmapped: too short | 3",
        "% of reads mapped to too many loci | 2%",
    ]
    (star_dir / "sample.star.log").write_text("\n".join(star_lines) + "\n")
    monkeypatch.setattr(sys, "argv", ["star", str(star_dir)])
    star_globals = runpy.run_module("scphylo.ul.star", run_name="__main__")
    star_data = star_globals["get_numreads_percmapped"](star_dir / "sample.star.log")
    assert all(value is not None for value in star_data.values())
    assert (star_dir / "_starinfo.csv").is_file()

    rsem_dir = tmp_path / "rsem"
    rsem_dir.mkdir()
    (rsem_dir / "sample.transcript.bam").touch()
    pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "expected_count": [1, 2],
            "TPM": [3, 4],
            "FPKM": [5, 6],
        }
    ).to_csv(rsem_dir / "sample.genes.results", sep="\t", index=False)
    monkeypatch.setattr(sys, "argv", ["rsem", str(rsem_dir)])
    runpy.run_module("scphylo.ul.rsem", run_name="__main__")
    expression = ad.read_h5ad(rsem_dir / "_expression.h5ad.gz")
    assert expression.shape == (1, 2)
    assert {"TPM", "FPKM"} <= set(expression.layers)

    readcount_dir = tmp_path / "readcount"
    readcount_dir.mkdir()
    lines = [
        "chr1\t1\tA\t0\tA:0\tC:0\tG:0\tT:0\tN:0",
        "chr1\t2\tA\t5\tA:5\tC:0\tG:0\tT:0\tN:0",
        "chr1\t3\tA\t5\tA:2\tC:3\tG:0\tT:0\tN:0",
    ]
    (readcount_dir / "cell.bamreadcount").write_text("\n".join(lines) + "\n")
    vcf_path = readcount_dir / "variants.csv"
    pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1", "chr1"],
            "POS": [1, 2, 3],
            "Allele": ["C", "C", "C"],
        }
    ).to_csv(vcf_path, index=False)
    monkeypatch.setattr(
        sys, "argv", ["bamreadcount", str(readcount_dir), str(vcf_path)]
    )
    bam_globals = runpy.run_module("scphylo.ul.bamreadcount", run_name="__main__")
    assert bam_globals["is_mut"](pd.Series({"COV": 0, "Allele": "C", "C": 0})) == (
        3,
        0,
    )
    readcounts = ad.read_h5ad(readcount_dir / "_bamreadcount.h5ad.gz")
    np.testing.assert_array_equal(readcounts.X, [[3, 0, 1]])


def test_cyvcf_conversion_script_with_in_memory_reader(tmp_path, monkeypatch):
    """Run the annotated-VCF converter against a two-sample fake reader."""
    annotation_fields = [
        "Allele",
        "Annotation",
        "Annotation_Impact",
        "Gene_Name",
        "Gene_ID",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "ERRORS / WARNINGS / INFO",
    ]
    annotation = "|".join(
        [
            "C",
            "missense_variant",
            "MODERATE",
            "GENE",
            "ENSG",
            "transcript",
            "ENST",
            "protein_coding",
            "1/1",
            "c.1A>C",
            "p.K1Q",
            "1/1",
            "1/1",
            "1/1",
            "0",
            "",
        ]
    )

    class Variant:
        def __init__(self, position, ann, *, filtered=False, variant=True):
            self.FILTER = "q10" if filtered else None
            self.is_snp = variant
            self.is_indel = False
            self.CHROM = "chr1"
            self.POS = position
            self.REF = "A"
            self.ALT = "C"
            self.start = position - 1
            self.end = position
            self.INFO = {"ANN": ann}
            self.gt_types = np.array([0, 1])
            self.gt_ref_depths = np.array([5, 3])
            self.gt_alt_depths = np.array([0, 2])

    class FakeVCF:
        def __init__(self, infile):
            self.samples = ["normal", "tumor"]

        def get_header_type(self, name):
            return {"Description": " | ".join(annotation_fields)}

        def __iter__(self):
            return iter(
                [
                    Variant(1, annotation),
                    Variant(2, None),
                    Variant(3, annotation, filtered=True),
                    Variant(4, annotation, variant=False),
                ]
            )

    fake_cyvcf2 = types.ModuleType("cyvcf2")
    fake_cyvcf2.VCF = FakeVCF
    monkeypatch.setitem(sys.modules, "cyvcf2", fake_cyvcf2)
    input_path = tmp_path / "sample.ann.vcf"
    input_path.touch()
    monkeypatch.setattr(sys, "argv", ["cyvcf", str(input_path)])
    runpy.run_module("scphylo.ul.cyvcf", run_name="__main__")

    converted = ad.read_h5ad(tmp_path / "_sample.h5ad.gz")
    assert converted.shape == (2, 2)
    assert {"genotype", "total", "mutant"} <= set(converted.layers)
    assert {"normal_TOTAL_READS", "tumor_MUTANT_READS"} <= set(converted.var)
    assert (tmp_path / "sample.ann.tsv").is_file()
