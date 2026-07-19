"""Fast unit tests for the command-line adapter modules."""

from importlib import import_module
from unittest.mock import Mock

import pandas as pd
import pytest

import scphylo as scp


@pytest.fixture(autouse=True)
def _restore_global_settings(monkeypatch):
    """Prevent command callbacks from leaking logging settings to later tests."""
    monkeypatch.setattr(scp.settings, "verbosity", scp.settings.verbosity)
    monkeypatch.setattr(scp.settings, "logfile", scp.settings.logfile)


def _invoke(command, *args):
    """Invoke a Click callback after narrowing its optional type."""
    callback = command.callback
    assert callback is not None
    return callback(*args)


@pytest.mark.parametrize(
    ("package", "expected"),
    [
        (
            "scphylo.commands.caller",
            [
                "rename",
                "sra",
                "fastp",
                "bwa",
                "star",
                "gatk",
                "hapcaller",
                "mutect2",
                "strelka2",
                "varscan2",
                "mpileup",
                "sequenza",
                "varfilter",
                "snpeff",
                "rsem",
                "bamreadcount",
                "bamsnap",
                "defuse",
                "vartrix",
                "pyclone",
            ],
        ),
        (
            "scphylo.commands.solver",
            [
                "booster",
                "phiscsb",
                "phiscsi",
                "scite",
                "scistree",
                "bnb",
                "huntress",
                "grmt",
                "sphyr",
                "gpps",
                "onconem",
                "siclonefit",
            ],
        ),
        ("scphylo.commands.methyl", ["build-crcs"]),
    ],
)
def test_command_groups_preserve_declared_order(package, expected):
    """Import every command and keep their intentionally curated help order."""
    module = import_module(package)
    group = next(
        value for name, value in vars(module).items() if name.startswith("cli_")
    )

    assert list(group.list_commands(None)) == expected
    assert _invoke(group) is None


@pytest.mark.parametrize(
    ("module_name", "command_name", "solver_name", "arguments", "expected"),
    [
        ("_bnb", "bnb", "bnb", ("real", 9), {"bounding": "real", "time_limit": 9}),
        (
            "_grmt",
            "grmt",
            "grmt",
            (0.1, 0.2, 3, 4),
            {"alpha": 0.1, "beta": 0.2, "n_iters": 3, "n_threads": 4},
        ),
        (
            "_huntress",
            "huntress",
            "huntress",
            (0.1, 0.2, 4),
            {"alpha": 0.1, "beta": 0.2, "n_threads": 4},
        ),
        (
            "_onconem",
            "onconem",
            "onconem",
            (0.1, 0.2),
            {"alpha": 0.1, "beta": 0.2},
        ),
        (
            "_phiscs",
            "phiscsb",
            "phiscsb",
            (0.1, 0.2),
            {"alpha": 0.1, "beta": 0.2},
        ),
        (
            "_phiscs",
            "phiscsi",
            "phiscsi",
            (0.1, 0.2, 9, 4),
            {"alpha": 0.1, "beta": 0.2, "time_limit": 9, "n_threads": 4},
        ),
        (
            "_scistree",
            "scistree",
            "scistree",
            (0.1, 0.2, 4),
            {"alpha": 0.1, "beta": 0.2, "n_threads": 4},
        ),
    ],
)
def test_simple_solver_commands_delegate(
    monkeypatch, tmp_path, module_name, command_name, solver_name, arguments, expected
):
    """Keep the thin solver commands covered without running their backends."""
    module = import_module(f"scphylo.commands.solver.{module_name}")
    command = getattr(module, command_name)
    source = tmp_path / "input.SC"
    source.touch()
    matrix = object()
    result = object()
    read = Mock(return_value=matrix)
    write = Mock()
    solve = Mock(return_value=result)
    monkeypatch.setattr(module.scp.io, "read", read)
    monkeypatch.setattr(module.scp.io, "write", write)
    monkeypatch.setattr(module.scp.tl, solver_name, solve)

    assert _invoke(command, str(source), *arguments) is None

    read.assert_called_once_with(str(source))
    solve.assert_called_once_with(matrix, **expected)
    write.assert_called_once_with(
        result, str(source.with_suffix(f".{command_name}.CFMatrix"))
    )


@pytest.mark.parametrize("alpha", [0, 0.1])
def test_gpps_command_normalizes_zero_alpha(monkeypatch, tmp_path, alpha):
    """Exercise both alpha paths in the GPPS adapter with a mocked solver."""
    module = import_module("scphylo.commands.solver._gpps")
    source = tmp_path / "input.SC"
    source.touch()
    solve = Mock(return_value="result")
    monkeypatch.setattr(module.scp.io, "read", Mock(return_value="matrix"))
    monkeypatch.setattr(module.scp.io, "write", Mock())
    monkeypatch.setattr(module.scp.tl, "gpps", solve)

    _invoke(module.gpps, str(source), alpha, 0.2, 1, 2, 3, 4, 5, 6)

    assert solve.call_args.kwargs["alpha"] == (1e-12 if alpha == 0 else alpha)


@pytest.mark.parametrize("alpha", [0, 0.1])
def test_sphyr_command_normalizes_zero_alpha(monkeypatch, tmp_path, alpha):
    """Exercise both alpha paths in the SPhyR adapter with a mocked solver."""
    module = import_module("scphylo.commands.solver._sphyr")
    source = tmp_path / "input.SC"
    source.touch()
    solve = Mock(return_value="result")
    monkeypatch.setattr(module.scp.io, "read", Mock(return_value="matrix"))
    monkeypatch.setattr(module.scp.io, "write", Mock())
    monkeypatch.setattr(module.scp.tl, "sphyr", solve)

    _invoke(module.sphyr, str(source), alpha, 0.2, 3, 4)

    assert solve.call_args.kwargs["alpha"] == (1e-12 if alpha == 0 else alpha)


@pytest.mark.parametrize("no_reconstruction", [False, True])
def test_booster_command_honors_reconstruction_flag(
    monkeypatch, tmp_path, no_reconstruction
):
    """Cover both output paths while replacing the expensive booster backend."""
    module = import_module("scphylo.commands.solver._booster")
    source = tmp_path / "input.SC"
    source.touch()
    solve = Mock(return_value="result")
    write = Mock()
    monkeypatch.setattr(module.scp.io, "read", Mock(return_value="matrix"))
    monkeypatch.setattr(module.scp.io, "write", write)
    monkeypatch.setattr(module.scp.tl, "booster", solve)

    _invoke(
        module.booster,
        str(source),
        0.1,
        0.2,
        "scite",
        "muts",
        2,
        3,
        0,
        1,
        5,
        6,
        7,
        None,
        True,
        False,
        False,
        no_reconstruction,
    )

    assert write.called is (not no_reconstruction)


@pytest.mark.parametrize("experiment", [False, True])
def test_scite_command_modes(monkeypatch, tmp_path, experiment):
    """Exercise ordinary and experiment orchestration without subprocesses."""
    module = import_module("scphylo.commands.solver._scite")
    source = tmp_path / "input.SC"
    source.touch()
    solve = Mock(return_value=("result", 1.0, 2.0, 0.3) if experiment else "result")
    monkeypatch.setattr(module.scp.io, "read", Mock(return_value="matrix"))
    monkeypatch.setattr(module.scp.io, "write", Mock())
    monkeypatch.setattr(module.scp.tl, "scite", solve)
    monkeypatch.setattr(module.scp.ul, "stat", Mock())
    monkeypatch.setattr(module, "delayed", lambda function: function)
    monkeypatch.setattr(module, "Parallel", lambda **_kwargs: lambda jobs: list(jobs))

    _invoke(module.scite, str(source), 0.1, 0.2, 10, 2, experiment, 4, 2, 5)

    assert solve.call_count == (3 if experiment else 1)


@pytest.mark.parametrize("experiment", [False, True])
def test_siclonefit_command_modes(monkeypatch, tmp_path, experiment):
    """Exercise ordinary and experiment orchestration without Java or joblib."""
    module = import_module("scphylo.commands.solver._siclonefit")
    source = tmp_path / "input.SC"
    source.touch()
    solve = Mock(return_value=("result", 1.0, True, 2.0) if experiment else "result")
    monkeypatch.setattr(module.scp.io, "read", Mock(return_value="matrix"))
    monkeypatch.setattr(module.scp.io, "write", Mock())
    monkeypatch.setattr(module.scp.tl, "siclonefit", solve)
    monkeypatch.setattr(module.scp.ul, "stat", Mock())
    monkeypatch.setattr(module, "delayed", lambda function: function)
    monkeypatch.setattr(module, "Parallel", lambda **_kwargs: lambda jobs: list(jobs))
    _invoke(
        module.siclonefit,
        str(source),
        0.1,
        0.2,
        10,
        2,
        1,
        experiment,
        4,
        2,
        5,
    )

    assert solve.call_count == (3 if experiment else 1)


def _stub_submission(monkeypatch, module):
    """Replace cluster submission and shared side effects in a caller module."""
    submit = Mock(return_value="job-id")
    monkeypatch.setattr(module, "write_cmds_get_main", submit, raising=False)
    if hasattr(module, "subprocess"):
        monkeypatch.setattr(module.subprocess, "getoutput", Mock(return_value="job-id"))
    if hasattr(module, "os"):
        monkeypatch.setattr(module.os, "system", Mock(return_value=0))
    monkeypatch.setattr(module.scp.ul, "mkdir", Mock())
    monkeypatch.setattr(module.scp.ul, "get_file", Mock(return_value="/tool.py"))
    monkeypatch.setattr(module.scp.logg, "info", Mock())
    monkeypatch.setattr(module.scp.logg, "error", Mock())
    return submit


def test_rename_command_previews_and_executes(monkeypatch, tmp_path):
    """Cover both naming schemes and preview/execute behavior."""
    module = import_module("scphylo.commands.caller._1rename")
    files = [tmp_path / f"sample_A_{mate}.fastq.gz" for mate in (1, 2)]
    for path in files:
        path.touch()
    rename = Mock()
    monkeypatch.setattr(module.os, "rename", rename)
    monkeypatch.setattr(module.scp.logg, "info", Mock())

    _invoke(module.rename, str(tmp_path), "0,1", False)
    _invoke(module.rename, str(tmp_path), ",", True)

    assert rename.call_count == 2
    assert rename.call_args_list[0].args[1].endswith("cell1_1.fastq.gz")


def test_sra_command_builds_layouts_and_checks_logs(monkeypatch, tmp_path):
    """Cover paired, single, renamed, duplicate, and log-check paths."""
    module = import_module("scphylo.commands.caller._1sra")
    _stub_submission(monkeypatch, module)
    table = tmp_path / "runs.csv"
    pd.DataFrame(
        [
            {"Library Name": "paired", "Run": "SRR1", "LibraryLayout": "PAIRED"},
            {"Library Name": "same", "Run": "same", "LibraryLayout": "SINGLE"},
            {"Library Name": "same", "Run": "SRR3", "LibraryLayout": "SINGLE"},
        ]
    ).to_csv(table, index=False)

    _invoke(module.sra, str(table), str(tmp_path), False)
    log = tmp_path / "job.o"
    log.write_text("header\n(sample)\n")
    monkeypatch.setattr(module.glob, "glob", Mock(return_value=[str(log)]))
    monkeypatch.setattr(module.subprocess, "getoutput", Mock(return_value="(sample)"))
    _invoke(module.sra, str(table), str(tmp_path), True)

    assert module.scp.logg.error.called
    assert module.scp.logg.info.called


@pytest.mark.parametrize(
    ("module_name", "command_name", "arguments", "discovery"),
    [
        ("_2fastp", "fastp", ("in", "out", "1:00", "10", None), "paired"),
        (
            "_3bwa",
            "bwa",
            ("in", "out", "hg19", "1:00", "20", "2", None, True),
            "paired",
        ),
        (
            "_5mutect2",
            "mutect2",
            ("out", "normal", "hg19", "1:00", "20", None),
            "samples",
        ),
        (
            "_5strelka2",
            "strelka2",
            ("out", "normal", "hg19", "1:00", "20", None),
            "samples",
        ),
        (
            "_5sequenza",
            "sequenza",
            ("out", "normal", "hg19", ("1", "2"), ("20", "30"), None),
            "samples",
        ),
        (
            "_5varscan2",
            "varscan2",
            ("out", "normal", "hg19", "1:00", "20", None),
            "samples",
        ),
    ],
)
def test_sample_based_caller_workflows(
    monkeypatch, module_name, command_name, arguments, discovery
):
    """Build representative jobs for callers while keeping execution mocked."""
    module = import_module(f"scphylo.commands.caller.{module_name}")
    submit = _stub_submission(monkeypatch, module)
    if discovery == "paired":
        monkeypatch.setattr(
            module.scp.ul,
            "is_paired",
            Mock(
                return_value=[
                    {"sample": "paired", "is_paired": True},
                    {"sample": "single", "is_paired": False},
                ]
            ),
        )
    else:
        monkeypatch.setattr(
            module.scp.ul,
            "get_samples_df",
            Mock(return_value=pd.DataFrame({"sample": ["tumor"]})),
        )

    _invoke(getattr(module, command_name), *arguments)

    assert submit.called
    submitted = submit.call_args_list[-1].args[0]
    assert submitted["cmd"].str.contains("Done!").all()


def test_star_builds_all_five_mocked_stages(monkeypatch):
    """Cover STAR's paired/unpaired and PDX command-generation branches."""
    module = import_module("scphylo.commands.caller._3star")
    submit = _stub_submission(monkeypatch, module)
    monkeypatch.setattr(
        module.scp.ul,
        "is_paired",
        Mock(
            return_value=[
                {"sample": "paired", "is_paired": True},
                {"sample": "single", "is_paired": False},
            ]
        ),
    )
    fake_fastq = Mock()
    fake_fastq.__enter__ = Mock(return_value=[b"@read\n", b"ACGT\n"])
    fake_fastq.__exit__ = Mock(return_value=False)
    monkeypatch.setattr(module.gzip, "open", Mock(return_value=fake_fastq))

    _invoke(
        module.star, "in", "out", "hg38", ("1", "2"), ("20", "30"), "10", None, True
    )

    assert submit.call_count == 5
    assert [call.args[1] for call in submit.call_args_list] == [
        "star-1of5",
        "star-2of5",
        "star-3of5",
        "star-4of5",
        "star-5of5",
    ]


@pytest.mark.parametrize(
    ("module_name", "command_name", "arguments", "pattern"),
    [
        (
            "_4gatk",
            "gatk",
            ("out", "hg19", "rna", "1:00", "20", None),
            "sample.mapped.bam",
        ),
        (
            "_4rsem",
            "rsem",
            ("out", "hg38", "paired", ("1", "2"), ("20", "30"), None),
            "sample.transcript.bam",
        ),
        (
            "_5mpileup",
            "mpileup",
            ("out", "m10x", "1:00", "20", None),
            "sample.markdup_bqsr.bam",
        ),
        (
            "_6defuse",
            "defuse",
            ("in", "out", "mm10", "1:00", "20", "2", None),
            "sample_R1_001.fastq",
        ),
        ("_6snpeff", "snpeff", ("out", "hg19", "1:00", "20", None), "sample.vcf"),
        ("_6varfilter", "varfilter", ("out", "hg38", "1:00", "20", None), "sample.vcf"),
        (
            "_pseudobulk",
            "pseudobulk",
            ("out", "mm10", ("1", "2"), ("20", "30"), None),
            "sample.markdup_bqsr.bam",
        ),
    ],
)
def test_file_discovery_caller_workflows(
    monkeypatch, module_name, command_name, arguments, pattern
):
    """Build jobs from a discovered input file without running shell commands."""
    module = import_module(f"scphylo.commands.caller.{module_name}")
    submit = _stub_submission(monkeypatch, module)
    monkeypatch.setattr(module.glob, "glob", Mock(return_value=[f"/data/{pattern}"]))

    _invoke(getattr(module, command_name), *arguments)

    assert submit.called


def test_gatk_and_rsem_alternate_modes(monkeypatch):
    """Cover the DNA and single-end alternatives in two pipeline adapters."""
    gatk = import_module("scphylo.commands.caller._4gatk")
    _stub_submission(monkeypatch, gatk)
    monkeypatch.setattr(
        gatk.glob, "glob", Mock(return_value=["/data/sample.mapped.bam"])
    )
    _invoke(gatk.gatk, "out", "mm10", "dna", "1:00", "20", None)

    rsem = import_module("scphylo.commands.caller._4rsem")
    _stub_submission(monkeypatch, rsem)
    monkeypatch.setattr(
        rsem.glob, "glob", Mock(return_value=["/data/sample.transcript.bam"])
    )
    _invoke(rsem.rsem, "out", "mm10", "single", ("1", "2"), ("20", "30"), None)


def test_haplotype_caller_reference_and_data_modes(monkeypatch):
    """Cover human/RNA and mouse/DNA chromosome job generation."""
    module = import_module("scphylo.commands.caller._5hapcaller")
    _stub_submission(monkeypatch, module)
    monkeypatch.setattr(
        module.glob, "glob", Mock(return_value=["/data/sample.markdup_bqsr.bam"])
    )

    _invoke(module.hapcaller, "out", "hg19", "rna", ("1", "2"), ("20", "30"), None)
    _invoke(module.hapcaller, "out", "mm10", "dna", ("1", "2"), ("20", "30"), None)


def test_bamreadcount_builds_input_and_jobs(monkeypatch, tmp_path):
    """Exercise read-count staging with a tiny real CSV and mocked jobs."""
    module = import_module("scphylo.commands.caller._6bamreadcount")
    submit = _stub_submission(monkeypatch, module)
    outdir = tmp_path / "output"
    outdir.mkdir()
    (outdir / "_tmp").mkdir()
    variants = tmp_path / "variants.csv"
    pd.DataFrame({"CHROM": ["1"], "POS": [2], "Allele": ["A"]}).to_csv(
        variants, index=False
    )
    monkeypatch.setattr(
        module.glob,
        "glob",
        Mock(return_value=[str(outdir / "sample.markdup_bqsr.bam")]),
    )

    _invoke(
        module.bamreadcount,
        str(outdir),
        str(variants),
        "hg19",
        ("1", "2"),
        ("20", "30"),
        None,
    )

    assert submit.call_count == 2
    assert (outdir / "_tmp" / "bamreadcount.input").read_text().strip() == "1\t2\t2"


def test_pyclone_reads_tumor_content(monkeypatch, tmp_path):
    """Build a PyClone job from its small Sequenza metadata file."""
    module = import_module("scphylo.commands.caller._6pyclone")
    submit = _stub_submission(monkeypatch, module)
    sample_dir = tmp_path / "sample.sequenza"
    sample_dir.mkdir()
    (sample_dir / "sample_cellularity.txt").write_text("0.75\n")
    monkeypatch.setattr(module.glob, "glob", Mock(return_value=[str(sample_dir)]))

    _invoke(module.pyclone, str(tmp_path), "1:00", "20", None)

    assert "--tumour_contents 0.75" in submit.call_args.args[0].loc[0, "cmd"]


def test_bamsnap_and_vartrix_submit_without_external_tools(monkeypatch, capsys):
    """Cover the two single-job adapters without invoking their executables."""
    bamsnap = import_module("scphylo.commands.caller._6bamsnap")
    _stub_submission(monkeypatch, bamsnap)
    _invoke(bamsnap.bamsnap, "hg19", "sample.bam", "chr1:2")
    assert "bamsnap" in capsys.readouterr().out

    vartrix = import_module("scphylo.commands.caller._6vartrix")
    submit = _stub_submission(monkeypatch, vartrix)
    monkeypatch.setattr(vartrix.scp.ul, "dir_base", Mock(return_value=("out", "x")))
    _invoke(
        vartrix.vartrix,
        "sample.bam",
        "barcodes.tsv",
        "calls.vcf",
        "m10x",
        "1:00",
        "20",
        None,
    )
    assert submit.called
