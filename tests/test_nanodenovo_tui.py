"""Automated tests for NanoDeNovo TUI, workflow engine, CLI flags, and pipeline edge cases."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "pipeline") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pipeline"))

from nanodenovo_tui.cli import main as cli_main  # noqa: E402
from nanodenovo_tui.deps import run_doctor  # noqa: E402
from nanodenovo_tui.rich_tui import render_dashboard_text  # noqa: E402
from nanodenovo_tui.workflow import WorkflowEngine, build_quoted_wsl_command, win_to_wsl  # noqa: E402
import process_sequences  # noqa: E402
import tfold_ag_dock_kwargs  # noqa: E402


def _make_mini_repo(tmp_path: Path) -> Path:
    """Create a minimal copy of required replay artifacts inside tmp_path without mutating the live repo."""
    copy_dirs = ("config", "docs", "pipeline", "tools", "reports", "figs", "results")
    for d in copy_dirs:
        src = REPO_ROOT / d
        if src.exists():
            shutil.copytree(src, tmp_path / d)
    # Copy lightweight specific artifacts referenced by S1..S7 replay_inputs/outputs
    rel_files = [
        "AUDIT.md",
        "1_virtual_screening/scaffolds.zip",
        "1_virtual_screening/proteinmpnn/dataset_20_stats.log",
        "1_virtual_screening/proteinmpnn/dataset_30_stats.log",
        "1_virtual_screening/proteinmpnn/parsed_sequences_20.txt",
        "1_virtual_screening/proteinmpnn/parsed_sequences_30.txt",
        "1_virtual_screening/proteinmpnn/folded_sequences_20.csv",
        "1_virtual_screening/proteinmpnn/folded_sequences_30.csv",
        "1_virtual_screening/tfold_ab/confident_proteins_20.csv",
        "1_virtual_screening/tfold_ab/confident_proteins_30.csv",
        "1_virtual_screening/usalign/v_20_usalign_rmsd_1.0_tm_0.95.txt",
        "1_virtual_screening/usalign/v_30_usalign_rmsd_1.0_tm_0.95.txt",
        "1_virtual_screening/confident_complexes/0389_top.2.20_8e8z_C_complex.pdb",
        "2_sequence_and_structure_assessement/openmm_em/openmm_minimizer.py",
        "2_sequence_and_structure_assessement/openmm_em/0389_top.2.20_8e8z_C_complex_em.pdb",
        "2_sequence_and_structure_assessement/openmm_em/0743_top.1.30_8e8z_C_complex_em.pdb",
    ]
    for rel in rel_files:
        src = REPO_ROOT / rel
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    # Copy any top-level dirs expected by run_doctor
    (tmp_path / "nanodenovo_tui").mkdir(parents=True, exist_ok=True)
    return tmp_path


def test_dry_run_and_replay_all_stages(tmp_path: Path) -> None:
    mini = _make_mini_repo(tmp_path)
    engine = WorkflowEngine(mini)

    dry = engine.run_all(mode="dry-run")
    assert dry["ok"] is True
    assert len(dry["records"]) == 8
    assert [r["stage"] for r in dry["records"]] == [f"S{i}" for i in range(1, 9)]

    rep = engine.run_all(mode="replay")
    assert rep["ok"] is True
    assert len(rep["records"]) == 8
    assert engine.ledger_path.exists()

    v = engine.verify_ledger()
    assert v["ok"] is True
    assert v["stages_checked"] == 8

    # Resume should skip matching S1..S8 stages
    resumed = engine.run_all(mode="replay", resume=True)
    assert resumed["ok"] is True
    assert all(r.get("resumed") is True for r in resumed["records"])


def test_fail_closed_missing_and_tampered_artifact(tmp_path: Path) -> None:
    mini = _make_mini_repo(tmp_path)
    engine = WorkflowEngine(mini)
    assert engine.run_all(mode="replay")["ok"] is True

    # Tamper with an output referenced in the ledger and ensure verify_ledger fails closed
    target = mini / "results" / "top3_20.tsv"
    target.write_text("tampered\n", encoding="utf-8")
    v = engine.verify_ledger()
    assert v["ok"] is False
    assert any("top3_20.tsv" in m for m in v["mismatches"])

    # Remove required input and verify run_stage fails closed
    missing_in = mini / "1_virtual_screening" / "scaffolds.zip"
    missing_in.unlink()
    s1 = engine.run_stage("S1", mode="replay")
    assert s1["status"] == "failed"


def test_cli_doctor_instructions_verify_and_selftest(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    mini = _make_mini_repo(tmp_path)
    assert cli_main(["--repo", str(mini), "--doctor", "--json"]) == 0
    doc = json.loads(capsys.readouterr().out)
    assert doc["ok"] is True

    assert cli_main(["--repo", str(mini), "--instructions", "--json"]) == 0
    ins = json.loads(capsys.readouterr().out)
    assert len(ins["sections"]) >= 2

    assert cli_main(["--repo", str(mini), "--headless", "--mode", "replay", "--json"]) == 0
    head = json.loads(capsys.readouterr().out)
    assert head["ok"] is True

    assert cli_main(["--repo", str(mini), "--verify-ledger", "--json"]) == 0
    ver = json.loads(capsys.readouterr().out)
    assert ver["ok"] is True

    assert cli_main(["--repo", str(mini), "--selftest"]) == 0


def test_nontty_rich_tui_render_and_cli_fallback(tmp_path: Path) -> None:
    mini = _make_mini_repo(tmp_path)
    engine = WorkflowEngine(mini)
    frame = render_dashboard_text(engine, mode="dry-run", ascii_only=True)
    for idx in range(1, 9):
        assert f"S{idx}" in frame
    assert "NanoDeNovo Workflow Console" in frame

    # Subprocess with closed stdin confirms non-TTY `--ui rich` exits 0 without blocking
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "nanodenovo.py"), "--repo", str(mini), "--ui", "rich"],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=15,
        check=False,
    )
    assert proc.returncode == 0
    assert "S1" in proc.stdout and "S8" in proc.stdout


def test_wsl_quoting_with_spaces_and_cyrillic() -> None:
    path = Path(r"C:\work\Desktop Рабочий\NanoDeNovo\tool check.py")
    wsl_p = win_to_wsl(path) if sys.platform == "win32" else str(path)
    cmd = build_quoted_wsl_command(["python3", wsl_p, "--flag"])
    assert cmd[:4] == ["wsl", "-e", "bash", "-lc"]
    assert "'" in cmd[4]
    assert "tool check.py" in cmd[4]


def test_pipeline_multiline_fasta_and_ascending_nll_tiebreak(tmp_path: Path) -> None:
    fasta = tmp_path / "sampled.fa"
    fasta.write_text(
        ">native, score=0.500\n"
        "AAAA\n"
        ">candA, sample=2, score=1.250\n"
        "EVQLV\n"
        "ESGGG\n"
        ">candB, sample=1, score=1.250\n"
        "EVQLVESGGX\n"
        ">candC, sample=3, score=1.100\n"
        "QVQLQ\n",
        encoding="utf-8",
    )
    recs = process_sequences.read_sampled_records(fasta)
    # Native record without sample= must be excluded; multiline sequence must join cleanly
    assert len(recs) == 3
    assert (1.250, 2, "EVQLVESGGG") in recs

    top = process_sequences.select_top_unique(recs, top=3)
    # Lowest NLL first (1.100), then deterministic sample_id ascending tie-break (sample=1 before sample=2)
    assert top[0] == (1, 1.100, 3, "QVQLQ")
    assert top[1] == (2, 1.250, 1, "EVQLVESGGX")
    assert top[2] == (3, 1.250, 2, "EVQLVESGGG")


def test_pipeline_tfold_ag_dock_antigen_line_first(tmp_path: Path) -> None:
    pdb = tmp_path / "complex.pdb"
    # Deliberately place Antigen Predicted lDDT-Ca score BEFORE Predicted lDDT-Ca score
    pdb.write_text(
        "REMARK 250 Antigen Predicted lDDT-Ca score: 91.25\n"
        "REMARK 250 Predicted ipTM score: 0.78\n"
        "REMARK 250 Predicted pTM score: 0.84\n"
        "REMARK 250 Predicted lDDT-Ca score: 88.50\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 88.50           C\n",
        encoding="utf-8",
    )
    stats = tfold_ag_dock_kwargs.parse_pdb_remarks(pdb)
    assert stats is not None
    assert stats["ag_lddt_ca"] == pytest.approx(91.25)
    assert stats["lddt_ca"] == pytest.approx(88.50)
    assert stats["iptm"] == pytest.approx(0.78)
    assert stats["ptm"] == pytest.approx(0.84)
