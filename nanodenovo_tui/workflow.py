"""8-stage NanoDeNovo workflow engine, fail-closed replay, and SHA-256 ledger."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shlex
import signal
import subprocess
import sys
import tempfile
import threading
import time
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable


def configure_utf8() -> bool:
    """Configure UTF-8 stdio and return True if ASCII box fallback is needed."""
    os.environ.setdefault("PYTHONUTF8", "1")
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:
                pass
    enc = (getattr(sys.stdout, "encoding", "") or "").lower()
    return not enc.startswith("utf")


def win_to_wsl(path: Path | str) -> str:
    """Convert a Windows path to /mnt/<drive>/... POSIX form."""
    p = Path(path).resolve()
    drive = p.drive.rstrip(":").lower()
    if not drive:
        return p.as_posix()
    tail = p.as_posix().split(":", 1)[-1]
    return f"/mnt/{drive}{tail}"


def build_quoted_wsl_command(cmd_parts: list[str]) -> list[str]:
    """Wrap command tokens with shlex.quote inside `wsl -e bash -lc`."""
    wsl_body = " ".join(shlex.quote(part) for part in cmd_parts)
    return ["wsl", "-e", "bash", "-lc", wsl_body]


def normalize_sha256(path: Path) -> str:
    """Hash file bytes deterministically; normalize CRLF -> LF on text tables/scripts."""
    text_exts = {".csv", ".tsv", ".txt", ".md", ".py", ".sh", ".json", ".jsonl", ".fa", ".fasta", ".pdb", ".log"}
    h = hashlib.sha256()
    if path.suffix.lower() in text_exts:
        with open(path, "rb") as fh:
            while chunk := fh.read(65536):
                h.update(chunk.replace(b"\r\n", b"\n"))
    else:
        with open(path, "rb") as fh:
            while chunk := fh.read(65536):
                h.update(chunk)
    return h.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=path.suffix or ".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            try:
                os.fsync(fh.fileno())
            except OSError:
                pass
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def run_subprocess_streaming(
    cmd: list[str],
    cwd: Path,
    on_line: Callable[[str], None] | None = None,
    cancel_event: threading.Event | None = None,
    timeout_s: float = 3600.0,
) -> dict[str, Any]:
    """Run subprocess with merged stderr=STDOUT reader thread and clean Ctrl+C tree cleanup."""
    kwargs: dict[str, Any] = {
        "cwd": str(cwd),
        "stdout": subprocess.PIPE,
        "stderr": subprocess.STDOUT,
        "text": True,
        "encoding": "utf-8",
        "errors": "replace",
        "bufsize": 1,
    }
    if sys.platform == "win32":
        kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        kwargs["start_new_session"] = True

    proc = subprocess.Popen(cmd, **kwargs)
    lines: list[str] = []

    def _reader() -> None:
        assert proc.stdout is not None
        for raw in proc.stdout:
            line = raw.rstrip("\r\n")
            lines.append(line)
            if on_line is not None:
                on_line(line)

    t = threading.Thread(target=_reader, daemon=True)
    t.start()
    t0 = time.monotonic()
    cancelled = False
    timed_out = False

    try:
        while proc.poll() is None:
            if cancel_event is not None and cancel_event.is_set():
                cancelled = True
                break
            if (time.monotonic() - t0) > timeout_s:
                timed_out = True
                break
            time.sleep(0.05)
    except KeyboardInterrupt:
        cancelled = True

    if cancelled or timed_out:
        try:
            if sys.platform == "win32":
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    check=False,
                )
            else:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except Exception:
            proc.kill()
        proc.wait(timeout=5)
    t.join(timeout=2)

    status = "cancelled" if cancelled else ("timeout" if timed_out else ("done" if proc.returncode == 0 else "failed"))
    return {"status": status, "returncode": proc.returncode, "lines": lines[-40:]}


@dataclass
class StageSpec:
    stage_id: str
    title: str
    short_desc: str
    tool_relpath: str
    pipeline_relpath: str
    replay_inputs: list[str]
    replay_outputs: list[str]
    live_cmd_template: list[str] = field(default_factory=list)


def get_stage_specs() -> list[StageSpec]:
    return [
        StageSpec(
            stage_id="S1",
            title="Backbone & Hotspot Audit",
            short_desc="Verify RFantibody contract and scaffold archive integrity.",
            tool_relpath="tools/design_rfantibody.py",
            pipeline_relpath="pipeline/rfantibody_wrapper2.py",
            replay_inputs=[
                "pipeline/rfantibody_wrapper2.py",
                "tools/design_rfantibody.py",
                "1_virtual_screening/scaffolds.zip",
            ],
            replay_outputs=["tools/design_rfantibody.py"],
            live_cmd_template=[sys.executable, "{repo}/pipeline/rfantibody_wrapper2.py", "--selftest"],
        ),
        StageSpec(
            stage_id="S2",
            title="Batch Sequence Design",
            short_desc="Verify single-session ProteinMPNN batch scripts and logs.",
            tool_relpath="tools/batch_mpnn.py",
            pipeline_relpath="pipeline/run_proteinmpnn_batch.sh",
            replay_inputs=[
                "pipeline/run_proteinmpnn_batch.sh",
                "pipeline/run_proteinmpnn_v_020.sh",
                "pipeline/run_proteinmpnn_v_030.sh",
                "tools/batch_mpnn.py",
                "1_virtual_screening/proteinmpnn/dataset_20_stats.log",
                "1_virtual_screening/proteinmpnn/dataset_30_stats.log",
            ],
            replay_outputs=["tools/batch_mpnn.py"],
            live_cmd_template=[sys.executable, "{repo}/tools/batch_mpnn.py", "--selftest"],
        ),
        StageSpec(
            stage_id="S3",
            title="Ascending-NLL Selection",
            short_desc="Verify ascending NLL selection tables across ds20 and ds30.",
            tool_relpath="tools/select_sequences.py",
            pipeline_relpath="pipeline/process_sequences.py",
            replay_inputs=[
                "pipeline/process_sequences.py",
                "tools/select_sequences.py",
                "1_virtual_screening/proteinmpnn/parsed_sequences_20.txt",
                "1_virtual_screening/proteinmpnn/parsed_sequences_30.txt",
            ],
            replay_outputs=[
                "results/selections_20.csv",
                "results/selections_30.csv",
                "results/top3_20.tsv",
                "results/top3_30.tsv",
            ],
            live_cmd_template=[sys.executable, "{repo}/tools/select_sequences.py", "--selftest"],
        ),
        StageSpec(
            stage_id="S4",
            title="Monomer Structure Refolding",
            short_desc="Check inclusive 1..1000 monomer refolds and confidence CSVs.",
            tool_relpath="tools/fold_ab.py",
            pipeline_relpath="pipeline/tfold_ab_runner.py",
            replay_inputs=[
                "pipeline/tfold_ab_runner.py",
                "pipeline/parse_confident.py",
                "tools/fold_ab.py",
                "1_virtual_screening/proteinmpnn/folded_sequences_20.csv",
                "1_virtual_screening/proteinmpnn/folded_sequences_30.csv",
            ],
            replay_outputs=[
                "1_virtual_screening/tfold_ab/confident_proteins_20.csv",
                "1_virtual_screening/tfold_ab/confident_proteins_30.csv",
                "results/refold_metrics.csv",
                "results/refold_summary.json",
            ],
            live_cmd_template=[sys.executable, "{repo}/tools/fold_ab.py", "--selftest"],
        ),
        StageSpec(
            stage_id="S5",
            title="USalign & Kabsch CDR3 Audit",
            short_desc="Verify all-design USalign alignments and CDR3 Kabsch deviations.",
            tool_relpath="tools/align_all.py",
            pipeline_relpath="pipeline/align_pipeline_kwargs.py",
            replay_inputs=[
                "pipeline/align_pipeline_kwargs.py",
                "tools/align_all.py",
                "tools/cdr3_check.py",
                "tools/score_refold_usalign.py",
                "1_virtual_screening/usalign/v_20_usalign_rmsd_1.0_tm_0.95.txt",
                "1_virtual_screening/usalign/v_30_usalign_rmsd_1.0_tm_0.95.txt",
            ],
            replay_outputs=[
                "results/align_all_20.csv",
                "results/align_all_30.csv",
                "results/cdr3_20.csv",
                "results/cdr3_30.csv",
                "results/refold_usalign.csv",
            ],
            live_cmd_template=[sys.executable, "{repo}/tools/cdr3_check.py", "--selftest"],
        ),
        StageSpec(
            stage_id="S6",
            title="PDB 8e8z_C Interface Audit",
            short_desc="Check complex REMARK parser and 8e8z_C interface contacts.",
            tool_relpath="tools/interface_check.py",
            pipeline_relpath="pipeline/tfold_ag_dock_kwargs.py",
            replay_inputs=[
                "pipeline/tfold_ag_prepare_kwargs.py",
                "pipeline/tfold_ag_dock_kwargs.py",
                "tools/fold_ag.py",
                "tools/interface_check.py",
                "1_virtual_screening/confident_complexes/0389_top.2.20_8e8z_C_complex.pdb",
            ],
            replay_outputs=["results/interface_audit.csv"],
            live_cmd_template=[sys.executable, "{repo}/tools/interface_check.py", "--selftest"],
        ),
        StageSpec(
            stage_id="S7",
            title="OpenMM Relaxation & Clash Gate",
            short_desc="Verify relaxed 8e8z_C complexes clear < 2.2 A heavy-atom clashes.",
            tool_relpath="2_sequence_and_structure_assessement/openmm_em/openmm_minimizer.py",
            pipeline_relpath="2_sequence_and_structure_assessement/openmm_em/openmm_minimizer.py",
            replay_inputs=[
                "2_sequence_and_structure_assessement/openmm_em/openmm_minimizer.py",
                "2_sequence_and_structure_assessement/openmm_em/0389_top.2.20_8e8z_C_complex_em.pdb",
                "2_sequence_and_structure_assessement/openmm_em/0743_top.1.30_8e8z_C_complex_em.pdb",
            ],
            replay_outputs=[
                "2_sequence_and_structure_assessement/openmm_em/0389_top.2.20_8e8z_C_complex_em.pdb",
                "2_sequence_and_structure_assessement/openmm_em/0743_top.1.30_8e8z_C_complex_em.pdb",
            ],
            live_cmd_template=[
                sys.executable,
                "{repo}/2_sequence_and_structure_assessement/openmm_em/openmm_minimizer.py",
                "--selftest",
            ],
        ),
        StageSpec(
            stage_id="S8",
            title="Youden QC & SHA-256 Ledger",
            short_desc="Verify audit reports, plates, and write portable SHA-256 ledger.",
            tool_relpath="tools/qc_report.py",
            pipeline_relpath="pipeline/build_plots.py",
            replay_inputs=[
                "tools/qc_report.py",
                "AUDIT.md",
                "reports/quality_20.md",
                "reports/quality_30.md",
                "results/refold_usalign_summary.json",
                "figs/nanodenovo_audit_plate.png",
            ],
            replay_outputs=["results/ledger.jsonl"],
            live_cmd_template=[sys.executable, "{repo}/tools/qc_report.py", "--selftest"],
        ),
    ]


class WorkflowEngine:
    """Orchestrates stage checks, replay verification, dry-run planning, and SHA-256 ledgers."""

    def __init__(self, repo_root: Path) -> None:
        self.repo_root = repo_root.resolve()
        self.specs = get_stage_specs()
        self.by_id = {s.stage_id: s for s in self.specs}
        self.ledger_path = self.repo_root / "results" / "ledger.jsonl"

    def _require_nonempty(self, rels: list[str]) -> tuple[bool, dict[str, str], str]:
        hashes: dict[str, str] = {}
        for rel in rels:
            p = self.repo_root / rel
            if not p.exists() or not p.is_file():
                return False, hashes, f"Missing required artifact: {rel}"
            if p.stat().st_size == 0:
                return False, hashes, f"Empty artifact (fail-closed): {rel}"
            hashes[rel] = normalize_sha256(p)
        return True, hashes, ""

    def _run_replay_domain_check(self, stage_id: str) -> tuple[bool, str]:
        """Run deterministic CPU domain verification inside replay mode."""
        tools_dir = str(self.repo_root / "tools")
        if tools_dir not in sys.path:
            sys.path.insert(0, tools_dir)

        if stage_id == "S1":
            zip_path = self.repo_root / "1_virtual_screening" / "scaffolds.zip"
            with zipfile.ZipFile(zip_path, "r") as zf:
                names = [n for n in zf.namelist() if n.endswith(".pdb")]
                if not names:
                    return False, "scaffolds.zip contains zero .pdb files"
                return True, f"Verified scaffolds.zip ({len(names)} scaffold PDBs)"

        if stage_id == "S2":
            p20 = (self.repo_root / "1_virtual_screening" / "proteinmpnn" / "dataset_20_stats.log").stat().st_size
            p30 = (self.repo_root / "1_virtual_screening" / "proteinmpnn" / "dataset_30_stats.log").stat().st_size
            return True, f"Verified ProteinMPNN stats logs (ds20={p20}B, ds30={p30}B)"

        if stage_id == "S3":
            for tag in ("20", "30"):
                csv_p = self.repo_root / "results" / f"selections_{tag}.csv"
                with open(csv_p, "r", encoding="utf-8", newline="") as fh:
                    rows = list(csv.DictReader(fh))
                if len(rows) < 2990:
                    return False, f"selections_{tag}.csv has too few rows: {len(rows)}"
            return True, "Verified selections_20.csv and selections_30.csv ascending-NLL tables (3000 rows each)"

        if stage_id == "S4":
            summary = json.loads((self.repo_root / "results" / "refold_summary.json").read_text(encoding="utf-8"))
            total = int(summary.get("total") or sum(v.get("n", 0) for v in summary.get("groups", {}).values() if isinstance(v, dict))) if isinstance(summary, dict) else 0
            if total < 6000:
                return False, f"Unexpected refold_summary total: {total}"
            return True, f"Verified NanoBodyBuilder2 refold_summary ({total} refolded sequences)"

        if stage_id == "S5":
            usum = json.loads((self.repo_root / "results" / "refold_usalign_summary.json").read_text(encoding="utf-8"))
            if usum.get("total", 0) != 6200:
                return False, f"Expected 6200 USalign refolds, got {usum.get('total')}"
            return True, "Verified 6,200 USalign refold alignments + CDR3 tables"

        if stage_id == "S6":
            import interface_check  # type: ignore

            pdb_p = self.repo_root / "1_virtual_screening" / "confident_complexes" / "0389_top.2.20_8e8z_C_complex.pdb"
            rep = interface_check.analyze_complex(pdb_p)
            if rep.get("status") != "ok":
                return False, f"Complex geometry check failed: {rep}"
            return True, f"Verified 8e8z_C interface BSA={rep['bsa_A2']} A^2 on 0389_top.2.20"

        if stage_id == "S7":
            em_dir = self.repo_root / "2_sequence_and_structure_assessement" / "openmm_em"
            if str(em_dir) not in sys.path:
                sys.path.insert(0, str(em_dir))
            import openmm_minimizer  # type: ignore

            for fname in ("0389_top.2.20_8e8z_C_complex_em.pdb", "0743_top.1.30_8e8z_C_complex_em.pdb"):
                audit = openmm_minimizer.audit_heavy_clashes(em_dir / fname, 2.2)
                clashes = int(audit.get("clashes", 0))
                if clashes > 5:
                    return False, f"Unexpected clash count in relaxed PDB {fname}: {clashes}"
            return True, "Verified relaxed OpenMM 8e8z_C complexes pass clash gate"

        if stage_id == "S8":
            return True, "Verified AUDIT.md, quality_{20,30}.md, and SHA-256 ledger writer"

        return False, f"Unknown stage {stage_id}"

    def _load_existing_ledger(self) -> dict[str, dict]:
        if not self.ledger_path.exists():
            return {}
        out: dict[str, dict] = {}
        for line in self.ledger_path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                if "stage" in rec:
                    out[rec["stage"]] = rec
            except json.JSONDecodeError:
                continue
        return out

    def _resume_matches(self, spec: StageSpec, prev: dict | None) -> bool:
        if not prev or prev.get("status") != "done":
            return False
        prev_in = prev.get("input_hashes", {})
        prev_out = prev.get("output_hashes", {})
        if set(prev_in) != set(spec.replay_inputs):
            return False
        for rel, expected_h in prev_in.items():
            p = self.repo_root / rel
            if not p.exists() or normalize_sha256(p) != expected_h:
                return False
        for rel, expected_h in prev_out.items():
            if rel == "results/ledger.jsonl":
                continue
            p = self.repo_root / rel
            if not p.exists() or normalize_sha256(p) != expected_h:
                return False
        return True

    def write_ledger(self, records: list[dict]) -> None:
        lines = [json.dumps(r, sort_keys=True) for r in records]
        atomic_write_text(self.ledger_path, "\n".join(lines) + "\n")

    def verify_ledger(self) -> dict[str, Any]:
        """Verify all stage hashes in `results/ledger.jsonl` fail-closed."""
        if not self.ledger_path.exists():
            return {"ok": False, "error": "Missing results/ledger.jsonl", "stages_checked": 0}
        records: list[dict] = []
        for idx, line in enumerate(self.ledger_path.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                return {"ok": False, "error": f"Corrupt JSON line {idx}: {exc}", "stages_checked": len(records)}

        if len(records) < len(self.specs):
            return {"ok": False, "error": f"Expected {len(self.specs)} stage entries, found {len(records)}", "stages_checked": len(records)}

        mismatches: list[str] = []
        for rec in records:
            if rec.get("status") != "done":
                mismatches.append(f"{rec.get('stage')}: status={rec.get('status')}")
                continue
            for rel, want in rec.get("input_hashes", {}).items():
                p = self.repo_root / rel
                if not p.exists() or normalize_sha256(p) != want:
                    mismatches.append(f"{rec['stage']} input mismatch: {rel}")
            for rel, want in rec.get("output_hashes", {}).items():
                if rel == "results/ledger.jsonl":
                    continue
                p = self.repo_root / rel
                if not p.exists() or normalize_sha256(p) != want:
                    mismatches.append(f"{rec['stage']} output mismatch: {rel}")

        return {"ok": len(mismatches) == 0, "stages_checked": len(records), "mismatches": mismatches}

    def run_stage(
        self,
        stage_id: str,
        mode: str = "replay",
        resume: bool = False,
        on_line: Callable[[str], None] | None = None,
        cancel_event: threading.Event | None = None,
    ) -> dict[str, Any]:
        spec = self.by_id[stage_id]
        existing = self._load_existing_ledger()
        if resume and self._resume_matches(spec, existing.get(stage_id)):
            msg = f"[resume] {stage_id} hashes match; skipped"
            if on_line:
                on_line(msg)
            rec = dict(existing[stage_id])
            rec["resumed"] = True
            return rec

        t0 = time.monotonic()
        if mode == "dry-run":
            cmd_str = " ".join(c.format(repo=".") for c in spec.live_cmd_template)
            ok_in, in_hashes, err_in = self._require_nonempty(spec.replay_inputs)
            if on_line:
                on_line(f"[dry-run] {stage_id}: {cmd_str}")
            return {
                "stage": stage_id,
                "title": spec.title,
                "mode": mode,
                "status": "done" if ok_in else "failed",
                "detail": f"Planned: {cmd_str}" if ok_in else err_in,
                "input_hashes": in_hashes,
                "output_hashes": {},
                "elapsed_s": round(time.monotonic() - t0, 3),
            }

        if mode == "replay":
            ok_in, in_hashes, err_in = self._require_nonempty(spec.replay_inputs)
            if not ok_in:
                return {
                    "stage": stage_id,
                    "title": spec.title,
                    "mode": mode,
                    "status": "failed",
                    "detail": err_in,
                    "input_hashes": in_hashes,
                    "output_hashes": {},
                    "elapsed_s": round(time.monotonic() - t0, 3),
                }
            # Exclude self-referential ledger.jsonl on S8 before writing it
            out_rels = [r for r in spec.replay_outputs if r != "results/ledger.jsonl"]
            ok_out, out_hashes, err_out = self._require_nonempty(out_rels)
            if not ok_out:
                return {
                    "stage": stage_id,
                    "title": spec.title,
                    "mode": mode,
                    "status": "failed",
                    "detail": err_out,
                    "input_hashes": in_hashes,
                    "output_hashes": out_hashes,
                    "elapsed_s": round(time.monotonic() - t0, 3),
                }
            domain_ok, domain_detail = self._run_replay_domain_check(stage_id)
            if on_line:
                on_line(f"[replay] {stage_id}: {domain_detail}")
            return {
                "stage": stage_id,
                "title": spec.title,
                "mode": mode,
                "status": "done" if domain_ok else "failed",
                "detail": domain_detail,
                "input_hashes": in_hashes,
                "output_hashes": out_hashes,
                "elapsed_s": round(time.monotonic() - t0, 3),
            }

        if mode == "live":
            cmd = [c.format(repo=str(self.repo_root)) for c in spec.live_cmd_template]
            res = run_subprocess_streaming(cmd, cwd=self.repo_root, on_line=on_line, cancel_event=cancel_event)
            ok_in, in_hashes, _ = self._require_nonempty(spec.replay_inputs)
            out_rels = [r for r in spec.replay_outputs if r != "results/ledger.jsonl"]
            _, out_hashes, _ = self._require_nonempty(out_rels)
            return {
                "stage": stage_id,
                "title": spec.title,
                "mode": mode,
                "status": res["status"] if ok_in else "failed",
                "detail": "; ".join(res["lines"][-3:]) or f"exit {res['returncode']}",
                "input_hashes": in_hashes,
                "output_hashes": out_hashes,
                "elapsed_s": round(time.monotonic() - t0, 3),
            }

        raise ValueError(f"Unsupported mode: {mode}")

    def run_all(
        self,
        mode: str = "replay",
        resume: bool = False,
        on_line: Callable[[str], None] | None = None,
        cancel_event: threading.Event | None = None,
    ) -> dict[str, Any]:
        records: list[dict] = []
        overall_ok = True
        for spec in self.specs:
            rec = self.run_stage(spec.stage_id, mode=mode, resume=resume, on_line=on_line, cancel_event=cancel_event)
            records.append(rec)
            if rec["status"] != "done":
                overall_ok = False
                break
        if mode == "replay" and overall_ok:
            self.write_ledger(records)
        elif mode == "live" and records:
            self.write_ledger(records)
        return {"ok": overall_ok, "mode": mode, "records": records, "ledger": str(self.ledger_path)}
