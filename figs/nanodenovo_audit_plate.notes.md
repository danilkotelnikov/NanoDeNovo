<!-- version: 0.3.0-a43ea1-4 | date: 15.XI.2026 -->
# NanoDeNovo Audit Plate (`nanodenovo_audit_plate`)

## Figure Summary

This plate audits VHH designs against PDB entry 8e8z chain C (`PDB 8e8z_C`). Panels cover ProteinMPNN ranking, USalign refolding (`n=6200`), pTM cutoffs, and interface clashes. Sorting lower negative log-likelihood (`NLL`) recovers `5.6%` placed CDR3 loops (`<= 3.0` Å) in both datasets. Zero refolded designs meet both `TM >= 0.9` and `RMSD < 1.5` Å, while five relaxed complexes eliminate all `< 2.2` Å clashes.

## Specification and Automated Checks

- **Canvas and fonts**: The canvas measures `178.1 mm` by `122.1 mm`. All labels use Arial (`ArialMT`, `Arial-BoldMT`) at `6.55 pt` to `8.00 pt`, inside a `0.45 mm` `#6E7781` frame.
- **Palette**: Colors use `#1D397B` (repaired / positive metric), `#317F45` (relaxed / false positives), `#9A5600` (baseline bug / clashes / missed), and `#6E7781`. Minimum colorblind delta is `CIEDE2000 = 10.7`.
- **Molecullage QA (`15.XI.2026`)**: Normal plate QA passed (`ok: true`, `ink: 7.42%`, `cells: 8.75%`, `band: 0%`). Blank share is `2.93%`, type floor is `6.55 pt`, word count is `161`, palette score is `10.7`, overplate share is `0%`, and box share is `0%`. Open-paper plots omit box fill, keeping ink at `7.42%` under the project normal plate gate.
- **GitHub edge contrast (`ghfit.py --no-log`)**: The PNG passed all nine light, dark, and dimmed preview checks (`worst 3.69:1`).
- **Blind read (`15.XI.2026`, Round 3)**: An independent `research` subagent reviewed the updated scratch PNG against `qa/blind_read.md` and returned `PASS`. It summarized that the sort fix lowers `Mean sequence NLL` below `Orig bug mean` (`1.026` and `1.079`) and reaches `5.6%` CDR3 recovery. Raising `pTM cutoff` from `0.80` to `0.85` drops `Youden J` from `0.698` to `0.540`, and all five relaxed models sit at `2.42–2.66` Å with zero clashes.

## Claim Boundary

All metrics represent static sequence scoring, NanoBodyBuilder2 predictions aligned with USalign against designed VHH scaffold backbones for `PDB 8e8z_C`, and local energy minimization audits. Wet-lab binding and developability remain untested.

## Data Provenance

- **Panel a**: Paired `NLL` distributions come from `AUDIT.md` and selection tables across `ds20` (`n=3000`) and `ds30` (`n=3000`). The plot compares repaired top-3 selections, worst selections, and original bug means (`1.026` and `1.079`).
- **Panel b**: USalign TM-score summaries and placed CDR3 fractions come from `results/refold_usalign_summary.json` and `refold_summary.json`. Counts cover `corr_20` (`145/2608 = 5.6%`), `orig_20` (`2/95 = 2.1%`), `corr_30` (`148/2665 = 5.6%`), and `orig_30` (`4/93 = 4.3%`).
- **Panel c**: Cutoff counts across `pTM` `0.70`, `0.75`, `0.80`, and `0.85` come from `reports/quality_20.md`. Evaluation uses `ds20` designs (`pLDDT > 0.90`, `n=2997`) against refold truth (`TM >= 0.90` and `RMSD <= 2.0` Å).
- **Panel d**: Heavy-atom clash counts (`< 2.2` Å) and minimum distances across `n=11` models come from `results/interface_audit.csv`. Unrelaxed Rosetta models show `18–59` clashes at `0.54–1.07` Å, while five relaxed models show `0` clashes at `2.42–2.66` Å.

## References

- Kotelnikov, S. et al. Structure-based analysis involving PDB entry `8e8z` chain `C` (`DOI: 10.3390/ijms26199262`, CC BY 4.0).
- Watson, J. L. et al. De novo design of protein structure and function with RFdiffusion (`DOI: 10.1038/s41586-023-06415-8`).
- Dauparas, J. et al. Robust deep learning-based protein sequence design using ProteinMPNN (`DOI: 10.1126/science.add2187`).
- Zhang, C. et al. US-align: universal structure alignment of monomeric and complex biomolecules (`DOI: 10.1038/s41592-022-01585-1`).
- Eastman, P. et al. OpenMM 8: molecular dynamics simulation with machine learning potentials (`DOI: 10.1021/acs.jpcb.3c06662`).
