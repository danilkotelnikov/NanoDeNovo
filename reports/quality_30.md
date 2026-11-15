<!-- version: w06 | date: 14.XI.2026 -->
# Quality report, dataset 30

## What this report measures

The screen under test keeps a design when pLDDT is above 0.9 and pTM is above 0.85.
The reference standard calls a design correct when its refolded model matches its backbone. The threshold is TM-score at or above 0.9 and RMSD at or below 2.0 angstrom.

The reference is a calculation, not a laboratory measurement. So the share below is a self-consistency failure rate among the designs the filter kept. It is a false discovery rate. It is not a false positive rate, because that would need the rejected designs as well.

## Integrity

- Folded designs in the table: 2997
- Repeated rows: 0
- 1 backbone(s) absent up to 1000; absent: 1000
- 3000 designs expected, 2997 present (3 absent)

## Confidence filter

- Designs that pass: 1417
- Designs that fail: 1580
- pLDDT, all designs: n=2997, mean=0.8927, median=0.9205, sd=0.1042, min=0.3005, max=0.9808
- pTM, all designs: n=2997, mean=0.8087, median=0.8457, sd=0.1295, min=0.1422, max=0.9173
- pLDDT, designs that pass: n=1417, mean=0.9528, median=0.9566, sd=0.0172, min=0.9023, max=0.9808
- pTM, designs that pass: n=1417, mean=0.8846, median=0.8869, sd=0.0178, min=0.8501, max=0.9173

## Designs that pass the filter but fail the structure check

- Designs that pass the filter and also have a structural measurement: 1417
- Correct by the reference standard: 1378
- Failures: 39
- Failure share of the measured set: 0.028
- TM-score of the measured set: n=1417, mean=0.9591, median=0.9659, sd=0.0216, min=0.8550, max=0.9896
- RMSD of the measured set: n=1417, mean=0.9558, median=0.9000, sd=0.3082, min=0.4000, max=2.3100

### How much the structural bar decides

A nanobody has a long and flexible third loop. A strict bar fails a design whose framework is right and whose loop moved, so the bar is shown at four settings.

| Bar | TM-score at or above | RMSD at or below | Failures | Share |
|---|---|---|---|---|
| strict | 0.95 | 1.0 | 495 | 0.349 |
| strict match, relaxed RMSD | 0.95 | 2.0 | 344 | 0.243 |
| common design bar | 0.9 | 2.0 | 39 | 0.028 |
| loose | 0.85 | 3.0 | 0 | 0.000 |

The ten weakest failures at the chosen bar:

| Design | pLDDT | pTM | TM-score | RMSD |
|---|---|---|---|---|
| 0551_top.3.30 | 0.9173 | 0.8629 | 0.8550 | 1.66 |
| 0378_top.3.30 | 0.9048 | 0.8536 | 0.8586 | 2.06 |
| 0070_top.1.30 | 0.9031 | 0.8520 | 0.8627 | 2.07 |
| 0002_top.1.30 | 0.9023 | 0.8576 | 0.8635 | 1.30 |
| 0076_top.3.30 | 0.9078 | 0.8548 | 0.8679 | 1.62 |
| 0277_top.1.30 | 0.9497 | 0.8919 | 0.8692 | 1.95 |
| 0224_top.2.30 | 0.9148 | 0.8538 | 0.8696 | 1.81 |
| 0447_top.3.30 | 0.9053 | 0.8587 | 0.8707 | 1.91 |
| 0300_top.2.30 | 0.9253 | 0.8522 | 0.8741 | 2.17 |
| 0378_top.1.30 | 0.9164 | 0.8561 | 0.8781 | 2.05 |

## False negatives

Every one of the 1580 rejected designs has a structural measurement, so the filter can be scored in both directions. A miss is a design that the reference calls correct and the filter rejected.

| Bar | Pass | Fail | Missed | Rej bad | FDR | Miss rate | FPR | Accuracy |
|---|---|---|---|---|---|---|---|---|
| strict (0.95/1.0) | 922 | 495 | 184 | 1396 | 0.349 | 0.166 | 0.262 | 0.369 |
| strict match, relaxed RMSD (0.95/2.0) | 1073 | 344 | 260 | 1320 | 0.243 | 0.195 | 0.207 | 0.445 |
| common design bar (0.9/2.0) | 1378 | 39 | 921 | 659 | 0.028 | 0.401 | 0.056 | 0.767 |
| loose (0.85/3.0) | 1417 | 0 | 1199 | 381 | 0.000 | 0.458 | 0.000 | 0.873 |

Failures among kept is the false discovery rate. Miss rate is the share of correct designs that the filter threw away. Wrong designs kept is the false positive rate. All three are measured against structural self-consistency, which is a calculation and not a laboratory result.

Separation at the chosen bar: a correct design scores higher pLDDT than a wrong one with probability 0.938, and higher pTM with probability 0.913. A value of 0.5 means the score tells nothing; 1.0 means it separates them perfectly. The whole chain is dominated by the conserved framework, so this value is generous. Read the third-loop check for the harder test.

## Threshold choice

Each row shows how many measured designs a threshold pair keeps, and how many of those are false positives.
Every threshold pair is scored on all designs. Rows are ordered by Youden's J, the share of correct designs kept minus the share of wrong designs kept. J is 0 for a filter that keeps everything, however common correct designs are.

| pLDDT above | pTM above | Kept | False positives | Missed | Failure share | Miss rate | J |
|---|---|---|---|---|---|---|---|
| 0.88 | 0.7 | 2230 | 149 | 218 | 0.067 | 0.095 | 0.692 |
| 0.88 | 0.75 | 2230 | 149 | 218 | 0.067 | 0.095 | 0.692 |
| 0.9 | 0.7 | 1921 | 83 | 461 | 0.043 | 0.201 | 0.681 |
| 0.9 | 0.75 | 1921 | 83 | 461 | 0.043 | 0.201 | 0.681 |
| 0.9 | 0.8 | 1914 | 83 | 468 | 0.043 | 0.204 | 0.678 |
| 0.88 | 0.8 | 2127 | 134 | 306 | 0.063 | 0.133 | 0.675 |
| 0.85 | 0.8 | 2166 | 156 | 289 | 0.072 | 0.126 | 0.651 |
| 0.8 | 0.8 | 2167 | 157 | 289 | 0.072 | 0.126 | 0.649 |
| 0.92 | 0.7 | 1508 | 24 | 815 | 0.016 | 0.355 | 0.611 |
| 0.92 | 0.75 | 1508 | 24 | 815 | 0.016 | 0.355 | 0.611 |
| 0.92 | 0.8 | 1508 | 24 | 815 | 0.016 | 0.355 | 0.611 |
| 0.85 | 0.75 | 2461 | 248 | 86 | 0.101 | 0.037 | 0.607 |

The study's own cut (pLDDT above 0.9, pTM above 0.85) scores J = 0.544. It keeps 1417 designs, keeps 39 wrong ones and misses 921.

The best row sits on the edge of the grid. A wider grid could do better, so this table shows a direction and not an optimum.

J ignores what the later stages cost. A looser filter sends more designs to docking and simulation. Pick a row with that cost in mind.

## Third-loop check

A whole-chain match is dominated by the conserved framework, so most designs pass it. This check fits the framework alone and measures how far the third loop (CDR3) moved. The loop is found by a sequence heuristic, so some designs cannot be placed.

- Designs placed: 2694 of 2997
- Kept designs placed: 1400 of 1417
- Rejected designs placed: 1294 of 1580

| Measure | Kept | Rejected |
|---|---|---|
| Median framework RMSD (angstrom) | 0.77 | 0.93 |
| Median CDR3 RMSD (angstrom) | 1.15 | 3.49 |
| Median CDR3 length (residues) | 9 | 12 |

Loop length matters. If kept designs have shorter loops than rejected ones, the filter favours short loops. Long loops that many antigens need are dropped before docking.

| CDR3 <= | Pass | Fail | Missed | Rej bad | FDR | Miss rate | FPR | Accuracy |
|---|---|---|---|---|---|---|---|---|
| 2.0 | 1054 | 346 | 288 | 1006 | 0.247 | 0.215 | 0.256 | 0.498 |
| 3.0 | 1187 | 213 | 518 | 776 | 0.152 | 0.304 | 0.215 | 0.633 |
| 4.0 | 1268 | 132 | 756 | 538 | 0.094 | 0.374 | 0.197 | 0.751 |

Separation at 3 angstrom: a design with a good loop has a higher pLDDT with probability 0.851. It has a higher pTM with probability 0.813.

The placed designs are mostly the kept ones. Designs that cannot be placed are mostly rejected, so the miss counts here are a floor for the rejected side.

## Cloud Run refolding of corrected sequences

A Cloud Run L4 GPU job refolded 3000 corrected designs and 100 paired original controls with NanoBodyBuilder2.
Corrected designs achieved a mean USalign TM-score of 0.6958 and mean whole-chain RMSD of 2.7444 angstrom. Alignment used designed VHH scaffold backbones for PDB 8e8z_C (DOI 10.3390/ijms26199262).
Paired original controls achieved a mean USalign TM-score of 0.6936 and mean whole-chain RMSD of 2.7252 angstrom.
Zero designs out of 3100 passed the common 2.0 angstrom (TM >= 0.90), < 1.5 angstrom, or strict 1.0 angstrom self-consistency bars.
The third loop was placed in 2665 corrected designs, and 148 (5.55%) passed the 3.0 angstrom CDR3 bar.
Paired original controls achieved 4 passes out of 93 placed loops (4.30%).
