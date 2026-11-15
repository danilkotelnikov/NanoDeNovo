<!-- version: w06 | date: 14.XI.2026 -->
# Quality report, dataset 20

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

- Designs that pass: 1382
- Designs that fail: 1615
- pLDDT, all designs: n=2997, mean=0.8870, median=0.9196, sd=0.1079, min=0.3103, max=0.9813
- pTM, all designs: n=2997, mean=0.8019, median=0.8443, sd=0.1325, min=0.1537, max=0.9168
- pLDDT, designs that pass: n=1382, mean=0.9482, median=0.9504, sd=0.0166, min=0.9011, max=0.9813
- pTM, designs that pass: n=1382, mean=0.8801, median=0.8803, sd=0.0161, min=0.8501, max=0.9168

## Designs that pass the filter but fail the structure check

- Designs that pass the filter and also have a structural measurement: 1382
- Correct by the reference standard: 1340
- Failures: 42
- Failure share of the measured set: 0.030
- TM-score of the measured set: n=1382, mean=0.9574, median=0.9638, sd=0.0240, min=0.8431, max=0.9903
- RMSD of the measured set: n=1382, mean=0.9686, median=0.8900, sd=0.3480, min=0.3900, max=2.1900

### How much the structural bar decides

A nanobody has a long and flexible third loop. A strict bar fails a design whose framework is right and whose loop moved, so the bar is shown at four settings.

| Bar | TM-score at or above | RMSD at or below | Failures | Share |
|---|---|---|---|---|
| strict | 0.95 | 1.0 | 528 | 0.382 |
| strict match, relaxed RMSD | 0.95 | 2.0 | 418 | 0.302 |
| common design bar | 0.9 | 2.0 | 42 | 0.030 |
| loose | 0.85 | 3.0 | 1 | 0.001 |

The ten weakest failures at the chosen bar:

| Design | pLDDT | pTM | TM-score | RMSD |
|---|---|---|---|---|
| 0070_top.1.20 | 0.9088 | 0.8649 | 0.8431 | 2.19 |
| 0148_top.3.20 | 0.9054 | 0.8679 | 0.8500 | 2.06 |
| 0575_top.1.20 | 0.9017 | 0.8528 | 0.8561 | 1.81 |
| 0070_top.3.20 | 0.9040 | 0.8572 | 0.8597 | 2.05 |
| 0258_top.3.20 | 0.9020 | 0.8508 | 0.8659 | 1.57 |
| 0636_top.2.20 | 0.9355 | 0.8654 | 0.8685 | 1.76 |
| 0318_top.2.20 | 0.9272 | 0.8536 | 0.8688 | 1.49 |
| 0532_top.3.20 | 0.9310 | 0.8685 | 0.8707 | 1.61 |
| 0916_top.1.20 | 0.9280 | 0.8508 | 0.8720 | 1.86 |
| 0932_top.1.20 | 0.9344 | 0.8632 | 0.8743 | 1.95 |

## False negatives

Every one of the 1615 rejected designs has a structural measurement, so the filter can be scored in both directions. A miss is a design that the reference calls correct and the filter rejected.

| Bar | Pass | Fail | Missed | Rej bad | FDR | Miss rate | FPR | Accuracy |
|---|---|---|---|---|---|---|---|---|
| strict (0.95/1.0) | 854 | 528 | 384 | 1231 | 0.382 | 0.310 | 0.300 | 0.413 |
| strict match, relaxed RMSD (0.95/2.0) | 964 | 418 | 435 | 1180 | 0.302 | 0.311 | 0.262 | 0.467 |
| common design bar (0.9/2.0) | 1340 | 42 | 908 | 707 | 0.030 | 0.404 | 0.056 | 0.750 |
| loose (0.85/3.0) | 1381 | 1 | 1191 | 424 | 0.001 | 0.463 | 0.002 | 0.858 |

Failures among kept is the false discovery rate. Miss rate is the share of correct designs that the filter threw away. Wrong designs kept is the false positive rate. All three are measured against structural self-consistency, which is a calculation and not a laboratory result.

Separation at the chosen bar: a correct design scores higher pLDDT than a wrong one with probability 0.942, and higher pTM with probability 0.905. A value of 0.5 means the score tells nothing; 1.0 means it separates them perfectly. The whole chain is dominated by the conserved framework, so this value is generous. Read the third-loop check for the harder test.

## Threshold choice

Each row shows how many measured designs a threshold pair keeps, and how many of those are false positives.
Every threshold pair is scored on all designs. Rows are ordered by Youden's J, the share of correct designs kept minus the share of wrong designs kept. J is 0 for a filter that keeps everything, however common correct designs are.

| pLDDT above | pTM above | Kept | False positives | Missed | Failure share | Miss rate | J |
|---|---|---|---|---|---|---|---|
| 0.9 | 0.7 | 1929 | 86 | 405 | 0.045 | 0.180 | 0.705 |
| 0.9 | 0.75 | 1929 | 86 | 405 | 0.045 | 0.180 | 0.705 |
| 0.88 | 0.7 | 2199 | 157 | 206 | 0.071 | 0.092 | 0.699 |
| 0.88 | 0.75 | 2199 | 157 | 206 | 0.071 | 0.092 | 0.699 |
| 0.9 | 0.8 | 1901 | 83 | 430 | 0.044 | 0.191 | 0.698 |
| 0.88 | 0.8 | 2063 | 147 | 332 | 0.071 | 0.148 | 0.656 |
| 0.85 | 0.75 | 2422 | 253 | 79 | 0.104 | 0.035 | 0.627 |
| 0.85 | 0.7 | 2458 | 263 | 53 | 0.107 | 0.024 | 0.625 |
| 0.92 | 0.7 | 1489 | 24 | 783 | 0.016 | 0.348 | 0.620 |
| 0.92 | 0.75 | 1489 | 24 | 783 | 0.016 | 0.348 | 0.620 |
| 0.92 | 0.8 | 1489 | 24 | 783 | 0.016 | 0.348 | 0.620 |
| 0.85 | 0.8 | 2113 | 183 | 318 | 0.087 | 0.141 | 0.614 |

The study's own cut (pLDDT above 0.9, pTM above 0.85) scores J = 0.540. It keeps 1382 designs, keeps 42 wrong ones and misses 908.

The best row sits on the edge of the grid. A wider grid could do better, so this table shows a direction and not an optimum.

J ignores what the later stages cost. A looser filter sends more designs to docking and simulation. Pick a row with that cost in mind.

## Third-loop check

A whole-chain match is dominated by the conserved framework, so most designs pass it. This check fits the framework alone and measures how far the third loop (CDR3) moved. The loop is found by a sequence heuristic, so some designs cannot be placed.

- Designs placed: 2638 of 2997
- Kept designs placed: 1357 of 1382
- Rejected designs placed: 1281 of 1615

| Measure | Kept | Rejected |
|---|---|---|
| Median framework RMSD (angstrom) | 0.73 | 0.89 |
| Median CDR3 RMSD (angstrom) | 1.44 | 3.12 |
| Median CDR3 length (residues) | 9 | 12 |

Loop length matters. If kept designs have shorter loops than rejected ones, the filter favours short loops. Long loops that many antigens need are dropped before docking.

| CDR3 <= | Pass | Fail | Missed | Rej bad | FDR | Miss rate | FPR | Accuracy |
|---|---|---|---|---|---|---|---|---|
| 2.0 | 891 | 466 | 459 | 822 | 0.343 | 0.340 | 0.362 | 0.512 |
| 3.0 | 1028 | 329 | 617 | 664 | 0.242 | 0.375 | 0.331 | 0.624 |
| 4.0 | 1152 | 205 | 767 | 514 | 0.151 | 0.400 | 0.285 | 0.727 |

Separation at 3 angstrom: a design with a good loop has a higher pLDDT with probability 0.806. It has a higher pTM with probability 0.725.

The placed designs are mostly the kept ones. Designs that cannot be placed are mostly rejected, so the miss counts here are a floor for the rejected side.

## Cloud Run refolding of corrected sequences

A Cloud Run L4 GPU job refolded 3000 corrected designs and 100 paired original controls with NanoBodyBuilder2.
Corrected designs achieved a mean USalign TM-score of 0.6880 and mean whole-chain RMSD of 2.7961 angstrom. Alignment used designed VHH scaffold backbones for PDB 8e8z_C (DOI 10.3390/ijms26199262).
Paired original controls achieved a mean USalign TM-score of 0.6834 and mean whole-chain RMSD of 2.8013 angstrom.
Zero designs out of 3100 passed the common 2.0 angstrom (TM >= 0.90), < 1.5 angstrom, or strict 1.0 angstrom self-consistency bars.
The third loop was placed in 2608 corrected designs, and 145 (5.56%) passed the 3.0 angstrom CDR3 bar.
Paired original controls achieved 2 passes out of 95 placed loops (2.11%).
