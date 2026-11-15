# Code and result audit, NanoDeNovo

Audit date: 14.XI.2026. Every finding below was checked against the original files on disk,
while corrected wrappers now live under `pipeline/` and `tools/`.

The published paper is Int. J. Mol. Sci. 26(19):9262, DOI 10.3390/ijms26199262.

## Summary

Thirteen findings. One of them changed which designs the study carried forward.
Four of them stop the committed scripts from reproducing the committed results.
Two are about the study design rather than the code. The rest are smaller.

| # | Where | Severity | Effect |
|---|---|---|---|
| 1 | `process_020.py`, `process_030.py` | Changes results | Kept the three weakest sequences of each batch of 128 |
| 2 | `align_pipeline.py` | Latent | Comment and code disagree on the match score; no number changed |
| 3 | `tfold_ab_20.py`, `tfold_ab_30.py` | Loses data | Stopped at backbone 999, so three designs per dataset are absent |
| 4 | `run_proteinmpnn_v_030.sh` | Stops the run | Writes to one folder and reads from another |
| 5 | `process_0*.py` vs `tfold_ab_*.py` | Stops the run | Writes to `tfold_020`, next step reads `tfold_20` |
| 6 | `align_pipeline.py` | Cannot reproduce | Writes to `usalign`, but the results on disk are in three other folders |
| 7 | `rfantibody_wrapper2.py` | Stops at once | The same option is declared twice, so the script cannot start |
| 8 | `tfold_ab_*.py` | Silent corruption | Adds rows to the table on every run, so a second run doubles them |
| 9 | `tfold_ab_*.py` | Fragile | Reads the scores from fixed line numbers of the structure file |
| 10 | Study design | Changes the reading | Survivors only were checked. Now measured: the filter drops 30 to 40 in 100 designs with a good structure, and favours short loops |
| 11 | `run_proteinmpnn_v_0*.sh` | Slow | Loads the model again for each of 1000 backbones |
| 12 | Study design | Weakens the claim | Only 1 of 838 interfaces reaches a confidence of 0.8 |
| 13 | `mmpbsa.in` | Weakens the claim | Energies have no entropy term, one run each, and the winner changes with the method |

## 1. The sequence choice was inverted

`process_020.py` line 51 and `process_030.py` line 51:

```python
sorted_entries = sorted(parsed_entries, key=lambda x: -x[0])
top_3 = sorted_entries[:3]
```

ProteinMPNN writes `score=` in every sequence header. That score is the mean
negative log-likelihood of the sequence. A low score means the model finds the
sequence likely on that backbone. The minus sign sorts the scores downwards, so
`top_3` holds the three highest scores, which are the three least likely
sequences.

Evidence from the files. Backbone 0001 of dataset 20 has 128 sampled sequences.
`pmpnn_020/seqs/0001_20.t.0.15.fa` contains samples down to `score=0.8832`.
`parsed_sequences_20.txt` records the chosen three as 0.9751, 0.9745 and 0.9717,
which are the three highest values in the file.

Size of the effect, measured over all 2000 batches:

| Dataset | Mean score of the chosen three | Mean score of the best three | Difference |
|---|---|---|---|
| 20 | 1.0257 | 0.8979 | 0.1278 |
| 30 | 1.0793 | 0.9489 | 0.1304 |

The difference is per residue. The designs are about 117 residues long.

Everything after this step used the chosen sequences: the refolding, the
confident sets, the docking, the simulations and the energy scores. The designs
in the paper are therefore not the designs the method was meant to produce.
They are still real sequences with real measured scores. The claim that needs
care is the claim that they are the best the pipeline found.

Fix: `tools/select_sequences.py`. It sorts upwards, ignores the input record,
and writes a record of every choice with its rank inside its batch.

## 2. The match score is mislabelled, and it changed no number

`align_pipeline.py`, `parse_usalign_output`:

```python
tm_matches = re.findall(r'TM-score=\s*([\d.]+)', output)
# Use the last TM-score which is normalized by Structure_2 (scaffold)
tm = tm_scores[-1] if len(tm_scores) > 0 else None
```

USalign is called as `USalign <scaffold> <design>`, so the scaffold is
structure 1 and the design is structure 2. USalign prints the score normalised
by structure 1 first and by structure 2 second. The code takes the last value,
which is normalised by the design, while the comment says it is normalised by
the scaffold. The code and the comment disagree.

I first rated this as a change to the numbers. The data say otherwise. The two
scores differ only when the design and the backbone differ in length. All 5994
designs have exactly the length of their backbone, and the largest gap between
the two scores is 0.0. So no published number depends on this. The defect is
latent. It will bite the first time a design has a different length.

One more fact from the disk. The saved text reports, such as
`0004_usalign.txt`, list the design as structure 1. The committed script lists
the scaffold first. So the script in the folder is not the code that wrote those
reports. That supports finding 6.

Fix: `tools/align_all.py` puts the design first, reads both scores by name and
keeps both in the table.

## 3. The folding loop stopped one backbone early

`tfold_ab_20.py` line 23: `for i in range(1, 1000)`. In Python that counts 1 to
999. `process_020.py` writes sequences for backbones 1 to 1000.

Evidence: `folded_sequences_20.csv` holds 2997 design rows. Three backbones
times 999 is 2997. Backbone 1000 is absent from both datasets.

## 4. The second launcher writes and reads different folders

`run_proteinmpnn_v_030.sh`:

```bash
--out_folder ./pmpnn_030
...
mv ./pmpnn_30/seqs/gen_${num}.fa ./pmpnn_30/seqs/${num}_30.t.0.15.fa
```

The output goes to `pmpnn_030` and the rename reads `pmpnn_30`, which does not
exist. As committed, the rename fails on every one of the 1000 turns. The
dataset 20 launcher does not have this defect.

## 5. The folder name changes between two steps

`process_020.py` writes to `tfold_020/`. `tfold_ab_20.py` reads from
`tfold_20/`. On disk only `tfold_20` exists. The same mismatch is in the
dataset 30 pair.

## 6. The alignment step cannot reproduce the alignment results

`align_pipeline.py` writes into a folder named `usalign` and filters at
RMSD below 2.0 with a match score above 0.90. The results on disk are in
`usalign_rmsd_1.0_tm_0.95`, `usalign_rmsd_1.5_tm_0.95` and
`usalign_rmsd_2.0_tm_0.90`. The committed script produces only the third of
these three. The thresholds that produced the other two are not in the code.

`align_pipeline_kwargs.py` takes the thresholds as options and is the version
that matches the folders. Keep that one and delete the fixed-threshold copy.

## 7. The RFantibody wrapper cannot start

`rfantibody_wrapper2.py` declares `--config` twice:

```python
initial_parser.add_argument('--config', help='Path to config JSON file')
...
parser.add_argument('--config', help='...', action='store_true')
```

The second declaration is on the main parser, which also inherits nothing from
the first, so the option means two different things. The main parser then
receives `remaining_args`, which no longer holds `--config`. Two further
problems sit next to this one. `required=not config` makes every argument
required whenever no configuration file is given, which contradicts the
defaults. And `Tee` claims live logging, but `subprocess.run` returns only
after the command ends, so the log appears in one block at the end.

## 8. The score table grows on every run

`tfold_ab_*.py` writes the header only when the file is absent, then appends a
row for each design. A second run on the same folder appends a second copy of
every row. The current tables hold no repeated rows, so this has not happened
yet.

## 9. The scores are read from fixed line numbers

`tfold_ab_*.py`:

```python
plddt = float(lines[1].split(':')[-1].strip())
ptm = float(lines[3].split(':')[-1].strip())
```

Any change in the header of the structure file moves the values and the script
reads the wrong number without complaining. `tfold_ag_dock_kwargs.py` does this
correctly: it searches for the `REMARK 250` lines by name. Use that approach.

## 10. The filter was never tested against anything independent

`align_pipeline.py` reads `confident_proteins_20.txt`, which holds only the
designs that already passed the confidence filter. So the structure check was
run on the survivors alone.

Evidence: `folded_sequences_20.csv` has 2997 designs, 1382 of them pass the
filter, and `usalign_20.csv` has exactly 1382 rows.

The study therefore could not say how many good designs the filter threw away.
This audit closes that gap. `tools/align_all.py` aligned all 5994 designs, the
rejected ones included. The new run reproduces every old row: 1382 of 1382 in
dataset 20 and 1417 of 1417 in dataset 30, with the match score equal to five
decimals and the RMSD equal to two. The old tables are therefore sound as
tables. Then `tools/qc_report.py` scored the filter in both directions.

### The whole-chain reference

The reference is structural self-consistency: a design is correct when its
refolded model matches its own backbone. Dataset 20, at four bars:

| Bar | Match >= | RMSD <= | Pass rate | Pass | Fail | Missed | FDR | Miss rate |
|---|---|---|---|---|---|---|---|---|
| strict | 0.95 | 1.0 | 0.413 | 854 | 528 | 384 | 0.382 | 0.310 |
| strict match, relaxed RMSD | 0.95 | 2.0 | 0.467 | 964 | 418 | 435 | 0.302 | 0.311 |
| common design bar | 0.90 | 2.0 | 0.750 | 1340 | 42 | 908 | 0.030 | 0.404 |
| loose | 0.85 | 3.0 | 0.858 | 1381 | 1 | 1191 | 0.001 | 0.463 |

Dataset 30 at the common bar: 1378 kept and correct, 39 kept and wrong, 921
missed, 0.028 failures among kept, 0.401 miss rate.

Read the third column first. At the common bar, 75 in 100 of all designs pass.
That bar is a weak test. A whole-chain match is dominated by the conserved
nanobody framework, and almost any nanobody-like model matches it. So the 3 in
100 failures at that bar flatter the filter. The strict bar is the opposite
case: 38 in 100 of the kept designs fail it, and 31 in 100 of the correct
designs are lost. No single bar describes the filter. Quote the bar with every
share.

### The third loop

The decision that matters is the third loop (CDR3), the long flexible one that
touches the antigen. `tools/cdr3_check.py` fits the framework alone and
measures how far the loop moved. It finds the loop with a sequence heuristic, so
10 to 12 designs in 100 cannot be placed. Most of those are rejected designs.

| | Dataset 20 kept | Dataset 20 rejected | Dataset 30 kept | Dataset 30 rejected |
|---|---|---|---|---|
| Median loop deviation (angstrom) | 1.44 | 3.12 | 1.15 | 3.49 |
| Median loop length (residues) | 9 | 12 | 9 | 12 |

At a loop deviation of 3 angstrom or less, 62 in 100 designs pass. Then 24 in
100 of the kept designs in dataset 20 fail, and 15 in 100 in dataset 30. The
miss rate is 38 in 100 and 30 in 100. So the filter does keep better loops. It
also keeps many designs whose loop moved, and it throws away about a third of
the designs whose loop held. The score separates good loops from bad loops with
probability 0.81 for pLDDT in dataset 20 and 0.85 in dataset 30. The whole-chain
figure of 0.94 was too generous.

The filter also favours short loops. Kept designs have a median loop of 9
residues and rejected designs 12. Many antigens need a long loop. A filter that
drops long loops before docking narrows what the pipeline can find.

### What drives the misses

A threshold sweep on all designs, at the common whole-chain bar, shows the
trade. Keep pLDDT above 0.90 and move pTM from above 0.85 to above 0.70. The
misses fall from 908 to 405. The wrong designs kept rise from 42 to 86. The
filter then keeps 1929 designs, not 1382. Youden's J, the share of correct
designs kept minus the share of wrong designs kept, rises from 0.540 to 0.705.
So the pTM cut-off costs most of the misses. The best row sits on the edge of
the grid, so this shows a direction and not an optimum, and a looser filter
sends more designs to the costly stages.

### Limits

- The reference is a calculation, not a laboratory result. The refolded model,
  pLDDT and pTM all come from one folding model. A low-confidence model that
  lands near the backbone may be a lucky guess. The miss rate is therefore a
  count of designs the reference calls correct, and it does not show that they
  would fold in a cell.
- The loop is placed by a heuristic, and the loop bar of 3 angstrom is a choice.
  The unplaced designs are mostly rejected ones, so the miss counts are a floor.
- The designs measured here are the ones finding 1 produced, the three least
  likely sequences of each batch. A corrected selection could change every rate.
  Testing that needs the folding model and a graphics card.
- Whether a rejected design with a good loop also binds the antigen is not
  measured. Nobody has docked the rejected designs.

## 11. The model is loaded 1000 times

`run_proteinmpnn_v_0*.sh` starts `protein_mpnn_run.py` once per backbone, so
the weights are read from disk and moved to the card 1000 times per dataset.
ProteinMPNN accepts a folder of structures in one call. One call per dataset
removes 999 model loads.

## 12. The predicted interfaces are weak

The complex stage predicted each nanobody together with the antigen chain
`8e8z_C` and reported an interface confidence, ipTM, for each one.

| Dataset | Complexes | Mean interface confidence | At or above 0.7 | At or above 0.8 |
|---|---|---|---|---|
| 20 | 838 | 0.297 | 11 (1.3 %) | 1 (0.1 %) |
| 30 | 905 | 0.310 | 16 (1.8 %) | 2 (0.2 %) |

Almost nine complexes in ten sit below 0.5. The subset the study called
confident has a mean of 0.52 and a lowest value of 0.36.

A predicted interface is normally trusted from about 0.6 upwards, and treated
as probably wrong below that. The exact floor for a nanobody and antigen pair
is not settled, and interface predictors are known to score true antibody
interfaces low, so this is a caution rather than a verdict.

The practical point stands. The complexes that went forward to the simulations
and the energy scoring carry interface confidences in a range where the
interface geometry itself may be wrong. A binding free energy computed on a
wrong interface is a precise number about the wrong arrangement. Any claim from
that stage should carry the interface confidence of the complex it came from.

## 13. The energy scores cannot rank the designs

Where: `md/<complex>/gromacs/UNIGBSA/FINAL_RESULTS_MMPBSA.dat` and
`mmpbsa.in` in each of the five simulated complexes.

Four facts from those files.

- The three entropy options in `mmpbsa.in` are all 0, and no result file has an
  entropy line. So the totals are not binding free energies. They run from
  about -13 to -70 kcal/mol. A binder with an affinity of 1 nM has a free energy
  near -12 kcal/mol at room temperature, from RT times ln(1e-9). The limit
  applies only to totals without conformational entropy, and: a review of small-ligand
  binding, with an emphasis on calibration, lists that lack among several crude
  approximations, as an expert opinion, not a new measurement
  [consensus:w:cb9ce8cb-30cb-4142-af1a-60d861167b90#1].
- Each complex has one trajectory. The spread between frames is 9 to 20
  kcal/mol. The three candidates differ by 11 kcal/mol in the PB totals: -37.66,
  -40.11 and -48.62. With no repeat run, nobody knows the run-to-run spread.
- The two solvation methods disagree. The same complex differs by up to 25
  kcal/mol between GB and PB, and the winner changes with the method.
- The unhumanised designs score stronger than their humanised versions in both
  methods, in both pairs that were simulated. Humanising costs 13 to 52
  kcal/mol in these totals. That is more than the 11 kcal/mol that separates
  the three candidates.

| Complex (folder) | GB total | PB total |
|---|---|---|
| `scFv-0389-304-6` | -65.35 | -70.45 |
| `scFv-0389-304-6-humanized` | -12.93 | -37.66 |
| `scFv-0389-459-5` | -51.48 | -67.05 |
| `scFv-0389-459-5-humanized` | -37.96 | -40.11 |
| `scFv-0743-166-6` | -32.10 | -48.62 |

Among the three humanised candidates, GB ranks `scFv-0389-459-5-humanized` first
and PB ranks `scFv-0743-166-6` first. The study reported the PB order.

What the literature supports, read as abstracts only. On 46 protein-protein
complexes, the best MM/GBSA setting gave a correlation of -0.647 with measured
affinity [doi:10.1039/c6cp03670h#1]. The authors chose that setting from several
force fields, solvation models and dielectric constants on the same data. The
complexes are not nanobody complexes. The correlation holds across different
complexes, not for ranking among close variants of one binder.

On 43 protein-protein systems with ZDOCK decoys, with structures relaxed or
simulated by the authors' protocol and abstract read only, MM/GBSA rescoring
separated the correct structures from the decoys better than the docking score
did [doi:10.1039/c6cp03670h#2].

I found no study in five sources that ranks close nanobody designs by this
energy. With pose quality judged against a reference structure, and clearest at
elevated temperature, a 2026 paper reports that heated molecular dynamics with
replicates separates near-native antibody-antigen complexes from incorrect ones
by the stability of the interface
[consensus:w:71d027b2-a0f8-4dd4-9034-0a8a4a922625#2]. That is a check that a
pose holds, and not an energy.

So the evidence supports simulation for one job here: checking that a predicted
interface holds. That support is weak, because the 2026 paper used heated
simulation and the study ran an unheated one. The energies carry some signal,
since both methods rank each parent above its humanised version. I found no
support for ranking close nanobody designs by them, and none for quoting them as
binding free energies. Report the old totals as scores without entropy.

I did not check whether the five trajectories kept their interface. The folders
hold RMSF files and no interface RMSD, and that is the measure the stability job
needs.

I suggest a cheaper step before any simulation: co-folding confidence, clash and
contact checks against the parent pose, and fast rescoring. The real test is a
binding assay. Keeping this stage requires an experimental affinity benchmark. See the Nano-2 specification.

## What is fixed in this folder

| Finding | State |
|---|---|
| 1 | Fixed in `tools/select_sequences.py`, with a check that fails if the sort turns round again |
| 13 | Written up from the old files. No code change. The decision on the stage is open |
| 3, 10, 12 | Measured. `tools/align_all.py` aligns every design, `tools/cdr3_check.py` measures the third loop, and `tools/qc_report.py` scores the filter in both directions |
| 2 | Latent. Fixed in `tools/align_all.py`, which reads both scores by name |
| 4, 5, 6, 7, 8, 9, 11 | Written up here with the correct form; the rewritten steps are staged |

## The independent check of this audit

An independent cross-family review checked the findings and new code on 14.XI.2026. It
agreed that the sequence score direction is as stated and that the pattern
`\bscore=` does not match inside `global_score=`. It changed three things.

- The share of failures is a false discovery rate, not a false positive rate.
  The wording above and in the tool was corrected.
- The strict structural bar is too tight for a nanobody with a flexible third
  loop. The four-bar table above exists because of this, and the headline
  reading changed from "a third of the designs fail" to "3 in 100 fail at the
  bar the field normally uses".
- The selection code could keep the same sequence twice when sampling returned
  it twice. The code now keeps one copy, and the check covers it.

The review marked its own statements about field-standard bars and interface
confidence floors as unverified memory. Those two numbers still need a dated
search before they go into anything published.

A second independent review, on the same day, covered the new alignment of every
design, the confusion counts, the licence reading and the conclusions drawn
from them. It found the code correct and changed four things.

- The common bar passes 75 in 100 of all designs, so it is a weak test, and the
  whole-chain match hides the third loop. The third-loop check above exists
  because of this, and the claim that the filter "holds up" was withdrawn.
- F1 orders thresholds badly when most designs are correct: keeping every
  design scores an F1 of 0.857 and beats the study's filter at 0.74. The sweep
  now uses Youden's J, which is 0 for a filter that keeps everything.
- The line "at the strict bar the numbers swap" overstated the contrast. It
  was removed.
- A project's own statement of its licence does not confirm a commercial right.
  The registry now records whether the licence text was read, and lists the six
  instruments that rest on a statement alone.

## GPU refolding and structural audit

The corrected sequences underwent live refolding on a Cloud Run L4 GPU instance.
The job ran NanoBodyBuilder2 across 6000 corrected designs and 200 paired controls.
Each design was refolded without knowledge of the target backbone.
Harvested C-alpha traces were aligned with USalign against their designed VHH scaffold backbones for PDB 8e8z_C (DOI 10.3390/ijms26199262).

Corrected designs in dataset 20 achieved a mean USalign TM-score of 0.6880 and mean whole-chain RMSD of 2.7961 angstrom.
Paired original controls in dataset 20 achieved a mean USalign TM-score of 0.6834 and mean whole-chain RMSD of 2.8013 angstrom.
Corrected designs in dataset 30 achieved a mean USalign TM-score of 0.6958 and mean whole-chain RMSD of 2.7444 angstrom.
Paired original controls in dataset 30 achieved a mean USalign TM-score of 0.6936 and mean whole-chain RMSD of 2.7252 angstrom.
Zero designs out of 6200 passed the common 2.0 angstrom (TM >= 0.90), < 1.5 angstrom, or strict 1.0 angstrom self-consistency bars.
Third-loop placement within 3.0 angstrom reached 145 of 2608 placed loops (5.56%) in dataset 20 and 148 of 2665 (5.55%) in dataset 30.
Paired controls reached 2 of 95 placed loops (2.11%) in dataset 20 and 4 of 93 (4.30%) in dataset 30.
Because NanoBodyBuilder2 predicts coordinates from sequence alone, whole-chain USalign metrics stayed close between corrected selections and paired controls.
