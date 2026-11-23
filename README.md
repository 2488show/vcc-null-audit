# vcc-null-audit

**Empirical null distributions and an adversarial audit of the metrics used to score
single-cell perturbation-prediction models (Arc Institute Virtual Cell Challenge).**

> Status: **work in progress** — stage 0 and the metric audit (round 1 plus two follow-up
> checks) are done. Not yet done: multiple seeds, a full-data reference, validation/test
> targets, and a comparison with real model submissions.

## The question

Virtual-cell models are scored against two reference points: 0 = "predict the control
mean", 1 = "a split-half replicate". The 0 point is a *single prediction with no
variance*, while a real perturbation with no effect still produces a noisy measured
change. This project asks:

1. **What does "no effect" look like as a distribution?** We build an empirical null
   from control cells, matched on cell count, and compare it to the analytic
   noise correction that the 2026 metrics already use (`expr_mse_unbiased`).
2. **Which official metrics can be raised by predictions that contain no biology?**
   Seven adversarial predictions (no-change, read thinning, global-mean scaling,
   matched noise, shuffled labels, random rotations, and "knock down only the target
   gene") are scored with the official `cell-eval2` scorer.

Prior work showed that simple baselines are hard to beat (see `Dar-kSun/baseline-first`
and the 2026 challenge wrap-up). This project is about the *measuring stick* rather
than the models.

## Stage 0 results (2025 training split)

Data: 221,273 H1 hESC cells × 18,080 genes; 150 CRISPRi targets + 38,176
non-targeting controls; 48 batches.

| Check | Pre-registered prediction | Result |
|---|---|---|
| Hand-computed `expr_mae` vs `cell-eval2` | match | match to 2.8e-17 |
| Targets not expressed in H1 controls | enough to build a second null | **0 of 150** — prediction void |
| Variance explained by the shared response direction | 40–60% | **12.4%** — prediction wrong |
| Does the empirical null remove the cell-count confound? | yes | Spearman with cell count −0.545 → **+0.074** |
| Same, replicated on validation / test | yes | −0.610 → −0.243 (p=0.09) / −0.490 → +0.053 |

Why the two predictions failed (both are part of the result):

- The training split contains **no** ribosome, spliceosome, proteasome, Pol II or
  Integrator targets, so the expected shared "growth-arrest" axis is weak by
  construction.
- Raw effect size mostly reflects **how many cells a perturbation has**: median ‖δ‖
  falls from 10.25 (243 cells) to 2.76 (2,282 cells). Without the null, this would
  have read as biology.

Null-distribution checks:

- Fit `null² = a/n + b`; held-out error +0.29% / −0.41%.
- Resampling from a finite control pool inflates the null by up to **31%** at large n;
  the fitted, theoretical and measured pool term agree and are subtracted.
- Every perturbation sits above the null (min signal-to-noise 1.15). **Batch** was
  ruled out (each perturbation spans a median 48/48 batches). The **on-target
  knockdown itself** explains it: all 150 targets are repressed (median 5.9% expression
  remaining), and the target's own gene carries a median 10% of ‖δ‖². Excluding it
  lowers the minimum to 1.04.

All 300 targets across the three splits are expressed in H1 controls, and none of the
three splits contains ribosome, Pol II, proteasome or Integrator targets.

The three released files share **the same 38,176 control cells** (identical barcodes and
counts). The dataset therefore has 414,694 unique cells, not the 491,046 obtained by
adding the files.

## Metric audit, round 1

Eleven predictions that contain no (or the wrong) biology were scored with
`cell-eval2 --preset vcc2026` on the official scale (0 = generic-response baseline,
1 = split-half replicate). Reference: all 150 training perturbations, 50 cells each,
1,000 controls — **a subsample, not the official full-data score**; single seed.

| Prediction | PDS | MSE | NMAE | FIDELITY | REACH | JACCARD | AVG |
|---|---|---|---|---|---|---|---|
| Control bootstrap (no change) | 0.10 | 0 | −0.10 | −0.04 | −0.04 | 0.02 | −0.01 |
| Reads thinned to 25% | 0.10 | 0 | −0.16 | **0.61** | −0.03 | −0.09 | **0.07** |
| Multiplicative noise | −0.03 | 0 | −0.11 | **0.54** | −0.05 | −0.11 | 0.04 |
| Real effects, shuffled labels | 0.03 | 0 | −0.09 | **0.57** | −0.10 | −0.10 | 0.05 |
| Mean response ×1 (sanity check) | −0.01 | 0 | −0.02 | 0.04 | −0.01 | −0.08 | −0.01 |
| Mean response ×2 | 0.05 | 0 | 0.05 | **0.36** | 0.11 | −0.13 | **0.07** |
| Mean response + target knock-down | −0.06 | 0 | −0.04 | 0.02 | −0.05 | −0.02 | −0.02 |

- **Refuted hypothesis:** knocking down the target gene from the label alone earns
  nothing — the 2026 metrics exclude the on-target gene, as documented.
- **FIDELITY rewards indiscriminate calling.** `direction_fidelity_yield_raw` is
  `k / max(n_pred, n_conf)`, where `k` counts sign matches over every gene the
  prediction calls significant (the reference only needs a non-zero log fold change)
  and is not chance-corrected. Calling thousands of genes yields a raw value near 0.5
  (baseline 0.03, replicate 0.82), i.e. ~0.6 on the scale.
- **The MSE clamp hides the cost.** These predictions have 1.6–2.9× the baseline's raw
  MSE, but the score is clamped at 0, so their average ends up positive.
- PDS and JACCARD are **not reported**: across three independent seeds, a no-change
  prediction's PDS spans +0.10 / +0.05 / −0.20, and the junk predictions' JACCARD spans
  0.16–0.18 — both as large as any effect they could show. See the seed table below.

### Follow-up checks

**The chance-corrected metric penalises the same predictions.** `cell-eval2` also ships
`direction_fidelity_yield` (chance-corrected); the `vcc2026` preset scores the `_raw`
variant. On the same reference (median 17 significant genes per real perturbation):

| Prediction | Genes called significant (median) | Raw (scored) | Corrected |
|---|---|---|---|
| Control bootstrap | 0 | 0.00 | +0.04 |
| Mean response ×1 | 0 | 0.06 | −0.11 |
| Multiplicative noise | 443 | 0.46 | −2.97 |
| Shuffled labels | 394 | 0.47 | −2.87 |
| Reads thinned to 25% | 4,641 | 0.51 | −3.02 |

**The effect persists with 4× the cells.** Scale rebuilt on 40 perturbations × 200 cells
(2,000 controls):

| Prediction | FIDELITY, 50 cells | FIDELITY, 200 cells | AVG, 200 cells |
|---|---|---|---|
| Reads thinned to 25% | 0.61 | 0.59 | +0.08 |
| Multiplicative noise | 0.54 | 0.58 | +0.11 |
| Shuffled labels | 0.57 | 0.43 | −0.01 |
| Mean response ×1 (sanity check) | 0.04 | 0.18 | +0.01 |
| Mean response + target knock-down | 0.02 | 0.14 | 0.00 |

### Three independent seeds

Each seed re-draws the reference subsample, rebuilds the scale (baseline emission seed and
anchor split seed) and regenerates the predictions. FIDELITY, `from_replicate`:

| Prediction | seed 0 | seed 1 | seed 2 | mean | range |
|---|---|---|---|---|---|
| Reads thinned to 25% | 0.610 | 0.632 | 0.632 | **0.62** | 0.022 |
| Multiplicative noise | 0.542 | 0.580 | 0.571 | **0.56** | 0.038 |
| Shuffled labels | 0.566 | 0.565 | 0.556 | **0.56** | 0.009 |
| Mean response ×1 (sanity check) | 0.037 | 0.023 | 0.014 | 0.03 | 0.023 |
| Mean response + target knock-down | 0.024 | 0.049 | 0.049 | 0.04 | 0.026 |
| Control bootstrap | −0.035 | −0.004 | −0.014 | −0.02 | 0.031 |

The FIDELITY result holds: 0.56–0.62 for the three junk predictions, each stable to within
0.04 across seeds. The sanity check stays near zero all three times, and the
target-knock-down prediction never beats the baseline (average −0.02 / −0.01 / +0.01).

Average score is positive but small — about +0.03 for the junk predictions, with a
per-prediction seed range of 0.03–0.06, so only its sign is worth reporting.

Two caveats from the 200-cell run. The sanity check drifts on FIDELITY (0.18 instead of
~0), probably because our integer re-emission of the mean response is called significant
more often than the official baseline's (not verified); read the two "mean response" rows
against that offset. The thinned and noisy predictions also beat the baseline on
JACCARD there (+0.16), but the seed table shows JACCARD moving by 0.16–0.18 on its own, so
that reversal cannot be separated from noise — whether JACCARD offsets the FIDELITY gain is
simply unresolved here.

### A real model for comparison

State ST (`arc-state` 0.11.1) trained from scratch on this data — 200 perturbations,
32,474 cells, 2,000 HVGs, `output_space=all`, no pretrained weights — then scored on a
**validation-split** reference (50 perturbations the ground-truth cells of which the model
never saw), alongside the same junk predictions on the same scale:

| Prediction | FIDELITY | REACH | PDS | AVG | Genes called significant (median) |
|---|---|---|---|---|---|
| **State ST, 40,000 steps** | **0.87** | **0.93** | 0.83 | **0.63** | **851** (mean 1,565) |
| State ST, 2,000 steps | 0.09 | 0.09 | 0.08 | 0.01 | 0 (mean 1.8) |
| Reads thinned to 25% | **0.64** | −0.12 | 0.16 | **0.07** | 4,641 |
| Shuffled labels | **0.53** | −0.04 | 0.11 | **0.06** | 394 |
| Multiplicative noise | **0.54** | −0.03 | −0.00 | **0.04** | 443 |
| Mean response ×1 (sanity check) | −0.02 | −0.07 | 0.05 | −0.01 | 0 |
| Control bootstrap | −0.01 | −0.09 | −0.03 | −0.04 | 0 |

Real data has a median of 15 significant genes per perturbation (mean 487).

**Training to the official 40,000 steps changes the picture** (local Apple-silicon GPU, 20 h,
~0.55 steps/s; the 2,000-step row above is the earlier, under-trained run). The fully trained
model scores 0.63 overall against 0.01, calls a mean of 1,565 genes significant (median 851;
real data: mean 487, median 15) and its FIDELITY rises from 0.09 to 0.87. That is the
prediction made from the 2,000-step run: "talking little" came from under-training, and once the
model calls hundreds of genes the metric rewards it. It now also beats all three junk
predictions on every headline score (0.63 vs 0.04–0.07).

Two things keep this from being a clean win. First, the model now **over-calls** by about 3×
(mean 1,565 vs 487), the same direction as the Challenge entries Arc described, and the
chance-corrected FIDELITY is **negative** (−1.84; the 2,000-step model: −0.63), so the corrected
metric still penalises the volume. Second, the sign agreement on called genes is 0.73 (junk
predictions: 0.51), which is real signal, but the scored FIDELITY mixes signal and volume and
this experiment cannot separate how much each contributes. So: a well-trained model beats noise
here, the gap between "model" and "junk" comes mainly from training, and the call-volume
mechanism of the metric still stands.

**Caveats for the 40,000-step numbers.** (i) The model emitted 1,813 predicted values above the
largest value in the training data (log1p 14.0 against 6.83), which made 380 cells exceed the
scorer's count limit; predictions were clipped to the training range (8,476 values clipped, with
`ST_CLIP=6.8272` in `scripts/29_score_state.py`) before scoring. (ii) One training run, one
seed, no repeat. (iii) Scores of PDS 0.83 and NMAE 1.09 are close to or above the replicate
level (1.0 on this scale); this is the cell-held-out setting, not the Challenge's unseen
perturbations, so the model has seen each of these perturbations in training.

The setting is also easier than the Challenge: State uses one-hot perturbation encoding and
cannot predict unseen perturbations without a gene-embedding file we could not reach, so the
held-out unit is cells, not perturbations.

Full notes: [`NOTES.md`](NOTES.md).

## Reproduce

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python \
    anndata numpy scipy pandas polars h5py cell-eval2
# data: gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/ -> data/raw/
python scripts/00_hard_gate.py          # needs the cell-eval2 repo's 5 MB tutorial subset in vendor/
python scripts/01_pass1_group_sums.py   # one streaming pass, ~15 s
python scripts/02_facts_1_and_3.py
python scripts/04_null_n2.py 160        # resumable; run again until all grid points are done
python scripts/05_null_model.py
python scripts/06_batch_matched.py
python scripts/07_exclude_target.py
```

Runs in under 3 GB RAM (streaming over the 15.5 GB CSR file).

## Data and tools

- Roohani et al. *Virtual Cell Challenge: Toward a Turing test for the virtual cell.* Cell (2025).
  Data via the Arc Virtual Cell Atlas.
- `cell-eval2` (Arc Institute) — official scorer.
- Replogle et al. *Mapping information-rich genotype-phenotype landscapes with
  genome-scale Perturb-seq.* Cell (2022).

## License

MIT
