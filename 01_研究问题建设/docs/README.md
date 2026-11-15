# SRS joint snapshot v2 — static code stage, 未运行

This directory is independent of historical static Gaussian MAT/checkpoints and
`chasi_paired` results. The formal v2 entry is `run_srs_joint.m`; the copied
`run_campaign.m` remains a legacy reference and must not label v2 output.
All commands below are **待用户运行**. No numerical code was executed here.

## What is implemented

`pilot/srs_joint_snapshot.m` defines a zero-mean joint Gaussian snapshot for
two observed UL ports mapped to rows 1 and 2 of a four-row DL channel. Each
TRP transmit column uses the same pilot mixing matrix M but independent
Gaussian columns. At measurement time, `Y=M H0+N`; the conditional mean is
`Q Mᴴ(M Q Mᴴ+σ²I)⁻¹Y`, and the full per-column posterior block matrix is
`Q-Q Mᴴ(M Q Mᴴ+σ²I)⁻¹M Q`. User cross blocks are retained. With independent
innovation of covariance Q, a snapshot at age Δ uses `a=exp(-Δ/Tc)`,
`H_true=a H0+sqrt(1-a²) innovation`, `H_hat=a H_hat0`, and
`C=a² C0+(1-a²)Q`. This exp correlation is a **project assumption**, not
Jakes and not a consistent multiple-time trajectory. Caller must supply
`H0`, innovation and Q from the same Gaussian prior. No CSI-RS observation
is fused. The two pilot indices belong to **port resources**; two orthogonal
port sequences per UE require 16 sequence resources if eight UE sequence
labels are used. The pilot index input makes that distinction explicit.

`python/srs_joint.py` is an array-level port of the MATLAB definition.
`pilot/export_joint_fixture.m` and `python/compare_frozen.py` compare exact
frozen arrays, not integer seeds. Both are pending user execution.
`prepare_srs_joint_input.m` generates all-user SRS observations on the copied
scene's geometry, schedule and two-TRP service masks. It discards the copied
scene's old H/C and diagnostics. All 210 UEs send SRS in every group while
only `resource.ueIds` are downlink scheduled. Port 1 uses pilot resources
1–8 and port 2 uses 9–16. The default prior gain now uses the source
scene's full UE×TRP long-term gain; the former distance exponent is available
only as `distance_power_law_ablation`. The data records the full SRS sender IDs,
pilot map, prior Rx blocks, noise, age and resource weights. If the source
scene supplies only one active service TRP, v2 adds the nearest other TRP
to enforce the fixed dual-TRP main experiment; this adapter choice needs
server audit.

The exact local/server command sequence and expected artifacts are in
`SERVER_RUNBOOK.md`. `compare_frozen.py` checks only the observation model;
`compare_solver_parity.py` is the six-branch WA/ADMM oracle comparison.
The propagation and conditional-MSE assumptions are audited in
`PROPAGATION_AND_MODEL_AUDIT.md` (all numerical checks NOT_RUN).

`run_srs_joint.m` calls the existing MATLAB matrix WMMSE–ADMM with six explicit
branches. `python/native_solver.py` is a copied matrix solver adapted by the
new `python/campaign.py` to the same branch semantics. Both use marginal C in
the receiver and precoder quadratic. Under independent TRP priors and noises,
same-UE cross-TRP conditional error blocks are zero; same-pilot cross-UE
blocks are nonzero but independent data symbols remove them from individual
UE second moment. The compact cross representation is priorRx + pilotId +
energy + noise; no dense all-UE cross tensor is stored for 210 UE runs.

The robust rate is `-log2 det(E)` from conditional expected MSE, a WMMSE
proxy. It is not claimed to equal expected true SE. PF aggregates each UE's
resource-weighted proxy rate then applies `sum log1p(R_u/floor)`; the local
quadratic uses its derivative. The six branches share initial precoders,
mask, rank, budget and frozen true channel. Each call performs bounded outer
work and writes identity-checked checkpoints. Five consecutive steps require
objective change ≤1e-5, transmit covariance change ≤1e-4 and ADMM residual
convergence at abs=1e-5, rel=1e-4. Iteration cap means unconverged.

## Reference and adaptation map

| Class | Reference/function | Project variable and implementation | Assumption |
|---|---|---|---|
| A | `matlab/pilot_mmse_from_channels.m`, author `main_DNN.py` lines 270–298 in primary repository | shared pilot observation and LMMSE gain; `srs_joint_snapshot.m` | source is single-link Rx covariance form; no KTH equation is claimed |
| B | joint Gaussian conditioning identity | full `crossUE`, 4×4 `C`, 2→4 mapping in MATLAB/Python | independent users and transmit columns; cross-TRP blocks absent |
| C | target scene `chasi_generate_scene.m` | exp age coefficient, 1 ms project TTI, selection rows | snapshot conditional model, not physical reciprocity calibration |
| D | robust advantage over nominal | pending six-branch campaign | hypothesis only |

The KTH thesis text was not located locally. Its §5.5, §5.6 and §10.1
remain unverified here; no empirical KTH coefficients have been used.

## Historical static input audit

`chasi_prepare_target_static_input.m` calls `chasi_generate_scene`, then
`chasi_target_static_scene_inputs.m` copies `scene.resources` and replaces
only `H`, `C`, `C_stat`. This leaves scene `age_ms` and `srs_diagnostic`
present despite the static label. `chasi_measure_target_group_static.m`
computes `pilotNoiseVariance=median(home prior gain)/10^(pilotSnrDb/10)`;
the scene generator's `srs_diagnostic.pilot_noise_variance` is produced by
its separate aged measurement path. They are different diagnostics. The
historical input is untouched. Old static `pilot_sumse/pilot_pf` set solver
covariance to zero. `chasi_run_target_static_batch.m` and the aged runner
now select `H_true` for `ideal_pf`; the source guard is
`pilot/test_branch_identity_static.m` (pending run).

## Model and solver limits to resolve before a formal campaign

* The current explicit all-user SRS sender assumption should be replaced if
  the actual SRS activity mask differs. The old scene's group-local set was
  not used as the SRS sender set.
* Cross-TRP blocks are zero only under this version's independent TRP prior
  and noise assumption. Correlated site shadowing in large-scale gain does
  not itself induce conditional small-scale error cross blocks.
* The copied solver's robust objective is a weighted expected MSE proxy, not
  an established exact expected true SE. Its `C_stat` enters receiver and
  quadratic update, but this must be checked on a frozen solver state.
* The copied legacy experiment uses `outerRelTol=1e-4` and
  `outerPrecoderRelTol=5e-3`; the v2 batch runner uses the new configured
  criteria. Numerical equivalence and convergence are pending server checks.

## Pending execution

Run in order on the server after uploading the whole `srs_joint_v2` directory:

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
matlab -batch "addpath('pilot'); validate_srs_joint_small"
matlab -batch "run_srs_joint('input')"
matlab -batch "run_srs_joint('export')"
python python/prepare_input.py output_srs_joint_v2/frozen_binary --report output_srs_joint_v2/input_import.json
python python/validate_contract.py output_srs_joint_v2/frozen_binary
python python/compare_frozen.py path/to/fixture.mat path/to/reference.mat
python python/campaign.py output_srs_joint_v2/frozen_binary output_srs_joint_v2/python_checkpoints --branch all
python python/status.py output_srs_joint_v2/python_checkpoints
python python/evaluate.py output_srs_joint_v2/frozen_binary output_srs_joint_v2/python_checkpoints output_srs_joint_v2/evaluation_python.json
python python/plot_zh.py output_srs_joint_v2/evaluation_python.json output_srs_joint_v2/plots_python --manifest output_srs_joint_v2/frozen_binary/manifest.json
```

The parity fixture/reference must separately be exported with
`pilot/export_joint_fixture.m` using explicit same arrays. For MATLAB solves,
repeat `matlab -batch "run_srs_joint('solve','all')"` until
`matlab -batch "run_srs_joint('status')"` reports convergence or cap;
then run `run_srs_joint('evaluate')` and `run_srs_joint('plot')`. In tmux,
`tmux new -s srs_joint` then run the solve command with `2>&1 | tee -a
output_srs_joint_v2/solve.log`; detach with Ctrl-B D, inspect with `tail -f`.
No scheduler or background process was started here.
