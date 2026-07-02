# Phase-1 Exchange-Symmetry Changes

This folder contains the copied and modified code for the first implementation
stage of monomer exchange symmetry for identical diatoms.

The original source files in the parent directory were not modified.

## Goal

Phase 1 validates the exchange-symmetric channel representation at the S-matrix
level only. Cross-section accumulation is intentionally not used here.

The code still builds the original ordered-channel matrices first, then projects
those matrices into the symmetric monomer-exchange basis:

```text
|a,b>_sym = (|a,b> + |b,a>) / sqrt(2)
```

If `a == b`, the symmetric channel has one component with coefficient `1`.

## Modified Files

- `module_skvp_AtomDiatom.f90`
  - Added exchange-symmetric basis bookkeeping:
    `nsym`, `n_sym_open`, `sym_quant_mat`, `sym_old_idx`, `sym_coeff`,
    `sym_ncomp`, `old_to_sym`, `sym_open_idx`,
    `old_open_pos_to_sym_open`.
  - Added projected matrices and S matrix:
    `mat_M_sym`, `mat_M0_sym`, `mat_M00_sym`, `mat_M10_sym`, `Smat_sym`.

- `skvp_AtomDiatom_sym_phase1.f90`
  - Builds the ordered `quant_mat` exactly as before.
  - Builds the symmetric monomer-exchange basis from ordered channels.
  - Builds the open-channel map in the symmetric basis after
    `solve_target_levels`.
  - Replaces the original cross-section call with an S-matrix-only phase-1
    validation path:
    `basic_aux_mat_calcul`, `potential_mat_calcul`, `make_scatt_mat`,
    ordered S solve, projection, symmetric S solve, comparison output.
  - Adds reusable S-matrix solver `SolveSMatrixGeneric`, following the same
    linear algebra structure as the original `PhaseShift`.
  - Adds output files documenting the channel map and S-matrix comparison.

- `Makefile`
  - Added separate build targets:
    `baseline` and `sym_phase1`.

## Important Interpretation

The new symmetric-basis probabilities should not be compared row-by-row to the
old ordered-channel `proba_all.dat`. A symmetric final channel may contain two
ordered final channels coherently, so the correct check is:

```text
ordered S matrix -> project amplitude into symmetric basis -> probability
```

This is what `phase1_smatrix_compare.dat` reports.

## Test Run

Input used:

- `input.nml` copied from the current working tree.

Commands used:

```sh
make baseline
./baseline > baseline_run.log 2>&1
make sym_phase1
./sym_phase1 > sym_phase1_run.log 2>&1
```

Observed result for the current input:

```text
old ncf      = 121
new nsym     = 66
old n_open   = 121
new n_open   = 66
projected N  = 4488
max |dP|     = 5.3290705182E-15
```

The probability sum from the projected ordered S matrix and the probability sum
from the symmetric-basis S matrix match to roundoff:

```text
sum_ordered_projected = 9.9014618385E-01
sum_sym               = 9.9014618385E-01
```

## Output Files

- `phase1_sym_channel_map.dat`: symmetric channel to ordered-channel components.
- `phase1_sym_open_map.dat`: symmetric open-channel map.
- `phase1_smatrix_compare.dat`: projected ordered probability vs symmetric
  probability for each symmetric open channel.
- `phase1_smatrix_summary.dat`: compact validation summary.
- `baseline_run.log`: baseline run log.
- `sym_phase1_run.log`: modified phase-1 run log.
