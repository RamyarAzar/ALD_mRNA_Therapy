# Scenario Equivalence — Stage-2 vs Stage-3 Baseline

**Objective:** Show that the Stage-3 Python core reproduces key outputs of the
Stage-2 Baseline scenario within predefined tolerances.

---

## Baseline Scenario Definition

- Platforms: LNP (primary), EV (optional booster).
- Route: IV to plasma.
- Dose schedule (example):
  - Day 0: 1× LNP dose
  - Day 1: FUS window (if applicable)
  - Day 2: 1× LNP dose
  - Day 2: EV booster (optional)
- FUS: single FUS window at 24 h (parameters set in `blockA`).

The exact configuration used for this Baseline is encoded in:

- `configs/example_params_mouse.yaml`
- `configs/example_doses_mouse.yaml`
- `configs/example_initial_state.yaml`

---

## Comparison Metrics (Example: Day 7 or Day 28)

At a chosen evaluation time (e.g., t = 168 h or 672 h):

| Metric                  | Stage-2 Reference | Stage-3 Core | Tolerance | Pass? |
|-------------------------|-------------------|-------------:|----------:|:-----:|
| ΔC26_plasma_rel        | ref_ΔC26_pl       | core_ΔC26_pl |    5%     |       |
| ΔC26_CNS_rel           | ref_ΔC26_CNS      | core_ΔC26_CNS|    5%     |       |
| NfL_mouse              | ref_NfL           | core_NfL     |   10%     |       |
| FA_mouse               | ref_FA            | core_FA      |   10%     |       |
| Rotarod_mouse          | ref_Rotarod       | core_Rotarod |   10%     |       |

- `Stage-2 Reference` values come from the locked Formula Pack / scenario tables.
- `Stage-3 Core` values come from `sample_run_derived.json`.

---

## Conclusion

When all metrics fall within the agreed tolerances, this document can be cited
as evidence that the Stage-3 implementation is **scenario-equivalent** to
the Stage-2 Baseline for the mouse branch.
