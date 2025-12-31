# ALD Dual-Platform mRNA Therapy Model — Stage 3 Core (Mouse Branch)

> **Goal (Stage 3):** Deliver a reproducible tool — and prove it matches the locked logic.

This repository implements the approved Stage-2 formulas for the **ALD dual-platform mRNA
therapy program** as a professional, reusable Python package.

The code in `ald_model/core.py` is a direct implementation of the Stage-2 Formula Pack and
Stage-1 architecture (Blocks A–E, mouse branch):

- **Block A — Delivery & Exposure**
  - Coarse-grained PBPK with 5 compartments: plasma, liver, periphery, CNS, CSF.
  - Saturable, surface-modified uptake into cell types (hepatocyte / glia, etc.).
  - Time-dependent blood–brain barrier flux with FUS window.

- **Block B — Expression & Peroxisomal Targeting**
  - Endosomal escape → cytosolic mRNA → ALDP protein.
  - Partitioning into peroxisomal vs. non-peroxisomal pools.
  - Dynamic peroxisomal insertion and degradation.
  - Effective peroxisomal β-oxidation capacity `Vβ(t)` as a weighted sum across cell types.

- **Block C (Mouse) — Disease / PD**
  - Very-long-chain fatty acid (VLCFA, C26:0) dynamics in plasma and CNS.
  - Inflammation in CNS and axonal injury.
  - Mouse functional endpoints: NfL, FA, Rotarod, derived algebraically from axon injury.

- **Block D — Safety & Immunogenicity**
  - Cytokine response, hepatic injury markers (ALT/AST), and anti-drug antibodies (ADA_k).
  - All driven by PBPK exposure and flexible driver functions (`f_cyt`, `f_hep`, `f_ADA`).

- **Block E — CMC / Batch Variability**
  - Batch-level CMC metrics `Q_{k,b}` mapped to algebraic factors:
    - `φ_k,b(Q)` for dose scaling (wired into bolus application).
    - `ψ_k,b(Q)` and `θ_mRNA,b(Q)` reserved for future use on `η_esc` and `kdeg_m`.

The model is packaged as a **long-term digital asset**:
you can run it locally or in the cloud, plug in your own parameter tables and dose schedules,
and generate reproducible simulation outputs and decision metrics.

---

## Scope of This Stage-3 Core

To make the scope explicit for reviewers:

- ✅ **Implemented**
  - Mouse branch (Blocks A–E) with coarse-grained PBPK.
  - BBB/FUS dynamics, saturable uptake, expression, peroxisomal targeting.
  - VLCFA / inflammation / axonal injury and mouse endpoints.
  - Safety block (cytokines, ALT/AST, ADA) driven by flexible exposure functions.
  - CMC dose scaling via `φ_k,b` in `BlockEParams`.
  - Robust piecewise ODE integration with bolus jumps (`solve_ivp`).

- ⚠️ **Not implemented in this core (reserved for future versions)**
  - Human projection branch (NfL_human, Loes score, etc.).
  - An explicit complement activation state (currently folded into `f_cyt` / `f_hep`).
  - Full wiring of `ψ_k,b` and `θ_mRNA,b` into `η_esc` and `kdeg_m`.
  - Scenario engine and VOI / power analysis (Blocks F & G) — to be built on top.

These constraints are deliberate so Stage-3 remains tightly matched to the approved Stage-2
Formula Pack and I/O contract, without silently extending the scope.

---

## Repository Layout

- `ald_model/core.py`  
  Main model implementation (`ALDModel`, `ModelConfig`, block parameter dataclasses,
  `DoseEvent`, simulation and derived endpoints).

- `configs/`  
  YAML templates for:
  - PBPK / PD / safety / CMC parameters (`example_params_mouse.yaml`)
  - Initial state (`example_initial_state.yaml`)
  - Dose and FUS schedule (`example_doses_mouse.yaml`)

- `examples/example_run_mouse.py`  
  Minimal runnable example:
  1. Load parameter and dose configs.
  2. Instantiate `ALDModel`.
  3. Integrate over time.
  4. Save state trajectories and derived outputs.

- `verification/`  
  Evidence that the code is faithful to the locked logic:
  - `TRACEABILITY_MATRIX.md`
  - `ACCEPTANCE_TESTS.md`
  - `SCENARIO_EQUIVALENCE.md`
  - `VERSION_LINKS.md`

- `outputs/`  
  Sample run results for a Baseline scenario (e.g., standard LNP dose + default FUS)
  will be placed here when you execute the example script.

---

## Stage-3 Inputs and Outputs

### Inputs (Stage-3 Requirements)

- **Sample dataset**  
  Synthetic or de-identified example parameter set and dose schedule, encoded as YAML/CSV
  in `configs/`. This defines:
  - PBPK parameters (volumes, flows, clearances).
  - PD parameters (production / clearance / Hill parameters).
  - Safety parameters (induction / decay rates and drivers).
  - CMC metrics for each batch (used by Block E).
  - Dose events and FUS timing.

- **Compute specs (optional)**  
  If you plan to run locally on a constrained machine, you may document your OS/CPU/RAM/GPU
  in `RUNBOOK.md` for reproducibility.

### Outputs (Stage-3 Deliverables)

- **Python Code Package / Repo (versioned)**  
  This repository itself, with semantic or git-based versioning.

- **Runbook + Config Templates**  
  See `RUNBOOK.md` for installation and execution, and `configs/` for plug-and-play
  config templates.

- **Sample Run Outputs**  
  In `outputs/`, we provide (after running `examples/example_run_mouse.py`):
  - `sample_run_times.csv`
  - `sample_run_states.npy`
  - `sample_run_derived.json` (ΔC26, NfL, FA, Rotarod, safety metrics).

- **Verification Evidence (“Proof of Match”)**
  - `TRACEABILITY_MATRIX.md`: Architecture → Formula Pack → Code components.
  - `ACCEPTANCE_TESTS.md`: Baseline pass/fail tests.
  - `SCENARIO_EQUIVALENCE.md`: Baseline scenario outputs vs. Stage-2 reference values.
  - `VERSION_LINKS.md`: Version IDs linking Stage-1 architecture ⇄ Stage-2 formulas ⇄ Stage-3 code.

---

## Relationship to Stage-1 and Stage-2

- **Stage-1 (Architecture)** defines the high-level block structure (A–G) and flows.
- **Stage-2 (Formula Pack & I/O Contract)** defines the dynamic equations and input/output
  schema. The current Stage-2 Formula Pack for this project is archived at:

  - DOI: `10.5281/zenodo.18093206`

- **Stage-3 (this repository)** is the code implementation of the mouse branch for Blocks A–E,
  aligned with those documents.

Version links and change control are documented in `verification/VERSION_LINKS.md`.

---

## License / Use

This Stage-3 package is provided primarily for **research and educational** purposes.

- You may:
  - Inspect, run, and adapt the code for internal, non-commercial research.
  - Use the example configs and outputs as templates for your own feasibility checks.

- You must **not**:
  - Redistribute or integrate this package into commercial tools or paid services
    without explicit permission.
  - Present this package as a stand-alone medical product or clinical decision tool.

For commercial use, redistribution, or integration into a paid deliverable, please contact:

- **Method2Model** — Method-to-Model computational studio  
- Email: `drazar@method2model.com`

By using this repository, you acknowledge that it is a **modeling prototype** only and
not a substitute for regulatory-grade safety or efficacy evaluation.
