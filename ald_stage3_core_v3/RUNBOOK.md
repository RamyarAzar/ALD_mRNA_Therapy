# RUNBOOK — ALD Stage-3 Core (Mouse Branch)

This runbook explains how to install the dependencies, configure the model, and
run a baseline simulation using the Stage-3 Python core.

---

## 1. Environment & Requirements

- Python ≥ 3.10
- Recommended packages:
  - `numpy`
  - `scipy`
  - `pyyaml`
  - `pydantic` (optional, for config validation)
- OS: Linux / macOS / Windows (any 64-bit system with sufficient RAM).

If you have strict compute constraints (e.g., running on a clinical workstation),
document them here:

- OS:
- CPU:
- RAM:
- GPU (if any):

---

## 2. Installation

From the root of the repository:

```bash
# 1. Create a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate      # on Windows: .venv\Scripts\activate

# 2. Install dependencies
pip install -r requirements.txt
```

Example `requirements.txt`:

```text
numpy
scipy
pyyaml
```

(You may add `pydantic` or other tooling as needed.)

---

## 3. Configuration Files

All model inputs are passed via configuration files in `configs/`.

### 3.1. Parameters (`example_params_mouse.yaml`)

This file defines:

- **ModelConfig**: platforms, compartments, cell types.
- **BlockAParams**: PBPK volumes, flows, clearance, BBB/FUS parameters, uptake parameters.
- **BlockBParams**: η_esc, α_kc,cell, k_deg, k_tr, peroxisomal insertion, Vβ parameters.
- **BlockCParamsMouse**: VLCFA production/clearance, inflammation and axon injury parameters.
- **BlockDParams**: safety induction/decay rates and driver function choices.
- **BlockEParams**: CMC metrics and mappings Φ, Ψ, Θ.

### 3.2. Initial State (`example_initial_state.yaml`)

Initial conditions for:

- Carrier concentrations `C_kc`
- Cumulative uptake `U_kc,cell`
- `m`, `P_tot`, `P_peri`
- `C26_plasma`, `C26_CNS`, `InflammCNS`, `AxonInjury`
- `Scytokine`, `S_ALT`, `S_AST`, `ADA_k`

### 3.3. Doses (`example_doses_mouse.yaml`)

Defines a list of dose events:

- `time` [hours]
- `platform` (e.g., `"LNP"` or `"EV"`)
- `compartment` (e.g., `"plasma"` for IV)
- `amount` (nominal dose per kg)
- `batch_index` (to link to CMC metrics)

Also contains `t0` and `tend` for the simulation window.

---

## 4. Running the Baseline Example

Once your configs are set:

```bash
python -m examples.example_run_mouse \
  --params configs/example_params_mouse.yaml \
  --initial configs/example_initial_state.yaml \
  --doses configs/example_doses_mouse.yaml
```

This will:

1. Load config files.
2. Build `ModelConfig` and block parameter objects.
3. Construct `ALDModel`.
4. Build initial state vector.
5. Run `.simulate((t0, tend), y0, doses)`.
6. Compute derived outputs.
7. Save results into `outputs/` as:

   - `sample_run_times.csv`
   - `sample_run_states.npy`
   - `sample_run_derived.json`

---

## 5. Interpreting the Outputs

- **Time series (`sample_run_times.csv`, `sample_run_states.npy`)**
  - You can reconstruct trajectories of any state from these arrays using the
    `StateIndexLayout` defined in `ald_model/core.py`.

- **Derived outputs (`sample_run_derived.json`)**
  - `DeltaC26_plasma_rel` and `DeltaC26_CNS_rel`  
    Relative change in VLCFA at final time compared to baseline.
  - `NfL_mouse`, `FA_mouse`, `Rotarod_mouse`  
    Mouse endpoints derived from axonal injury.

These are the primary metrics used in Stage-2 for Aim-2 and Aim-3-style questions
(VLCFA correction and functional benefit).

---

## 6. Verification Workflow

For Stage-3 acceptance:

1. **Traceability Matrix**  
   Use `verification/TRACEABILITY_MATRIX.md` to link each model block and formula
   to its implementation in `ald_model/core.py`.

2. **Acceptance Tests**  
   Implement and run the tests described in `verification/ACCEPTANCE_TESTS.md`.
   Record pass/fail status for each commit or release.

3. **Scenario Equivalence**  
   In `verification/SCENARIO_EQUIVALENCE.md`, document the comparison between
   Stage-2 reference scenarios and this Stage-3 baseline run.

4. **Version Links**  
   Maintain `verification/VERSION_LINKS.md` to keep Stage-1, Stage-2, and Stage-3
   documents tied together.

Once these documents are completed, this Stage-3 package is ready to be published
(e.g., Zenodo, internal registry, or client delivery).
