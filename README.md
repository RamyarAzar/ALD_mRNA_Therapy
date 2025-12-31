# ALD Dual-Platform mRNA Therapy — Stage 3 Core (Mouse Branch)

> **Status:** Research prototype • **Focus:** Mouse branch (Blocks A–E) of a dual-platform mRNA therapy model for X-linked Adrenoleukodystrophy (ALD).  
> **Owner:** Method2Model — Method-to-Model computational studio.

This repository contains the **Stage-3 Python implementation** of a mechanistic model
for a dual-platform (LNP + EV) mRNA therapy program in **X-linked Adrenoleukodystrophy (ALD)**.

The goal of Stage 3 is:

> **“Deliver a reproducible tool — and prove it matches the locked logic.”**

That means:

- The approved Stage-2 formulas and I/O contract are implemented as a **professional,
  reusable Python package**.
- We provide **evidence** that the code matches the Stage-1 architecture and Stage-2
  Formula Pack (traceability matrix, acceptance tests, scenario equivalence).

⚠️ **Important:** This repository is a **research modeling prototype**, not a medical
device, diagnostic, or clinical decision support system.

---

## Project Overview

The model represents a dual-platform mRNA program (LNP + EV) and follows the
architecture defined in **Stage 1** and formalized in the **Stage-2 Formula Pack**:

- **Block A — Delivery & Exposure**
  - Coarse-grained PBPK with 5 compartments: plasma, liver, periphery, CNS, CSF.
  - Saturable, surface-modified uptake into cell types (e.g., hepatocytes, glia).
  - Time-dependent blood–brain barrier flux with FUS window.

- **Block B — Expression & Peroxisomal Targeting**
  - Endosomal escape → cytosolic mRNA → ALDP protein.
  - Partitioning into peroxisomal vs. non-peroxisomal pools.
  - Dynamic peroxisomal insertion and turnover.
  - Effective peroxisomal β-oxidation capacity **Vβ(t)** as weighted sum of cell types.

- **Block C (Mouse PD) — Disease / Progression**
  - Very-long-chain fatty acid (VLCFA, C26:0) dynamics in plasma and CNS.
  - CNS inflammation and axonal injury.
  - Mouse functional endpoints: **NfL, FA, Rotarod**, derived algebraically from axon injury.

- **Block D — Safety & Immunogenicity**
  - Cytokine response, hepatic injury markers (ALT/AST), and anti-drug antibodies (ADA_k).
  - Flexible driver functions `f_cyt`, `f_hep`, `f_ADA` driven by PBPK exposure.

- **Block E — CMC / Batch Variability**
  - Batch-level CMC metrics `Q_{k,b}` mapped to algebraic factors:
    - `φ_k,b(Q)` for dose scaling (wired into bolus application).
    - `ψ_k,b(Q)` and `θ_mRNA,b(Q)` reserved placeholders for future tuning of `η_esc` and `kdeg_m`.

The Stage-3 core is designed to be:

- **Reproducible** — all inputs are configuration-driven (YAML).
- **Traceable** — every formula is linked to code via a traceability matrix.
- **Reviewable** — architecture ⇄ formulas ⇄ code are linked via version documents.

---

## Repository Structure

```text
ald_stage3_core/
  README.md           # High-level project description (this file or a close variant)
  RUNBOOK.md          # Installation and execution guide
  requirements.txt    # Minimal Python dependencies

  ald_model/
    core.py           # Main model implementation (Blocks A–E, mouse branch)
    __init__.py       # Public exports for the package

  configs/
    example_params_mouse.yaml     # PBPK, PD, safety, CMC parameter template
    example_initial_state.yaml    # Initial conditions
    example_doses_mouse.yaml      # Dose & FUS schedule

  examples/
    example_run_mouse.py          # End-to-end example script

  verification/
    TRACEABILITY_MATRIX.md        # Architecture → formulas → code mapping
    ACCEPTANCE_TESTS.md           # Test plan and status
    SCENARIO_EQUIVALENCE.md       # Stage-2 vs Stage-3 baseline comparison
    VERSION_LINKS.md              # Version IDs across Stage 1 / 2 / 3

  outputs/
    (generated outputs from example runs)
```

---

## Installation

We recommend Python ≥ 3.10 and a virtual environment.

```bash
git clone [https://github.com/RamyarAzar/ALD_mRNA_Therapy.git]
cd ALD_mRNA_Therapy

python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

Minimal `requirements.txt`:

```text
numpy
scipy
pyyaml
```

Feel free to add `pydantic` or other tooling for config validation and testing.

---

## Quickstart (Mouse Baseline Example)

1. Inspect / edit the example configs in `configs/`:

   - `example_params_mouse.yaml`
   - `example_initial_state.yaml`
   - `example_doses_mouse.yaml`

2. Run the example script:

```bash
python -m examples.example_run_mouse \
  --params configs/example_params_mouse.yaml \
  --initial configs/example_initial_state.yaml \
  --doses configs/example_doses_mouse.yaml
```

3. After running, you should see outputs in `outputs/`:

- `sample_run_times.csv` — time grid
- `sample_run_states.npy` — state trajectories
- `sample_run_derived.json` — derived metrics (ΔC26, NfL, FA, Rotarod, etc.)

These derived metrics correspond to the **Stage-2 endpoints** for mouse branch
(Aim-2 / Aim-3: VLCFA correction and functional benefit).

---

## Scope & Limitations

To avoid ambiguity, the scope of this Stage-3 core is intentionally narrow and explicit:

- ✅ **Implemented**
  - Mouse branch only (Blocks A–E).
  - Coarse-grained PBPK (5 compartments).
  - BBB/FUS dynamics and saturable uptake.
  - VLCFA, CNS inflammation, axonal injury, mouse endpoints.
  - Safety block (cytokines, ALT/AST, ADA_k).
  - CMC dose scaling via `φ_k,b(Q)`.

- ⚠️ **Not implemented in this repository**
  - Human projection branch (human NfL, Loes score, etc.).
  - Explicit complement activation state (currently folded into driver functions).
  - Full wiring of `ψ_k,b` and `θ_mRNA,b` into `η_esc` and `kdeg_m`.
  - Scenario engine, VOI analysis, and power simulations (Blocks F & G) — to be built on top.

This is deliberate: the current repository focuses on the faithfully implemented
core that can be reviewed, tested, and extended, without silently expanding scope.

---

## Verification Artifacts

To support **“Proof of Match”** between architecture, formulas, and code, we include:

- `verification/TRACEABILITY_MATRIX.md`  
  For each block and formula, shows exactly where it lives in `ald_model/core.py`.

- `verification/ACCEPTANCE_TESTS.md`  
  Test definitions (mass balance, monotonicity, CMC scaling, etc.) and pass/fail status.

- `verification/SCENARIO_EQUIVALENCE.md`  
  How to compare the Stage-3 baseline scenario outputs to the frozen Stage-2 references.

- `verification/VERSION_LINKS.md`  
  Links and version IDs for:
  - Stage-1 architecture
  - Stage-2 Formula Pack & I/O Contract
  - Stage-3 Python core (this repo)

These documents are meant to be **filled in** as part of internal or external review.

---

## How to Cite

If you use this code or concepts in a publication or internal report, please cite:

> Azar R., *ALD Dual-Platform mRNA Therapy — Stage-2 Formula Pack and Stage-3 Core Model.*  
> Method2Model, 2025. Zenodo, DOI: **10.5281/zenodo.18093206** (Stage-2 Formula Pack).  
> GitHub: `<[https://github.com/RamyarAzar/ALD_mRNA_Therapy.git]>` (Stage-3 Core).

You can also reference Method2Model as:

> Method2Model — Method-to-Model computational studio (https://method2model.com)


---

## License

### Short Version

- **Allowed:** research and educational use, inspection, and internal prototyping.
- **Not allowed without permission:** commercial use, redistribution as part of a paid product,
  or positioning this code as a clinical / regulatory-grade tool.

### Full License Text (Method2Model Research License 1.0)

Copyright (c) 2025 Method2Model and contributors.

Permission is hereby granted, free of charge, to any person or organization obtaining a copy
of this software and associated documentation files (the "Software"), to use, reproduce,
and modify the Software **solely for research and educational purposes**, subject to the
following conditions:

1. **Non-Commercial Use Only.**  
   The Software may not be used, in whole or in part, for commercial purposes without prior
   written permission from Method2Model. “Commercial purposes” include, but are not limited to,
   use in or as part of a paid product or service, consulting engagements, or any activity that
   generates direct or indirect revenue.

2. **No Clinical or Regulatory Claims.**  
   The Software is a research prototype and must not be used as a medical device, diagnostic,
   treatment recommendation system, or as the basis for regulatory submissions without
   appropriate validation and authorization. The authors and Method2Model make **no** warranties
   regarding fitness for clinical or regulatory use.

3. **Attribution.**  
   Any public use (e.g., publications, internal slides, non-commercial demos) must credit:
   “Method2Model — ALD Dual-Platform mRNA Therapy Stage-3 Core” and, where applicable,
   provide a link to the GitHub repository and relevant Zenodo record(s).

4. **Redistribution.**  
   Redistribution of modified or unmodified versions of the Software is permitted only if:
   - The use remains non-commercial, and  
   - This license text is included in full, and  
   - Any modifications are clearly indicated.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.

For commercial licensing or collaboration inquiries, please contact:

- **Method2Model — Method-to-Model computational studio**  
- Email: `drazar@method2model.com`

---

## Contact

For questions, collaboration, or commercial licensing:

- **Name:** Dr. Ramyar Azar  
- **Role:** Founder, Method2Model  
- **Email:** `drazar@method2model.com`  
- **Website:** https://method2model.com

If you open issues or pull requests on GitHub, please avoid including any
patient-identifiable or proprietary clinical data.
