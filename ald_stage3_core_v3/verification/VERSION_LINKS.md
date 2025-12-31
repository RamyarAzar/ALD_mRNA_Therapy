# Version Links — Architecture ⇄ Formulas ⇄ Code

This document ties together the version IDs of the three key layers:

1. **Stage-1 Architecture**
2. **Stage-2 Formula Pack + I/O Contract**
3. **Stage-3 Python Core (this repository)**

---

## Stage-1 Architecture

- Title: ALD Dual-Platform mRNA Therapy — Stage-1 Architecture  
- Version: v1.0  
- DOI / URL: <ADD_ZENODO_OR_INTERNAL_LINK>

---

## Stage-2 Formula Pack & I/O Contract

- Title: ALD mRNA Therapy — Stage-2 Formula Pack  
- Version: v1.0  
- DOI: `10.5281/zenodo.18093206`  
- I/O Contract PDF: <ADD_LINK_IF_SEPARATE>

---

## Stage-3 Python Core

- Repository: `ald_stage3_core_v2`  
- Version: v0.3.1  
- Core file: `ald_model/core.py`

---

## Change Control Notes

- v0.1.0 → First implementation of Blocks A–D (mouse), φ_k,b not wired.  
- v0.2.0 → Robust time-stepping and coarse-grained PBPK finalization.  
- v0.3.0 → `core_v3` initial packaging as Stage-3 core.  
- v0.3.1 → Repo v2: README License/Use section, RUNBOOK polish, configs + example glue tightened.

Any future changes to the model must update this file and the traceability matrix.
If changes touch the locked logic, a new Stage-2 / Stage-3 alignment review may be required.
