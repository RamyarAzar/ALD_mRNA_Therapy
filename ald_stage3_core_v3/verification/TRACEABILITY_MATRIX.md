# Traceability Matrix — Stage-3 Core (Mouse Branch)

| Block | Concept / Formula                          | Stage-1 / Stage-2 Source                       | Code Implementation (`ald_model/core.py`)           | Test ID |
|-------|--------------------------------------------|------------------------------------------------|------------------------------------------------------|--------|
| A     | PBPK compartments, flows, clearance       | Arch §A; Formula eq. (1)                       | `ALDModel._rhs`: `dC_dt` loop over compartments      | T-A01  |
| A     | BBB flux with FUS window                  | Arch §A; Formula eqs. (3)–(4)                  | `_P_BBB`, `_rhs` (flux term plasma↔CNS)             | T-A02  |
| A     | Saturable uptake `v_uptake`               | Formula eq. (5)                                | `_v_uptake`                                         | T-A03  |
| B     | mRNA input `k_in` from uptake             | Formula eq. (7)                                | `_rhs`: `kin_cell` (`η_esc * α * v_uptake`)         | T-B01  |
| B     | mRNA / P_tot ODEs                         | Formula eqs. (8)–(9)                           | `_rhs`: `dm_dt`, `dP_tot_dt`                        | T-B02  |
| B     | Peroxisomal insertion & P_peri dynamics   | Formula eqs. (10)–(11)                         | `_rhs`: `dP_peri_dt`                                | T-B03  |
| B     | Vβ(t) via PALDP_peri                      | Formula eqs. (13)–(14)                         | `_compute_Vbeta_per_cell`, `_compute_Vbeta_global`  | T-B04  |
| C(m)  | VLCFA plasma & CNS dynamics               | Formula eqs. (16)–(17)                         | `_rhs`: `dC26_plasma_dt`, `dC26_CNS_dt`             | T-C01  |
| C(m)  | Inflammation & axon injury                | PD spec; ϕ_infl, Faxon                         | `_rhs`: `dInflamm_dt`, `dAxon_dt`                   | T-C02  |
| C(m)  | Endpoints NfL, FA, Rotarod (mouse)        | Formula eqs. (21)–(23)                         | `ALDModel.compute_derived_outputs`                  | T-C03  |
| D     | Cytokine / ALT / AST ODEs                 | Formula Safety section                         | `_rhs`: `dS_cyt_dt`, `dS_ALT_dt`, `dS_AST_dt`       | T-D01  |
| D     | ADA_k ODEs                                | Formula Safety section                         | `_rhs`: `dADA_dt`                                   | T-D02  |
| E     | CMC mapping φ_k,b                         | Formula CMC section (Φ, Ψ, Θ)                  | `BlockEParams.effective_scale_and_escape`           | T-E01  |
| E/A   | Dose scaling with φ_k,b                   | CMC → PBPK coupling                            | `ALDModel._apply_bolus_doses_at_time` (Deff / V_c)  | T-E02  |

Notes:

- Human branch, complement state, و wiring کامل ψ/θ در این Stage-3 core پیاده نشده
  و در `README.md` به‌صورت صریح در Scope ذکر شده است.
