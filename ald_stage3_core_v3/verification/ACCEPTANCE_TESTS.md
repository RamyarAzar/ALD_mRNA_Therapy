# Acceptance Test Report — Stage-3 Core (Mouse Branch)

## Test Suite Overview

- **T-A**: PBPK / delivery consistency tests
- **T-B**: Expression & peroxisomal targeting tests
- **T-C**: VLCFA / PD / endpoints tests
- **T-D**: Safety & immunogenicity tests
- **T-E**: CMC dose scaling tests

## Summary Table

| Test ID | Description                                      | Status | Notes |
|---------|--------------------------------------------------|--------|-------|
| T-A01   | Mass balance in PBPK without clearance/uptake    |        |       |
| T-A02   | Symmetry of BBB flux with P_BBB(t)=const         |        |       |
| T-A03   | `v_uptake` → `U_kc,cell` consistency             |        |       |
| T-B01   | `k_in` computed from `v_uptake` and α, η_esc     |        |       |
| T-B02   | m/P_tot ODEs reduce to analytic solution (toy)   |        |       |
| T-B03   | P_peri equilibrium under constant P_tot          |        |       |
| T-C01   | Steady state VLCFA with Vβ=0 vs Vβ>0             |        |       |
| T-C02   | Monotonicity: higher C26_CNS → higher Inflamm    |        |       |
| T-C03   | Endpoints monotone in AxonInjury                 |        |       |
| T-D01   | Cytokine response to single bolus                |        |       |
| T-D02   | ADA_k response to repeated dosing                |        |       |
| T-E01   | φ_k,b = 1 leaves reference dose unchanged        |        |       |
| T-E02   | φ_k,b = 2 doubles ΔC at t_dose (within tol)      |        |       |

## Example: T-E02 (φ scaling)

- **Setup:**
  - Single platform, single compartment.
  - Two runs with identical settings except:
    - Run A: φ_k,b = 1
    - Run B: φ_k,b = 2

- **Expectation:**
  - بلافاصله بعد از زمان دوز `t_dose`، غلظت کمپارتمان در Run B حدوداً ۲ برابر Run A است
    (در حد تلورانس عددی solver).

Fill in the `Status` and `Notes` columns as you execute the tests.  
This document is part of the Stage-3 delivery as **Proof of Match**.
