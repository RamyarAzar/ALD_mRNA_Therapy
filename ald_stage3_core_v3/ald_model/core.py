"""
core_v3.py — ALD dual-platform mRNA therapy model (Stage-3 core, mouse branch)

This version is a *direct refinement* of core_v2 with the same block structure:

  - Block A: Delivery & Exposure (PBPK / BBB / uptake)
  - Block B: Expression & Peroxisomal Targeting (m, PALDP, Vβ)
  - Block C (mouse): Disease / PD (VLCFA, inflammation, axonal injury + endpoints)
  - Block D: Safety & Immunogenicity (cytokines, ALT/AST, ADA)
  - Block E: CMC / batch-level variability → φ, ψ, θ_mRNA (algebraic; φ wired-in)
  - Blocks F/G: Skeleton via Monte Carlo & decision functions (left external)

**Scope of core_v3 (explicit, matches your Stage-2 Formula Pack):**
  - PBPK is coarse-grained to 5 compartments: plasma, liver, periphery, CNS, CSF.
  - Mouse branch of Block C is fully implemented; human branch is NOT yet implemented.
  - Complement activation does NOT have a separate state; it is absorbed into drivers f_cyt / f_hep.
  - Block E: φ_k,b is fully wired into dose scaling; ψ_k,b and θ_mRNA,b are present but not yet
    plumbed into η_esc or kdeg_m at the ODE level (kept for future extension).
  - Vmax_C26 is kept as a parameter for I/O compatibility but *intentionally not used* in the mouse
    ODEs, exactly as in the Stage-2 Formula Pack (Vβ(t) already plays the “capacity” role).
  - The simulate() stepper has been hardened: no exact `t_stop in dose_times` check; all dose times
    are handled with a numerical tolerance, so floating-point round-off cannot silently skip a dose
    at the final time.

The block structure, state ordering, and overall architecture are preserved relative to v2; only
small internal refinements have been made (Vmax_C26 docs + robust time-stepping).

You still need to wire in *numerical values* and scenario engines outside this core.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
#  Model configuration / index sets
# ---------------------------------------------------------------------------

@dataclass
class ModelConfig:
    """
    Holds the discrete index sets (platforms, compartments, cell types).

    Defaults:
      platforms   = ("LNP", "EV")
      compartments= ("plasma", "liver", "periphery", "CNS", "CSF")
      cell_types  = ("hepatocyte", "oligodendrocyte", "astrocyte", "microglia")
    """

    platforms: Tuple[str, ...] = ("LNP", "EV")
    compartments: Tuple[str, ...] = ("plasma", "liver", "periphery", "CNS", "CSF")
    cell_types: Tuple[str, ...] = ("hepatocyte", "oligodendrocyte", "astrocyte", "microglia")

    def __post_init__(self) -> None:
        self.K: int = len(self.platforms)
        self.C: int = len(self.compartments)
        self.N_cell: int = len(self.cell_types)

        self.platform_index: Dict[str, int] = {name: i for i, name in enumerate(self.platforms)}
        self.compartment_index: Dict[str, int] = {name: i for i, name in enumerate(self.compartments)}
        self.cell_index: Dict[str, int] = {name: i for i, name in enumerate(self.cell_types)}

        # Convenience indices for plasma / CNS as used in Block A/C
        self.idx_plasma: int = self.compartment_index["plasma"]
        self.idx_CNS: int = self.compartment_index["CNS"]


# ---------------------------------------------------------------------------
#  State layout helper (keeps vectorized y aligned with blocks A–D)
# ---------------------------------------------------------------------------

@dataclass
class StateIndexLayout:
    """
    Defines slices into the 1D state vector y(t) for all dynamic states.

    State ordering:

      [A] PBPK (Block A)
        - C_kc:  shape (K, C)          → slice_C
        - U_kc_cell: (K, C, N_cell)    → slice_U

      [B] Expression (Block B)
        - m[cell]:          (N_cell,)  → slice_m
        - P_tot[cell]:      (N_cell,)  → slice_P_tot
        - P_peri[cell]:     (N_cell,)  → slice_P_peri
        (Vβ is algebraic, not part of y)

      [C] Disease / PD (mouse)
        - C26_plasma                     (scalar)
        - C26_CNS                        (scalar)
        - InflammCNS                     (scalar)
        - AxonInjury                     (scalar)

      [D] Safety & Immunogenicity
        - Scytokine                      (scalar)
        - S_ALT                          (scalar)
        - S_AST                          (scalar)
        - ADA_k for each platform k      (K-component vector)
    """

    K: int
    C: int
    N_cell: int

    def __post_init__(self) -> None:
        idx = 0

        # Block A: C_kc
        self.slice_C = slice(idx, idx + self.K * self.C)
        idx += self.K * self.C

        # Block A: U_kc_cell
        self.slice_U = slice(idx, idx + self.K * self.C * self.N_cell)
        idx += self.K * self.C * self.N_cell

        # Block B: m, P_tot, P_peri (per cell)
        self.slice_m = slice(idx, idx + self.N_cell)
        idx += self.N_cell

        self.slice_P_tot = slice(idx, idx + self.N_cell)
        idx += self.N_cell

        self.slice_P_peri = slice(idx, idx + self.N_cell)
        idx += self.N_cell

        # Block C: mouse VLCFA + inflammation + axonal injury
        self.idx_C26_plasma = idx
        idx += 1

        self.idx_C26_CNS = idx
        idx += 1

        self.idx_InflammCNS = idx
        idx += 1

        self.idx_AxonInjury = idx
        idx += 1

        # Block D: safety states
        self.idx_S_cytokine = idx
        idx += 1

        self.idx_S_ALT = idx
        idx += 1

        self.idx_S_AST = idx
        idx += 1

        self.slice_ADA = slice(idx, idx + self.K)
        idx += self.K

        self.n_states: int = idx

    # Small helper methods for unpacking

    def unpack_A(self, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Return (C_kc, U_kc_cell)."""
        C_flat = y[self.slice_C]
        U_flat = y[self.slice_U]
        C = C_flat.reshape(self.K, self.C)
        U = U_flat.reshape(self.K, self.C, self.N_cell)
        return C, U

    def unpack_B(self, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (m, P_tot, P_peri) per cell."""
        m = y[self.slice_m]
        P_tot = y[self.slice_P_tot]
        P_peri = y[self.slice_P_peri]
        return m, P_tot, P_peri

    def unpack_C_mouse(self, y: np.ndarray) -> Tuple[float, float, float, float]:
        C26_plasma = float(y[self.idx_C26_plasma])
        C26_CNS = float(y[self.idx_C26_CNS])
        Inflamm = float(y[self.idx_InflammCNS])
        Axon = float(y[self.idx_AxonInjury])
        return C26_plasma, C26_CNS, Inflamm, Axon

    def unpack_D(self, y: np.ndarray) -> Tuple[float, float, float, np.ndarray]:
        S_cyt = float(y[self.idx_S_cytokine])
        S_ALT = float(y[self.idx_S_ALT])
        S_AST = float(y[self.idx_S_AST])
        ADA = y[self.slice_ADA]
        return S_cyt, S_ALT, S_AST, ADA


# ---------------------------------------------------------------------------
#  Block A: Delivery & Exposure parameters
# ---------------------------------------------------------------------------

@dataclass
class BlockAParams:
    """
    Parameters for PBPK / delivery (Block A).

    Shapes follow Stage-2:

      - V_c:           (C,)         volumes
      - Q:             (C, C)       intercompartmental flows (row→col or col→row, as chosen consistently)
      - CL_c:          (C,)         clearance rates
      - A_BBB:         scalar       BBB area
      - P_baseline:    scalar       baseline P_BBB
      - M_FUS:         scalar       multiplier during FUS window
      - t_FUS:         scalar
      - dt_FUS:        scalar
      - Vmax_uptake:   (K, C, N_cell)   Vmax,k,c,cell
      - Km_uptake:     (K, C, N_cell)   Km,k,c,cell
      - theta_surface: (K, C, N_cell)   ϕ_surface(θ_surface, ligands)
    """

    V_c: np.ndarray          # (C,)
    Q: np.ndarray            # (C, C)
    CL_c: np.ndarray         # (C,)

    A_BBB: float
    P_baseline: float
    M_FUS: float
    t_FUS: float
    dt_FUS: float

    Vmax_uptake: np.ndarray  # (K, C, N_cell)
    Km_uptake: np.ndarray    # (K, C, N_cell)
    theta_surface: np.ndarray  # (K, C, N_cell)


# ---------------------------------------------------------------------------
#  Block B: Expression & peroxisomal targeting
# ---------------------------------------------------------------------------

@dataclass
class BlockBParams:
    """
    Parameters for Block B.

    ηesc and α_k,c,cell link uptake to mRNA input (eq. 7):

        k_in(t) = ηesc * sum_{k,c} α_{k,c,cell} * dU_{k,c,cell}/dt

    The rest follows Stage-2 eqs. 8–14.
    """

    eta_esc: np.ndarray              # (K,) endosomal escape efficiency per platform
    alpha_kc_cell: np.ndarray        # (K, C, N_cell) allocation factors (sum over cell ~ 1 per (k,c))

    kdeg_m: float                    # mRNA degradation
    ktr: float                       # translation m → P_tot
    kdeg_p: float                    # protein degradation (all pools)

    kins_max: float                  # maximum peroxisomal insertion rate
    K_ins: float                     # saturation vs P_tot
    C_PEX: float                     # PEX capacity
    kturnover_peri: float            # peroxisomal ALDP turnover

    Vbeta_WT: float                  # wild-type capacity
    PALDP_WT: float                  # reference peroxisomal ALDP
    K_beta: float                    # Hill half-saturation
    gamma_beta: float                # Hill exponent

    cell_weights: np.ndarray         # (N_cell,) weights to aggregate cell-level Vβ into global Vβ


# ---------------------------------------------------------------------------
#  Block C (mouse): Disease / PD
# ---------------------------------------------------------------------------

@dataclass
class BlockCParamsMouse:
    """
    Mouse branch parameters for Block C.

    NOTE on Vmax_C26:
      - Stage-2 lists Vmax,C26 as a parameter but the mouse VLCFA ODEs use Vβ(t) * h(C26)
        with h(C26) = C26 / (Km + C26), i.e. Vβ(t) already carries the “capacity” scaling.
      - To keep the I/O contract compatible, we keep Vmax_C26 as a parameter but DO NOT use it in
        the ODEs, mirroring the Stage-2 Formula Pack exactly.
    """

    # VLCFA production rates can be given as callables of time (and optionally phenotype),
    # but for this core we assume simple time-only callables.
    Pprod_plasma: Callable[[float], float]
    Pprod_CNS: Callable[[float], float]

    Km_C26: float
    Vmax_C26: float  # intentionally unused in ODEs; see note above

    kclear_plasma: float
    kclear_CNS: float
    k_cross: float

    kappa_infl: float
    kappa_resolve: float

    # ϕ_infl and Faxon are provided as callables to keep PD flexible:
    #   ϕ_infl(C26_CNS) ≥ 0; typically Hill-like
    #   Faxon(C26_CNS, Inflamm, AxonInjury) = d(AxonInjury)/dt
    phi_infl: Callable[[float], float]
    Faxon: Callable[[float, float, float], float]

    # Endpoint mapping parameters (mouse)
    NfL_baseline: float
    alpha_NfL: float

    FA_baseline: float
    alpha_FA: float

    Rotarod_baseline: float
    alpha_Rotarod: float


# ---------------------------------------------------------------------------
#  Block D: Safety & Immunogenicity
# ---------------------------------------------------------------------------

@dataclass
class BlockDParams:
    """
    Parameters and driver functions for Block D (safety).

    Drivers:

      E_cyt(t)  = f_cyt(t, C)              where C has shape (K, C)
      E_hep(t)  = f_hep(t, C, C26_pl, C26_CNS)
      E_ADA_k   = f_ADA_k(t, C)  (one function per platform)
    """

    k_cyt_ind: float
    k_cyt_decay: float

    k_ALT_ind: float
    k_ALT_decay: float

    k_AST_ind: float
    k_AST_decay: float

    k_ADA_ind: np.ndarray          # (K,)
    k_ADA_decay: np.ndarray        # (K,)

    # Safety exposure drivers
    f_cyt: Callable[[float, np.ndarray], float]
    f_hep: Callable[[float, np.ndarray, float, float], float]
    f_ADA: Sequence[Callable[[float, np.ndarray], float]]  # length K

    # Safe bounds (used later in Block G / decision layer)
    Esafe_cyt: float
    Esafe_ALT: float
    Esafe_AST: float


# ---------------------------------------------------------------------------
#  Block E: CMC / batch variability
# ---------------------------------------------------------------------------

@dataclass
class BlockEParams:
    """
    CMC block: algebraic maps from batch-level CMC metrics to effective factors.

    For each platform k and batch b, we observe a vector Q_{k,b} and compute:

      φ_k,b   = Φ_k(Q_{k,b})        → dose scaling
      ψ_k,b   = Ψ_k(Q_{k,b})        → effective η_esc factor (NOT yet wired in core)
      θ_mRNA,b= Θ_mRNA(Q_{k,b})     → effective kdeg_m etc. (NOT yet wired in core)

    Only φ_k,b is operationally used in this core_v3; ψ and θ are kept for future extensions.
    """

    # φ_k(Q), ψ_k(Q), θ_mRNA_k(Q) — indexed by platform index k
    phi_fns: Sequence[Callable[[np.ndarray], float]]
    psi_fns: Sequence[Callable[[np.ndarray], float]]
    theta_mrna_fns: Sequence[Callable[[np.ndarray], float]]

    # CMC metrics Q_{k,b}
    cmc_metrics: Mapping[Tuple[int, int], np.ndarray]  # key: (platform_index, batch_index)

    def effective_scale_and_escape(
        self,
        platform: int,
        batch_index: int,
    ) -> Tuple[float, float, Optional[float]]:
        """
        Return (φ_k,b, ψ_k,b, θ_mRNA,b) for a given platform and batch.

        If no CMC metrics are available for (k,b), returns (1.0, 1.0, None).
        """
        Q = self.cmc_metrics.get((platform, batch_index))
        if Q is None:
            return 1.0, 1.0, None

        phi = self.phi_fns[platform](Q) if self.phi_fns else 1.0
        psi = self.psi_fns[platform](Q) if self.psi_fns else 1.0
        theta = self.theta_mrna_fns[platform](Q) if self.theta_mrna_fns else None
        return float(phi), float(psi), (float(theta) if theta is not None else None)


# ---------------------------------------------------------------------------
#  Dose specification
# ---------------------------------------------------------------------------

@dataclass
class DoseEvent:
    """
    A single bolus dose event.

      time         : administration time t_j
      platform     : index k in config.platforms
      compartment  : index c_r in config.compartments (dosing compartment for this route)
      amount       : nominal D_{k,r}(t_j) per body weight
      batch_index  : batch b used to look up CMC features Q_{k,b}
    """

    time: float
    platform: int
    compartment: int
    amount: float
    batch_index: int


# ---------------------------------------------------------------------------
#  Core model
# ---------------------------------------------------------------------------

class ALDModel:
    """
    Core ALD dual-platform mRNA therapy model (mouse branch only).

    Blocks implemented:

      - A: PBPK delivery (C_kc, U_kc,cell)
      - B: Expression (m, P_tot, P_peri, Vβ algebraic)
      - C(mouse): VLCFA, CNS inflammation, axonal injury, NfL/FA/Rotarod
      - D: Cytokines, ALT, AST, ADA_k
      - E: CMC → φ_k,b (dose scaling) [ψ_k,b, θ_mRNA,b present but not wired yet]

    Human-projected Block C (NfL_human, Loes, etc.) and the full scenario / decision layers
    (Blocks F & G) are intentionally out of scope in this core_v3 and are meant to be built
    on top.
    """

    def __init__(
        self,
        config: ModelConfig,
        blockA: BlockAParams,
        blockB: BlockBParams,
        blockC_mouse: BlockCParamsMouse,
        blockD: BlockDParams,
        blockE: Optional[BlockEParams] = None,
    ) -> None:
        self.config = config
        self.A = blockA
        self.B = blockB
        self.Cm = blockC_mouse
        self.D = blockD
        self.E = blockE

        self.layout = StateIndexLayout(
            K=config.K,
            C=config.C,
            N_cell=config.N_cell,
        )

        # Quick sanity checks on shapes
        assert self.A.V_c.shape == (self.config.C,)
        assert self.A.CL_c.shape == (self.config.C,)
        assert self.A.Q.shape == (self.config.C, self.config.C)

        assert self.A.Vmax_uptake.shape == (self.config.K, self.config.C, self.config.N_cell)
        assert self.A.Km_uptake.shape == (self.config.K, self.config.C, self.config.N_cell)
        assert self.A.theta_surface.shape == (self.config.K, self.config.C, self.config.N_cell)

        assert self.B.eta_esc.shape == (self.config.K,)
        assert self.B.alpha_kc_cell.shape == (self.config.K, self.config.C, self.config.N_cell)
        assert self.B.cell_weights.shape == (self.config.N_cell,)

        assert self.D.k_ADA_ind.shape == (self.config.K,)
        assert self.D.k_ADA_decay.shape == (self.config.K,)
        assert len(self.D.f_ADA) == self.config.K

    # ------------------------------------------------------------------
    #  Initial state
    # ------------------------------------------------------------------

    def make_initial_state(
        self,
        C0: Optional[np.ndarray] = None,
        U0: Optional[np.ndarray] = None,
        m0: Optional[np.ndarray] = None,
        P_tot0: Optional[np.ndarray] = None,
        P_peri0: Optional[np.ndarray] = None,
        C26_plasma0: float = 0.0,
        C26_CNS0: float = 0.0,
        Inflamm0: float = 0.0,
        AxonInjury0: float = 0.0,
        Scytokine0: float = 0.0,
        S_ALT0: float = 0.0,
        S_AST0: float = 0.0,
        ADA0: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Construct an initial state vector y0 matching the layout.

        Any None entries are replaced with sensible defaults (zeros or Stage-2 baselines).
        """

        L = self.layout
        y0 = np.zeros(L.n_states, dtype=float)

        # Block A: C_kc
        if C0 is None:
            C0 = np.zeros((self.config.K, self.config.C), dtype=float)
        assert C0.shape == (self.config.K, self.config.C)
        y0[L.slice_C] = C0.ravel()

        # Block A: U_kc_cell
        if U0 is None:
            U0 = np.zeros((self.config.K, self.config.C, self.config.N_cell), dtype=float)
        assert U0.shape == (self.config.K, self.config.C, self.config.N_cell)
        y0[L.slice_U] = U0.ravel()

        # Block B: m, P_tot, P_peri
        if m0 is None:
            m0 = np.zeros(self.config.N_cell, dtype=float)
        if P_tot0 is None:
            P_tot0 = np.zeros(self.config.N_cell, dtype=float)
        if P_peri0 is None:
            P_peri0 = np.zeros(self.config.N_cell, dtype=float)

        assert m0.shape == (self.config.N_cell,)
        assert P_tot0.shape == (self.config.N_cell,)
        assert P_peri0.shape == (self.config.N_cell,)

        y0[L.slice_m] = m0
        y0[L.slice_P_tot] = P_tot0
        y0[L.slice_P_peri] = P_peri0

        # Block C (mouse)
        y0[L.idx_C26_plasma] = float(C26_plasma0)
        y0[L.idx_C26_CNS] = float(C26_CNS0)
        y0[L.idx_InflammCNS] = float(Inflamm0)
        y0[L.idx_AxonInjury] = float(AxonInjury0)

        # Block D
        y0[L.idx_S_cytokine] = float(Scytokine0)
        y0[L.idx_S_ALT] = float(S_ALT0)
        y0[L.idx_S_AST] = float(S_AST0)

        if ADA0 is None:
            ADA0 = np.zeros(self.config.K, dtype=float)
        assert ADA0.shape == (self.config.K,)
        y0[L.slice_ADA] = ADA0

        return y0

    # ------------------------------------------------------------------
    #  Block A helpers
    # ------------------------------------------------------------------

    def _P_BBB(self, t: float) -> float:
        """Time-dependent BBB permeability P_BBB(t) with FUS window."""
        if self.A.dt_FUS <= 0.0:
            return self.A.P_baseline
        if self.A.t_FUS <= t <= self.A.t_FUS + self.A.dt_FUS:
            return self.A.P_baseline * self.A.M_FUS
        return self.A.P_baseline

    def _v_uptake(self, C: np.ndarray) -> np.ndarray:
        """
        Compute v_uptake[k,c,cell](t) according to eq. (5):

            v = Vmax * C / (Km + C) * ϕ_surface(θ_surface)

        Here ϕ_surface is pre-encoded as theta_surface ∈ (0,1].
        """
        # Ensure non-negative concentrations for numerical stability
        C_pos = np.maximum(C, 0.0)

        # Broadcast C[k,c] over cell axis
        C_bc = C_pos[..., None]  # (K, C, 1) → will broadcast over N_cell

        num = self.A.Vmax_uptake * C_bc
        denom = self.A.Km_uptake + C_bc + 1e-12
        v = num / denom * self.A.theta_surface
        return v  # (K, C, N_cell)

    # ------------------------------------------------------------------
    #  Block B helpers
    # ------------------------------------------------------------------

    def _compute_Vbeta_per_cell(self, P_peri: np.ndarray) -> np.ndarray:
        """
        Compute Vβ_cell(t) for each cell using eqs. (13)-(14):

            x = PALDP_peri / PALDP_WT
            g(x) = x^γ / (x^γ + Kβ^γ)
            Vβ_cell = Vβ_WT * g(x)
        """
        x = np.maximum(P_peri, 0.0) / (self.B.PALDP_WT + 1e-12)
        x_gamma = np.power(x, self.B.gamma_beta)
        denom = x_gamma + np.power(self.B.K_beta, self.B.gamma_beta)
        g = np.where(denom > 0.0, x_gamma / denom, 0.0)
        Vbeta_cell = self.B.Vbeta_WT * g
        return Vbeta_cell  # (N_cell,)

    def _compute_Vbeta_global(self, P_peri: np.ndarray) -> float:
        """Aggregate cell-level Vβ into a single global Vβ(t) used in Block C."""
        Vbeta_cell = self._compute_Vbeta_per_cell(P_peri)
        return float(np.dot(self.B.cell_weights, Vbeta_cell))

    # ------------------------------------------------------------------
    #  Block C (mouse) helpers
    # ------------------------------------------------------------------

    def _h_C26(self, C26: float) -> float:
        """
        Michaelis–Menten form h(C26) = C26 / (Km + C26).

        NOTE: Vmax_C26 is *not* used here, by design; Vβ(t) (imported from Block B) plays the
        role of a capacity-like prefactor, matching Stage-2 equations exactly.
        """
        C26_pos = max(C26, 0.0)
        return C26_pos / (self.Cm.Km_C26 + C26_pos + 1e-12)

    # ------------------------------------------------------------------
    #  Right-hand side: dy/dt
    # ------------------------------------------------------------------

    def _rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        L = self.layout

        # Unpack blocks
        C, U = L.unpack_A(y)                    # C_kc, U_kc_cell
        m, P_tot, P_peri = L.unpack_B(y)        # per cell
        C26_plasma, C26_CNS, Inflamm, Axon = L.unpack_C_mouse(y)
        S_cyt, S_ALT, S_AST, ADA = L.unpack_D(y)

        # Allocate derivative vector
        dydt = np.zeros_like(y)

        # --- Block A: PBPK / delivery ---------------------------------
        V_c = self.A.V_c
        Q = self.A.Q
        CL_c = self.A.CL_c

        v_uptake = self._v_uptake(C)  # (K, C, N_cell)

        # Carrier concentrations C_kc
        dC_dt = np.zeros_like(C)

        K, Cn = self.config.K, self.config.C

        # Inflow / outflow / clearance
        for c in range(Cn):
            # inflow: sum_{c' != c} Q_{c'→c} / V_c[c] * C[k,c']
            inflow = np.sum((Q[:, c][None, :] / V_c[c]) * C, axis=1)  # (K,)
            # outflow: sum_{c'' != c} Q_{c→c''} / V_c[c] * C[k,c]
            outflow = np.sum(Q[c, :]) / V_c[c] * C[:, c]              # (K,)
            # clearance
            clearance = CL_c[c] / V_c[c] * C[:, c]                   # (K,)
            # uptake
            uptake_term = np.sum(v_uptake[:, c, :], axis=1) / V_c[c]  # (K,)

            dC_dt[:, c] += inflow - outflow - clearance - uptake_term

        # BBB flux between plasma and CNS (eqs. 3–4)
        idx_pl = self.config.idx_plasma
        idx_cns = self.config.idx_CNS
        P_BBB_t = self._P_BBB(t)
        flux = P_BBB_t * self.A.A_BBB * (C[:, idx_pl] - C[:, idx_cns])  # (K,)

        dC_dt[:, idx_pl] -= flux / V_c[idx_pl]
        dC_dt[:, idx_cns] += flux / V_c[idx_cns]

        # Write back dC_dt into dydt
        dydt[L.slice_C] = dC_dt.ravel()

        # Cumulative uptake U_kc_cell: dU/dt = v_uptake
        dydt[L.slice_U] = v_uptake.ravel()

        # --- Block B: Expression & peroxisomal targeting --------------
        # mRNA input rate kin(t) = sum_kc η_esc[k] * α_{k,c,cell} * v_uptake[k,c,cell]
        # v_uptake is already dU/dt, so we can use it directly.

        # Broadcast η_esc[k] over (c, cell)
        eta = self.B.eta_esc[:, None, None]  # (K,1,1)
        kin_cell = np.sum(eta * self.B.alpha_kc_cell * v_uptake, axis=(0, 1))  # (N_cell,)

        dm_dt = kin_cell - self.B.kdeg_m * m
        dP_tot_dt = self.B.ktr * m - self.B.kdeg_p * P_tot

        # Peroxisomal insertion (eq. 10) and P_peri dynamics (eq. 11)
        P_tot_pos = np.maximum(P_tot, 0.0)
        kins = (
            self.B.kins_max
            * P_tot_pos / (self.B.K_ins + P_tot_pos + 1e-12)
            * self.B.C_PEX / (self.B.C_PEX + P_tot_pos + 1e-12)
        )
        dP_peri_dt = kins - self.B.kturnover_peri * P_peri

        # Write back Block B derivatives
        dydt[L.slice_m] = dm_dt
        dydt[L.slice_P_tot] = dP_tot_dt
        dydt[L.slice_P_peri] = dP_peri_dt

        # --- Block C: Disease / PD (mouse) -----------------------------
        # Compute global Vβ(t) from Block B
        Vbeta_global = self._compute_Vbeta_global(P_peri)

        # VLCFA dynamics (eqs. 16–17)
        h_plasma = self._h_C26(C26_plasma)
        h_CNS = self._h_C26(C26_CNS)

        dC26_plasma_dt = (
            self.Cm.Pprod_plasma(t)
            - Vbeta_global * h_plasma
            - self.Cm.kclear_plasma * C26_plasma
        )

        dC26_CNS_dt = (
            self.Cm.Pprod_CNS(t)
            - Vbeta_global * h_CNS
            - self.Cm.kclear_CNS * C26_CNS
            + self.Cm.k_cross * (C26_plasma - C26_CNS)
        )

        # Inflammation & axonal injury
        dInflamm_dt = self.Cm.kappa_infl * self.Cm.phi_infl(C26_CNS) - self.Cm.kappa_resolve * Inflamm
        dAxon_dt = self.Cm.Faxon(C26_CNS, Inflamm, Axon)

        dydt[L.idx_C26_plasma] = dC26_plasma_dt
        dydt[L.idx_C26_CNS] = dC26_CNS_dt
        dydt[L.idx_InflammCNS] = dInflamm_dt
        dydt[L.idx_AxonInjury] = dAxon_dt

        # --- Block D: Safety & Immunogenicity --------------------------
        # Exposure drivers
        E_cyt = self.D.f_cyt(t, C)
        E_hep = self.D.f_hep(t, C, C26_plasma, C26_CNS)

        E_ADA = np.zeros(self.config.K, dtype=float)
        for k in range(self.config.K):
            E_ADA[k] = self.D.f_ADA[k](t, C)

        dS_cyt_dt = self.D.k_cyt_ind * E_cyt - self.D.k_cyt_decay * S_cyt
        dS_ALT_dt = self.D.k_ALT_ind * E_hep - self.D.k_ALT_decay * S_ALT
        dS_AST_dt = self.D.k_AST_ind * E_hep - self.D.k_AST_decay * S_AST

        dADA_dt = self.D.k_ADA_ind * E_ADA - self.D.k_ADA_decay * ADA

        dydt[L.idx_S_cytokine] = dS_cyt_dt
        dydt[L.idx_S_ALT] = dS_ALT_dt
        dydt[L.idx_S_AST] = dS_AST_dt
        dydt[L.slice_ADA] = dADA_dt

        return dydt

    # ------------------------------------------------------------------
    #  Dosing (Block A + E)
    # ------------------------------------------------------------------

    def _apply_bolus_doses_at_time(
        self,
        t: float,
        y: np.ndarray,
        dose_events: Sequence[DoseEvent],
    ) -> np.ndarray:
        """
        Apply all bolus doses scheduled at time t (within a small tolerance).

        For each dose event d at time t:

          - Look up φ_k,b from Block E (if present).
          - Compute effective dose Deff = φ_k,b * D_nominal.
          - Add Deff / V_c[comp] to C_{k,comp}.
        """
        if not dose_events:
            return y

        L = self.layout
        C, U = L.unpack_A(y)  # we will update C and write back

        V_c = self.A.V_c
        tol = 1e-9 * max(1.0, abs(t))

        for d in dose_events:
            if abs(d.time - t) > tol:
                continue

            phi_kb, psi_kb, theta_mrna_b = (1.0, 1.0, None)
            if self.E is not None:
                phi_kb, psi_kb, theta_mrna_b = self.E.effective_scale_and_escape(
                    d.platform, d.batch_index
                )

            Deff = phi_kb * d.amount

            c = d.compartment
            k = d.platform

            # Add Dirac-like jump: C_{k,c} ← C_{k,c} + Deff / V_c[c]
            C[k, c] += Deff / (V_c[c] + 1e-12)

            # NOTE: As per current scope, ψ_k,b and θ_mRNA,b are NOT yet wired into
            # η_esc or kdeg_m. They are carried at the CMC layer for future extensions.

        # Write updated C back into y
        y_new = y.copy()
        y_new[L.slice_C] = C.ravel()
        # U remains unchanged (no direct jump from dose)
        return y_new

    # ------------------------------------------------------------------
    #  Simulation (piecewise integration with robust dose handling)
    # ------------------------------------------------------------------

    def simulate(
        self,
        t_span: Tuple[float, float],
        y0: np.ndarray,
        doses: Sequence[DoseEvent],
        *,
        rtol: float = 1e-6,
        atol: float = 1e-9,
        max_step: Optional[float] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Integrate the ODE system over [t0, t_end] with instantaneous bolus doses.

        Key design point vs. core_v2:
          - No exact `t_stop in dose_times` check.
          - All dose times are handled with a tolerance, preventing floating-point
            round-off from silently skipping a final-time dose.

        Algorithm:
          1. Sort dose events by time, group them by unique times in (t0, t_end].
          2. Integrate from segment to segment using solve_ivp.
          3. At each dose time, apply all doses at that time as a jump in C_kc.
        """

        t0, t_end = float(t_span[0]), float(t_span[1])
        assert y0.shape == (self.layout.n_states,)

        doses_sorted = sorted(doses, key=lambda d: d.time)
        tol_time = 1e-9 * max(1.0, abs(t_end - t0))

        # Unique dose times within [t0, t_end] (with tolerance)
        unique_times: List[float] = []
        for d in doses_sorted:
            if d.time < t0 - tol_time or d.time > t_end + tol_time:
                continue
            if not unique_times:
                unique_times.append(d.time)
            else:
                if all(abs(d.time - t_exist) > tol_time for t_exist in unique_times):
                    unique_times.append(d.time)
        unique_times.sort()

        # Helper: filter dose events at a given time
        def _doses_at_time(t_event: float) -> List[DoseEvent]:
            return [d for d in doses_sorted if abs(d.time - t_event) <= tol_time]

        t_current = t0
        y_current = y0.copy()

        ts: List[float] = [t_current]
        ys: List[np.ndarray] = [y_current.copy()]

        def _integrate_segment(t_start: float, t_stop: float, y_start: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            if t_stop <= t_start:
                return np.array([t_start]), y_start[None, :]
            sol = solve_ivp(
                fun=self._rhs,
                t_span=(t_start, t_stop),
                y0=y_start,
                rtol=rtol,
                atol=atol,
                max_step=max_step if max_step is not None else np.inf,
                dense_output=False,
                vectorized=False,
            )
            if not sol.success:
                raise RuntimeError(f"ODESolver failed between t={t_start} and t={t_stop}: {sol.message}")
            return sol.t, sol.y.T  # (n_times,), (n_times, n_states)

        # Step through each dose time
        for t_dose in unique_times:
            # Integrate up to t_dose (or up to t_end if t_dose is beyond it)
            t_segment_end = min(max(t_dose, t_current), t_end)
            if t_segment_end - t_current > tol_time:
                t_seg, y_seg = _integrate_segment(t_current, t_segment_end, y_current)
                # Append, skipping duplicate first point
                ts.extend(list(t_seg[1:]))
                ys.extend(list(y_seg[1:]))
                t_current = float(t_seg[-1])
                y_current = y_seg[-1].copy()

            if t_dose > t_end + tol_time:
                break

            # Apply bolus at t_dose if within [t0, t_end]
            if t_dose >= t0 - tol_time and t_dose <= t_end + tol_time:
                y_current = self._apply_bolus_doses_at_time(t_dose, y_current, doses_sorted)
                t_current = t_dose
                ts.append(t_current)
                ys.append(y_current.copy())

        # Final segment from last time to t_end
        if t_end - t_current > tol_time:
            t_seg, y_seg = _integrate_segment(t_current, t_end, y_current)
            ts.extend(list(t_seg[1:]))
            ys.extend(list(y_seg[1:]))

        t_array = np.asarray(ts, dtype=float)
        y_array = np.vstack(ys)
        return t_array, y_array

    # ------------------------------------------------------------------
    #  Derived outputs (mouse endpoints + key PD metrics)
    # ------------------------------------------------------------------

    def compute_derived_outputs(
        self,
        t: np.ndarray,
        Y: np.ndarray,
    ) -> Dict[str, float]:
        """
        Compute scalar derived outputs from a time-course solution (mouse branch).

        Returns (at final time):
          - DeltaC26_plasma_rel
          - DeltaC26_CNS_rel
          - NfL_mouse
          - FA_mouse
          - Rotarod_mouse
        """
        L = self.layout

        C26_plasma_traj = Y[:, L.idx_C26_plasma]
        C26_CNS_traj = Y[:, L.idx_C26_CNS]
        Axon_traj = Y[:, L.idx_AxonInjury]

        C26_plasma_baseline = C26_plasma_traj[0]
        C26_CNS_baseline = C26_CNS_traj[0]

        C26_plasma_end = C26_plasma_traj[-1]
        C26_CNS_end = C26_CNS_traj[-1]
        Axon_end = Axon_traj[-1]

        # Relative changes (∆C / C_baseline)
        def _rel_change(final: float, baseline: float) -> float:
            if abs(baseline) < 1e-12:
                return 0.0
            return (final - baseline) / baseline

        dC26_plasma_rel = _rel_change(C26_plasma_end, C26_plasma_baseline)
        dC26_CNS_rel = _rel_change(C26_CNS_end, C26_CNS_baseline)

        # Mouse endpoints (eqs. 21–23, without explicit noise terms)
        NfL_mouse = self.Cm.NfL_baseline + self.Cm.alpha_NfL * Axon_end
        FA_mouse = self.Cm.FA_baseline - self.Cm.alpha_FA * Axon_end
        Rotarod_mouse = self.Cm.Rotarod_baseline - self.Cm.alpha_Rotarod * Axon_end

        return {
            "DeltaC26_plasma_rel": float(dC26_plasma_rel),
            "DeltaC26_CNS_rel": float(dC26_CNS_rel),
            "NfL_mouse": float(NfL_mouse),
            "FA_mouse": float(FA_mouse),
            "Rotarod_mouse": float(Rotarod_mouse),
        }
