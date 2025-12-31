import argparse
import json
from pathlib import Path

import numpy as np
import yaml

from ald_model.core import (
    ModelConfig,
    BlockAParams,
    BlockBParams,
    BlockCParamsMouse,
    BlockDParams,
    BlockEParams,
    ALDModel,
    DoseEvent,
)


def build_block_params_from_yaml(cfg: dict):
    # --- ModelConfig ---
    mc = cfg["model_config"]
    config = ModelConfig(
        platforms=tuple(mc["platforms"]),
        compartments=tuple(mc["compartments"]),
        cell_types=tuple(mc["cell_types"]),
    )

    # --- Block A ---
    ba = cfg["blockA"]
    blockA = BlockAParams(
        V_c=np.asarray(ba["V_c"], float),
        Q=np.asarray(ba["Q"], float),
        CL_c=np.asarray(ba["CL_c"], float),
        A_BBB=float(ba["A_BBB"]),
        P_baseline=float(ba["P_baseline"]),
        M_FUS=float(ba["M_FUS"]),
        t_FUS=float(ba["t_FUS"]),
        dt_FUS=float(ba["dt_FUS"]),
        Vmax_uptake=np.asarray(ba["Vmax_uptake"], float),
        Km_uptake=np.asarray(ba["Km_uptake"], float),
        theta_surface=np.asarray(ba["theta_surface"], float),
    )

    # --- Block B ---
    bb = cfg["blockB"]
    blockB = BlockBParams(
        eta_esc=np.asarray(bb["eta_esc"], float),
        alpha_kc_cell=np.asarray(bb["alpha_kc_cell"], float),
        kdeg_m=float(bb["kdeg_m"]),
        ktr=float(bb["ktr"]),
        kdeg_p=float(bb["kdeg_p"]),
        kins_max=float(bb["kins_max"]),
        K_ins=float(bb["K_ins"]),
        C_PEX=float(bb["C_PEX"]),
        kturnover_peri=float(bb["kturnover_peri"]),
        Vbeta_WT=float(bb["Vbeta_WT"]),
        PALDP_WT=float(bb["PALDP_WT"]),
        K_beta=float(bb["K_beta"]),
        gamma_beta=float(bb["gamma_beta"]),
        cell_weights=np.asarray(bb["cell_weights"], float),
    )

    # --- Block C (mouse) ---
    bc = cfg["blockC_mouse"]

    def Pprod_plasma(t: float) -> float:
        return float(bc["Pprod_plasma_const"])

    def Pprod_CNS(t: float) -> float:
        return float(bc["Pprod_CNS_const"])

    def phi_infl(C26_CNS: float) -> float:
        C_half = float(bc["infl_C_half"])
        gamma = float(bc["infl_gamma"])
        C_pos = max(C26_CNS, 0.0)
        num = C_pos**gamma
        den = C_half**gamma + num + 1e-12
        return num / den

    def Faxon(C26_CNS: float, Inflamm: float, AxonInjury: float) -> float:
        # Example mapping; in a production version you may parameterize this.
        return 0.01 * Inflamm - 0.001 * AxonInjury

    blockC_mouse = BlockCParamsMouse(
        Pprod_plasma=Pprod_plasma,
        Pprod_CNS=Pprod_CNS,
        Km_C26=float(bc["Km_C26"]),
        Vmax_C26=float(bc["Vmax_C26"]),
        kclear_plasma=float(bc["kclear_plasma"]),
        kclear_CNS=float(bc["kclear_CNS"]),
        k_cross=float(bc["k_cross"]),
        kappa_infl=float(bc["kappa_infl"]),
        kappa_resolve=float(bc["kappa_resolve"]),
        phi_infl=phi_infl,
        Faxon=Faxon,
        NfL_baseline=float(bc["NfL_baseline"]),
        alpha_NfL=float(bc["alpha_NfL"]),
        FA_baseline=float(bc["FA_baseline"]),
        alpha_FA=float(bc["alpha_FA"]),
        Rotarod_baseline=float(bc["Rotarod_baseline"]),
        alpha_Rotarod=float(bc["alpha_Rotarod"]),
    )

    # --- Block D ---
    bd = cfg["blockD"]

    def f_cyt(t: float, C: np.ndarray) -> float:
        # Example driver: total carrier in plasma
        return float(np.sum(C[:, config.compartment_index["plasma"]]))

    def f_hep(t: float, C: np.ndarray, C26_pl: float, C26_CNS: float) -> float:
        # Example driver: liver conc of primary platform
        return float(C[0, config.compartment_index["liver"]])

    f_ADA = []
    for k in range(len(config.platforms)):
        # capture k correctly
        def _f_ADA(t: float, C: np.ndarray, kk=k) -> float:
            return float(np.sum(C[kk, :]))

        f_ADA.append(_f_ADA)

    blockD = BlockDParams(
        k_cyt_ind=float(bd["k_cyt_ind"]),
        k_cyt_decay=float(bd["k_cyt_decay"]),
        k_ALT_ind=float(bd["k_ALT_ind"]),
        k_ALT_decay=float(bd["k_ALT_decay"]),
        k_AST_ind=float(bd["k_AST_ind"]),
        k_AST_decay=float(bd["k_AST_decay"]),
        k_ADA_ind=np.asarray(bd["k_ADA_ind"], float),
        k_ADA_decay=np.asarray(bd["k_ADA_decay"], float),
        f_cyt=f_cyt,
        f_hep=f_hep,
        f_ADA=f_ADA,
        Esafe_cyt=float(bd["Esafe_cyt"]),
        Esafe_ALT=float(bd["Esafe_ALT"]),
        Esafe_AST=float(bd["Esafe_AST"]),
    )

    # --- Block E (CMC) ---
    be = cfg.get("blockE", {})
    cmc_raw = be.get("cmc_metrics", {})
    cmc_metrics = {}
    for key, vec in cmc_raw.items():
        k_str, b_str = key.split("_")
        cmc_metrics[(int(k_str), int(b_str))] = np.asarray(vec, float)

    def default_phi(Q: np.ndarray) -> float:
        return 1.0

    def default_psi(Q: np.ndarray) -> float:
        return 1.0

    def default_theta(Q: np.ndarray) -> float:
        return 1.0

    phi_fns = [default_phi for _ in config.platforms]
    psi_fns = [default_psi for _ in config.platforms]
    theta_mrna_fns = [default_theta for _ in config.platforms]

    blockE = BlockEParams(
        phi_fns=phi_fns,
        psi_fns=psi_fns,
        theta_mrna_fns=theta_mrna_fns,
        cmc_metrics=cmc_metrics,
    )

    return config, blockA, blockB, blockC_mouse, blockD, blockE


def build_dose_events_from_yaml(dose_cfg: dict, config: ModelConfig):
    events = []
    for d in dose_cfg.get("doses", []):
        t = float(d["time"])
        plat_name = d["platform"]
        comp_name = d["compartment"]
        amount = float(d["amount"])
        batch_index = int(d["batch_index"])

        k = config.platform_index[plat_name]
        c = config.compartment_index[comp_name]

        events.append(
            DoseEvent(
                time=t,
                platform=k,
                compartment=c,
                amount=amount,
                batch_index=batch_index,
            )
        )

    t0 = float(dose_cfg["t0"])
    tend = float(dose_cfg["tend"])
    return (t0, tend), events


def build_initial_state_from_yaml(init_cfg: dict, model: ALDModel):
    C0 = np.asarray(init_cfg["C0"], float)
    U0 = np.asarray(init_cfg["U0"], float)
    m0 = np.asarray(init_cfg["m0"], float)
    P_tot0 = np.asarray(init_cfg["P_tot0"], float)
    P_peri0 = np.asarray(init_cfg["P_peri0"], float)

    y0 = model.make_initial_state(
        C0=C0,
        U0=U0,
        m0=m0,
        P_tot0=P_tot0,
        P_peri0=P_peri0,
        C26_plasma0=float(init_cfg["C26_plasma0"]),
        C26_CNS0=float(init_cfg["C26_CNS0"]),
        Inflamm0=float(init_cfg["Inflamm0"]),
        AxonInjury0=float(init_cfg["AxonInjury0"]),
        Scytokine0=float(init_cfg["Scytokine0"]),
        S_ALT0=float(init_cfg["S_ALT0"]),
        S_AST0=float(init_cfg["S_AST0"]),
        ADA0=np.asarray(init_cfg["ADA0"], float),
    )
    return y0


def main():
    parser = argparse.ArgumentParser(description="Run ALD Stage-3 mouse-core example.")
    parser.add_argument("--params", type=str, required=True, help="Path to params YAML")
    parser.add_argument("--initial", type=str, required=True, help="Path to initial state YAML")
    parser.add_argument("--doses", type=str, required=True, help="Path to doses YAML")
    parser.add_argument("--outdir", type=str, default="outputs", help="Output directory")
    args = parser.parse_args()

    params_path = Path(args.params)
    initial_path = Path(args.initial)
    doses_path = Path(args.doses)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with params_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    config, blockA, blockB, blockC_mouse, blockD, blockE = build_block_params_from_yaml(cfg)
    model = ALDModel(
        config=config,
        blockA=blockA,
        blockB=blockB,
        blockC_mouse=blockC_mouse,
        blockD=blockD,
        blockE=blockE,
    )

    with initial_path.open("r", encoding="utf-8") as f:
        init_cfg = yaml.safe_load(f)
    y0 = build_initial_state_from_yaml(init_cfg, model)

    with doses_path.open("r", encoding="utf-8") as f:
        dose_cfg = yaml.safe_load(f)
    (t0, tend), dose_events = build_dose_events_from_yaml(dose_cfg, config)

    t, Y = model.simulate((t0, tend), y0, dose_events)
    derived = model.compute_derived_outputs(t, Y)

    # Save outputs
    np.savetxt(outdir / "sample_run_times.csv", t, delimiter=",")
    np.save(outdir / "sample_run_states.npy", Y)
    with (outdir / "sample_run_derived.json").open("w", encoding="utf-8") as f:
        json.dump(derived, f, indent=2)

    print("Simulation finished.")
    print("Derived outputs:")
    print(json.dumps(derived, indent=2))


if __name__ == "__main__":
    main()
