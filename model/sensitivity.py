#!/usr/bin/env python3
"""Analisis de sensibilidad global (Sobol) sobre el modelo cinetico.

Las conclusiones del modelo salen de parametros de literatura, no de
mediciones de este proyecto. Este analisis responde cuanto de cada conclusion
depende de esa eleccion:

  - indice S1 (primer orden): cuanta varianza explica un parametro por si solo
  - indice ST (total): cuanta explica incluyendo sus interacciones
  - ST >> S1 indica que el parametro actua sobre todo por interaccion

Se evaluan tres salidas, una por conclusion del modelo:

  1. t_signal   -- tiempo hasta senal visible del sensor (v3)
  2. margin     -- margen de MazF libre bajo el umbral letal con el plasmido
                   retenido; negativo significa que el switch dispara solo
  3. kill_delay -- tiempo hasta la muerte tras perder el plasmido

Cada parametro se muestrea en un rango que refleja su incertidumbre real:
las vidas medias y las tasas de sintesis estan bien acotadas en literatura,
mientras que los umbrales de deteccion y de letalidad son estimaciones
gruesas y se muestrean en rangos mas amplios.

Uso:  .venv/bin/python model/sensitivity.py [--fast]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import kinetics as kn  # noqa: E402
from model.kinetics import KillSwitchConfig, SensorConfig  # noqa: E402

FIG_DIR = Path(__file__).resolve().parent / "figures"

# (parametro, factor de incertidumbre, justificacion del rango)
# El rango es [valor/factor, valor*factor]. Factor 2 = medio orden de
# magnitud a cada lado; factor 4 = parametro mal acotado.
UNCERTAINTY = [
    ("k_txn_max", 3.0, "fuerza del promotor, no medida en este chasis"),
    ("leak", 3.0, "expresion basal, muy variable entre sistemas"),
    ("k_tln", 3.0, "tasa de traduccion por mRNA"),
    ("d_mrna", 2.0, "vida media de mRNA, bien acotada"),
    ("d_prot", 2.0, "dilucion por division, depende del medio"),
    ("t_mat_chromo", 2.0, "maduracion de cromoproteina"),
    ("thr_visual", 4.0, "umbral visual, estimacion gruesa"),
    ("k_mazF", 3.0, "sintesis de toxina"),
    ("k_mazE", 3.0, "sintesis de antitoxina"),
    ("k_on", 4.0, "asociacion MazE-MazF, mal acotada"),
    ("d_mazE", 2.0, "MazE degradada por ClpAP"),
    ("d_mazF_ssrA", 2.0, "MazF con tag ssrA"),
    ("tox_lethal", 4.0, "umbral letal, estimacion gruesa"),
]

# Eficiencias del RBS de la v3, fijas: son propiedad del diseno, no del
# entorno, y se quieren medir los parametros biologicos por separado.
EFF_E = kn.rbs_efficiency_from_sd(7)
EFF_F = kn.rbs_efficiency_from_sd(6)


def problem() -> dict:
    names = [p for p, _, _ in UNCERTAINTY]
    bounds = []
    for name, factor, _ in UNCERTAINTY:
        v = kn.p(name)
        bounds.append([v / factor, v * factor])
    return {"num_vars": len(names), "names": names, "bounds": bounds}


def evaluate(row: np.ndarray, names: list[str]) -> tuple[float, float, float]:
    """Corre el modelo con un juego de parametros y devuelve las tres salidas."""
    saved = {n: kn.PARAMS[n] for n in names}
    for n, v in zip(names, row):
        u, src = kn.PARAMS[n][1], kn.PARAMS[n][2]
        kn.PARAMS[n] = (float(v), u, src)
    try:
        # 1. tiempo hasta senal visible del sensor v3
        cfg = SensorConfig(reporter="chromoprotein", readout="visual",
                           rbs_efficiency=EFF_F, inducer=10.0)
        t_sig = kn.time_to_signal(cfg, t_end=21600.0)
        t_sig = 360.0 if t_sig is None else t_sig   # censurado a 6 h

        # 2. margen bajo el umbral letal con el plasmido retenido
        ks = KillSwitchConfig(rbs_mazE=EFF_E, rbs_mazF=EFF_F, ssrA=True)
        f_ss = kn.simulate_killswitch(ks, t_end=21600.0, dense=False)[1][1, -1]
        margin = kn.p("tox_lethal") - f_ss

        # 3. retardo hasta la muerte tras perder el plasmido
        lost = KillSwitchConfig(rbs_mazE=EFF_E, rbs_mazF=EFF_F, ssrA=True,
                                plasmid_lost_at=3600.0)
        t_kill = kn.time_to_death(lost, t_end=21600.0, dense=False)
        t_kill = 360.0 if t_kill is None else t_kill
    finally:
        for n, val in saved.items():
            kn.PARAMS[n] = val
    return t_sig, margin, t_kill


OUTPUTS = [
    ("t_signal", "tiempo hasta senal visible (min)"),
    ("margin", "margen bajo el umbral letal (nM)"),
    ("kill_delay", "retardo hasta la muerte (min)"),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="menos muestras; para iterar, no para reportar")
    args = ap.parse_args()

    from SALib.analyze import sobol as sobol_analyze
    from SALib.sample import sobol as sobol_sample

    prob = problem()
    n = 128 if args.fast else 1024
    X = sobol_sample.sample(prob, n, calc_second_order=False)

    print("=" * 70)
    print("Analisis de sensibilidad global (Sobol)")
    print("=" * 70)
    print(f"{len(prob['names'])} parametros, {len(X)} evaluaciones del modelo")
    print("Rangos: [valor/factor, valor*factor], factor por incertidumbre")
    print()

    Y = np.empty((len(X), 3))
    for i, row in enumerate(X):
        Y[i] = evaluate(row, prob["names"])
        if (i + 1) % max(1, len(X) // 10) == 0:
            print(f"  {i + 1:>6}/{len(X)}")

    results = {}
    for j, (key, label) in enumerate(OUTPUTS):
        si = sobol_analyze.analyze(prob, Y[:, j], calc_second_order=False,
                                   print_to_console=False)
        results[key] = si
        section = f"{key} — {label}"
        print()
        print(section)
        print("-" * len(section))
        print(f"  media {Y[:, j].mean():.1f}   "
              f"p5 {np.percentile(Y[:, j], 5):.1f}   "
              f"p95 {np.percentile(Y[:, j], 95):.1f}")
        print()
        print(f"  {'parametro':<16}{'S1':>8}{'ST':>8}   interpretacion")
        order = np.argsort(-si["ST"])
        for idx in order:
            s1, st = si["S1"][idx], si["ST"][idx]
            if st < 0.02:
                continue
            note = ""
            if st > 0.3:
                note = "dominante"
            elif st > 0.1:
                note = "relevante"
            if st > 2.5 * max(s1, 0.01):
                note += " (via interaccion)" if note else "solo por interaccion"
            print(f"  {prob['names'][idx]:<16}{s1:>8.3f}{st:>8.3f}   {note}")

    make_figure(prob, results, Y)
    verdicts(Y)
    return 0


def make_figure(prob: dict, results: dict, Y: np.ndarray) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})

    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    for ax, (key, label) in zip(axes, OUTPUTS):
        si = results[key]
        order = np.argsort(si["ST"])
        keep = [i for i in order if si["ST"][i] >= 0.02]
        names = [prob["names"][i] for i in keep]
        st = [si["ST"][i] for i in keep]
        s1 = [max(si["S1"][i], 0) for i in keep]
        y = np.arange(len(keep))
        ax.barh(y, st, color="#c0392b", alpha=0.45, label="ST (total)")
        ax.barh(y, s1, color="#2471a3", height=0.55, label="S1 (primer orden)")
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=7.5)
        ax.set_xlabel("indice de Sobol")
        ax.set_title(label, fontsize=8.5)
        ax.set_xlim(0, 1)
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Que parametros gobiernan cada conclusion", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_sensitivity.png")
    plt.close(fig)

    # --- la tension entre las dos condiciones del kill switch ---
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    safe, kills = Y[:, 1] > 0, Y[:, 2] < 360
    groups = [
        (safe & kills, "#1d6b2b", "cumple ambas"),
        (safe & ~kills, "#e67e22", "seguro, no contiene"),
        (~safe & kills, "#c0392b", "contiene, dispara solo"),
        (~safe & ~kills, "#7f8c8d", "ninguna"),
    ]
    for mask, color, label in groups:
        if mask.sum():
            ax.scatter(Y[mask, 1], np.clip(Y[mask, 2], 0, 380), s=5,
                       alpha=0.35, color=color,
                       label=f"{label} ({mask.mean():.0%})")
    ax.axvline(0, c="k", lw=1, ls="--")
    ax.axhline(360, c="k", lw=1, ls=":")
    ax.set_ylim(-15, 420)
    ax.annotate("nunca alcanza nivel letal\n(censurado a 6 h)",
                xy=(ax.get_xlim()[1] * 0.72, 360),
                xytext=(ax.get_xlim()[1] * 0.50, 405),
                fontsize=7.5, color="#555", ha="center",
                arrowprops=dict(arrowstyle="->", color="#888", lw=1))
    ax.set_xlabel("margen bajo el umbral letal en reposo (nM)")
    ax.set_ylabel("retardo hasta la muerte tras el escape (min)")
    ax.set_title("Margen amplio en reposo = la toxina nunca llega a matar")
    ax.legend(fontsize=7, loc="upper left", markerscale=2.5)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "06_killswitch_tradeoff.png")
    plt.close(fig)
    print()
    print(f"  figuras escritas en {FIG_DIR.name}/")


def verdicts(Y: np.ndarray) -> None:
    """Robustez de las conclusiones del modelo bajo incertidumbre."""
    print()
    print("Robustez de las conclusiones")
    print("-" * 28)

    slow = Y[:, 0] > 30
    safe = Y[:, 1] > 0        # no dispara solo con el plasmido retenido
    kills = Y[:, 2] < 360     # mata tras perder el plasmido

    print(f"  1. el sensor tarda mas de 30 min:         {slow.mean():>6.1%}")
    print(f"  2. el kill switch no dispara solo:        {safe.mean():>6.1%}")
    print(f"  3. mata al perder el plasmido:            {kills.mean():>6.1%}")
    print()
    print("  Cerca del 100% = estructural, se sostiene en todo el rango de")
    print("  incertidumbre. Cerca del 50% = depende de los parametros y hay")
    print("  que medirlo.")

    print()
    print("  Las dos condiciones del kill switch estan en tension:")
    print(f"    cumple ambas:                           {(safe & kills).mean():>6.1%}")
    print(f"    seguro pero no contiene:                {(safe & ~kills).mean():>6.1%}")
    print(f"    contiene pero dispara solo:             {(~safe & kills).mean():>6.1%}")
    print(f"    ni una ni otra:                         {(~safe & ~kills).mean():>6.1%}")
    print()
    margin_ok = Y[safe & kills, 1]
    if margin_ok.size:
        print(f"    entre las que cumplen ambas, margen mediano"
              f" {np.median(margin_ok):.0f} nM")

    # El fallo dominante no es el disparo espurio sino la falta de muerte:
    # con margen alto la toxina nunca alcanza concentracion letal.
    m_fail = Y[safe & ~kills, 1]
    m_ok = Y[safe & kills, 1]
    print()
    print("    El modo de fallo dominante es NO matar, no disparar solo:")
    if m_fail.size and m_ok.size:
        print(f"      margen mediano cuando SI mata:   {np.median(m_ok):>6.0f} nM")
        print(f"      margen mediano cuando NO mata:   {np.median(m_fail):>6.0f} nM")
    print("    Un margen de seguridad amplio en reposo es exactamente lo que")
    print("    impide que la toxina llegue a niveles letales tras el escape:")
    print("    con promotores constitutivos, seguridad y contencion compiten")
    print("    por el mismo parametro.")


if __name__ == "__main__":
    raise SystemExit(main())
