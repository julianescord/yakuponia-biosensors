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


CACHE_DIR = Path(__file__).resolve().parent / ".cache"


def _cache_key(prob: dict, n: int) -> Path:
    """Identifica un barrido por sus rangos y su tamano de muestra.

    Cambiar un rango de UNCERTAINTY o el numero de muestras invalida la cache;
    cambiar solo el estilo de las figuras no.
    """
    import hashlib
    import json
    blob = json.dumps({"names": prob["names"], "bounds": prob["bounds"],
                       "n": n}, sort_keys=True).encode()
    return CACHE_DIR / f"sobol_{hashlib.sha256(blob).hexdigest()[:12]}.npz"


def sample_model(prob: dict, n: int, *, refresh: bool = False):
    """Evalua el modelo en la muestra de Sobol, con cache en disco.

    El barrido cuesta minutos; las figuras se reajustan en segundos. Sin cache,
    cada retoque visual pagaria el barrido entero de nuevo.
    """
    from SALib.sample import sobol as sobol_sample

    X = sobol_sample.sample(prob, n, calc_second_order=False)
    path = _cache_key(prob, n)

    if path.exists() and not refresh:
        data = np.load(path)
        if data["X"].shape == X.shape:
            print(f"  muestras reutilizadas de {path.name}"
                  f"  ({len(X)} evaluaciones)")
            return X, data["Y"]

    Y = np.empty((len(X), 3))
    for i, row in enumerate(X):
        Y[i] = evaluate(row, prob["names"])
        if (i + 1) % max(1, len(X) // 10) == 0:
            print(f"  {i + 1:>6}/{len(X)}")

    CACHE_DIR.mkdir(exist_ok=True)
    np.savez_compressed(path, X=X, Y=Y)
    print(f"  muestras guardadas en {path.name}")
    return X, Y


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="menos muestras; para iterar, no para reportar")
    ap.add_argument("--refresh", action="store_true",
                    help="ignora la cache y recorre el modelo")
    ap.add_argument("--figures-only", action="store_true",
                    help="solo redibuja desde la cache; falla si no existe")
    args = ap.parse_args()

    from SALib.analyze import sobol as sobol_analyze

    prob = problem()
    n = 128 if args.fast else 1024

    print("=" * 70)
    print("Analisis de sensibilidad global (Sobol)")
    print("=" * 70)
    print(f"{len(prob['names'])} parametros, muestra base {n}")
    print("Rangos: [valor/factor, valor*factor], factor por incertidumbre")
    print()

    if args.figures_only and not _cache_key(prob, n).exists():
        print("  no hay cache para esta configuracion; corre sin --figures-only")
        return 1

    X, Y = sample_model(prob, n, refresh=args.refresh)

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
    from matplotlib.gridspec import GridSpec
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    from model import figstyle as fs

    fs.apply_style()

    # ----------------------------------------------------------------
    # Fig 5 — indices de Sobol por conclusion
    # ----------------------------------------------------------------
    # Un panel por conclusion. Los paneles comparten eje x (0-1) y el mismo
    # grosor de barra: si un panel retiene menos parametros, sobra espacio
    # abajo en vez de engordar sus barras, porque el grosor no debe leerse
    # como importancia.
    per_panel = []
    for key, _ in OUTPUTS:
        si = results[key]
        keep = [i for i in np.argsort(si["ST"]) if si["ST"][i] >= 0.02]
        per_panel.append(keep)
    n_max = max(len(k) for k in per_panel)

    fig = plt.figure(figsize=(fs.COL_DOUBLE, 2.55))
    gs = GridSpec(1, 3, figure=fig, wspace=0.62,
                  left=0.105, right=0.975, top=0.815, bottom=0.275)
    axes = [fig.add_subplot(gs[0, j]) for j in range(3)]

    for ax, keep, (key, label) in zip(axes, per_panel, OUTPUTS):
        si = results[key]
        # las barras van de menor ST abajo a mayor arriba; todas apoyadas en
        # la base para que los tres paneles compartan linea de arranque
        y = np.arange(len(keep))
        ax.barh(y, [si["ST"][i] for i in keep], height=0.62,
                color=fs.C["st"], label="ST (con interacciones)", zorder=2)
        ax.barh(y, [max(si["S1"][i], 0) for i in keep], height=0.30,
                color=fs.C["s1"], label="S1 (efecto propio)", zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([prob["names"][i] for i in keep])
        ax.tick_params(axis="y", length=0, pad=1.5)
        # mismo grosor de barra en los tres: el eje se fija al panel mayor
        # el tope se fija al panel mas poblado: iguala el grosor de barra
        # entre paneles, que si no se leeria como importancia
        ax.set_ylim(-0.75, n_max - 0.25)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xlabel("indice de Sobol")
        ax.set_title(label, fontsize=6.8, color=fs.C["ink"], pad=3.5)

    fs.label_panels(axes, "abc", dx_pt=-40.0)
    # leyenda bajo el panel a, fuera del area de barras
    axes[0].legend(loc="upper left", bbox_to_anchor=(-0.01, -0.235),
                   fontsize=6.0, handlelength=1.2, ncol=2, columnspacing=1.4)
    fig.suptitle("Que parametros gobiernan cada conclusion del modelo",
                 fontsize=8.0, y=0.965)
    fs.save(fig, "05_sensitivity")

    # ----------------------------------------------------------------
    # Fig 6 — la tension entre seguridad y contencion
    # ----------------------------------------------------------------
    # El 59% de las muestras nunca alcanza nivel letal y quedaba aplastado
    # sobre la linea de censura en y=360. Se separa en dos paneles: el (a)
    # muestra solo las que SI matan, con eje real; el (b) cuantifica cuantas
    # caen en cada cuadrante, que es la conclusion de verdad.
    safe, kills = Y[:, 1] > 0, Y[:, 2] < 360

    fig = plt.figure(figsize=(fs.COL_DOUBLE, 2.7))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.45, 1.0], wspace=0.30,
                  left=0.078, right=0.985, top=0.80, bottom=0.165)
    ax_s = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # (a) solo las muestras no censuradas: aqui el eje y significa algo
    live = kills
    for mask, color, lab in ((live & safe, fs.C["ok"], "cumple ambas"),
                             (live & ~safe, fs.C["fail"], "dispara sola")):
        if mask.sum():
            ax_s.scatter(Y[mask, 1], Y[mask, 2], s=2.6, alpha=0.30,
                         color=color, linewidths=0, label=lab, rasterized=True)
    ax_s.axvline(0, color=fs.C["ink"], lw=0.7, ls="--", zorder=4)
    ax_s.set_xlabel("margen bajo el umbral letal en reposo (nM)")
    ax_s.set_ylabel("retardo hasta la muerte\ntras el escape (min)")
    ax_s.set_title("Muestras que si alcanzan nivel letal", fontsize=6.8, pad=3.5)
    ax_s.legend(loc="upper left", markerscale=3.2, fontsize=6.0,
                borderaxespad=0.25)
    # nota dentro del area del panel: fuera se saldria de la pagina
    ax_s.text(0.015, 0.055, "margen negativo = el switch dispara\n"
              "con el plasmido retenido", transform=ax_s.transAxes,
              fontsize=5.6, color=fs.C["mid"], va="bottom", linespacing=1.4)

    # (b) los cuatro cuadrantes, incluido el 59% censurado que el panel (a)
    # no puede mostrar
    cats = [
        ("cumple ambas", safe & kills, fs.C["ok"]),
        ("seguro, no contiene", safe & ~kills, fs.C["warn"]),
        ("contiene, dispara sola", ~safe & kills, fs.C["fail"]),
        ("ninguna", ~safe & ~kills, fs.C["light"]),
    ]
    ypos = np.arange(len(cats))[::-1]
    for yp, (lab, mask, color) in zip(ypos, cats):
        frac = mask.mean()
        ax_b.barh(yp, frac * 100, height=0.60, color=color, zorder=2)
        # una categoria vacia no dibuja barra: se marca el cero para que no
        # se confunda con un dato que falta
        txt = f"{frac:.0%}" if frac >= 0.005 else "0% (ninguna muestra)"
        ax_b.text(frac * 100 + 1.6, yp, txt, va="center",
                  fontsize=6.5, color=fs.C["ink"])
    ax_b.set_yticks(ypos)
    ax_b.set_yticklabels([c[0] for c in cats], fontsize=6.3)
    ax_b.tick_params(axis="y", length=0, pad=1.5)
    ax_b.set_xlim(0, 100)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_xlabel("muestras (%)")
    ax_b.set_title("El fallo dominante es no matar", fontsize=6.8, pad=3.5)
    ax_b.set_ylim(-0.7, len(cats) - 0.3)

    fs.label_panels([ax_s], "a", dx_pt=-30.0)
    fs.label_panels([ax_b], "b", dx_pt=-74.0)
    fig.suptitle("Un margen amplio en reposo impide que la toxina llegue "
                 "a matar", fontsize=8.0, y=0.965)
    fs.save(fig, "06_killswitch_tradeoff")

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
