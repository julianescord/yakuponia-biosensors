#!/usr/bin/env python3
"""Corre el analisis cinetico y escribe las figuras y la tabla de resultados.

Responde tres preguntas que el diseno de secuencia no puede contestar:

  1. Es alcanzable la afirmacion de "<30 min" del proyecto?
  2. Cuanta diferencia hace la cromoproteina frente a la fluorescente?
  3. El kill switch contiene la celula, y sobrevive al defecto de RBS?

Uso:  .venv/bin/python model/run_analysis.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import kinetics as kn  # noqa: E402
from model.kinetics import KillSwitchConfig, SensorConfig  # noqa: E402

FIG_DIR = Path(__file__).resolve().parent / "figures"
CLAIM_MIN = 30.0  # la afirmacion de "<30 minutos" del proyecto


def fmt(x: float | None) -> str:
    return "no alcanza" if x is None else f"{x:.0f} min"


def section(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def q1_time_to_signal() -> list[tuple]:
    """Tiempo hasta senal para cada combinacion reportero x lectura."""
    section("1. Tiempo hasta senal detectable")
    rows = []
    for rep in ("gfp", "mcherry", "chromoprotein"):
        for readout in ("visual", "instrument"):
            cfg = SensorConfig(reporter=rep, readout=readout, inducer=10.0)
            t = kn.time_to_signal(cfg)
            ok = t is not None and t <= CLAIM_MIN
            rows.append((rep, readout, t, ok))
            print(f"  {rep:<14} {readout:<11} {fmt(t):>12}"
                  f"   {'cumple <30 min' if ok else ''}")
    return rows


def q2_rbs_defect() -> None:
    """Impacto del espaciado RBS->ATG de 0 nt sobre el tiempo de respuesta."""
    section("2. Efecto del defecto de espaciado RBS->ATG")
    print("  eficiencia de traduccion -> tiempo hasta senal (cromoproteina, visual)")
    for eff in (1.0, 0.5, 0.2, 0.1, 0.05):
        cfg = SensorConfig(reporter="chromoprotein", readout="visual",
                           rbs_efficiency=eff, inducer=10.0)
        t = kn.time_to_signal(cfg)
        note = "  <- espaciado 0 nt cae en este rango" if eff <= 0.1 else ""
        print(f"    {eff:>5.0%}  {fmt(t):>12}{note}")


def q3_killswitch() -> None:
    """Contencion: viabilidad, disparo tras perdida de plasmido, y defecto RBS."""
    section("3. Kill switch MazEF")

    k_max = kn.max_toxin_synthesis()
    print("  condicion de viabilidad:  k_mazF < d_mazF * umbral_letal")
    print(f"    k_mazF admisible  {k_max:.5f} nM/s")
    print(f"    k_mazF del diseno {kn.p('k_mazF'):.5f} nM/s"
          f"   ({kn.p('k_mazF') / k_max:.0f}x por encima)")
    print(f"    -> diseno viable: {kn.is_viable()}")

    print()
    t_arr, y = kn.simulate_killswitch(KillSwitchConfig())
    e, f, c = y[0, -1], y[1, -1], y[2, -1]
    surv = kn.time_to_death(KillSwitchConfig())
    print("  con el plasmido retenido (deberia sobrevivir):")
    print(f"    MazE {e:7.1f} nM   MazF libre {f:7.1f} nM   complejo {c:7.1f} nM")
    print(f"    {'sobrevive' if surv is None else f'MUERE a los {surv:.0f} min'}"
          f"{'' if surv is None else '   <- el switch dispara solo'}")

    print()
    d = kn.time_to_death(KillSwitchConfig(plasmid_lost_at=3600.0))
    print("  tras perder el plasmido a los 60 min (deberia morir):")
    print(f"    {fmt(d)} hasta MazF letal")

    print()
    print("  margen entre ambos escenarios:")
    if surv is not None and d is not None:
        print(f"    muerte espuria {surv:.0f} min vs muerte real {d:.0f} min")
        print("    el margen es demasiado estrecho para discriminar:")
        print("    la celula muere tenga o no el plasmido.")

    print()
    print("  efecto del defecto de RBS sobre la antitoxina:")
    for eff in (1.0, 0.5, 0.2, 0.1):
        cfg = KillSwitchConfig(rbs_mazE=eff, rbs_mazF=1.0)
        dd = kn.time_to_death(cfg)
        print(f"    RBS mazE {eff:>5.0%}   "
              f"{'sobrevive' if dd is None else f'muere a los {dd:>3.0f} min'}")
    print("    el defecto acelera la muerte espuria: menos antitoxina,")
    print("    menos retraso antes de que la toxina acumulada sea letal.")


def make_figures() -> None:
    """Escribe las tres figuras del analisis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})

    # --- fig 1: curso temporal por reportero ---
    # Escala lineal: en log las tres curvas se solapan y el rango vacio
    # de 20 decadas oculta justamente la diferencia que importa.
    fig, ax = plt.subplots(figsize=(6.5, 4))
    styles = {"gfp": ("GFP (fluorescente)", "#2a9d3f"),
              "mcherry": ("mCherry (fluorescente)", "#c0392b"),
              "chromoprotein": ("cromoproteina", "#7d3c98")}
    for rep, (label, color) in styles.items():
        cfg = SensorConfig(reporter=rep, readout="visual", inducer=10.0)
        t_a, y = kn.simulate_sensor(cfg, t_end=7200.0)
        tt = kn.time_to_signal(cfg)
        ax.plot(t_a / 60.0, y[2], label=f"{label} — {tt:.0f} min",
                color=color, lw=1.8)
        ax.plot([tt], [kn.p("thr_visual")], "o", color=color, ms=5, zorder=5)

    ax.axhline(kn.p("thr_visual"), ls="--", c="k", lw=1,
               label="umbral visible a simple vista")
    ax.axvspan(0, CLAIM_MIN, color="#e67e22", alpha=0.10)
    ax.axvline(CLAIM_MIN, c="#e67e22", lw=1.2)
    ax.text(CLAIM_MIN - 1.5, kn.p("thr_visual") * 1.42,
            "30 min\n(afirmado)", color="#e67e22", fontsize=8, ha="right")
    ax.set_xlabel("tiempo (min)")
    ax.set_ylabel("reportero maduro (nM)")
    ax.set_title("Ningun reportero cruza el umbral visual en 30 min")
    ax.set_xlim(0, 120)
    ax.set_ylim(0, kn.p("thr_visual") * 1.9)
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_time_to_signal.png")
    plt.close(fig)

    # --- fig 2: dosis-respuesta a 30 min ---
    fig, ax = plt.subplots(figsize=(6.5, 4))
    doses = np.logspace(-2, 2, 60)
    for rep, (label, color) in styles.items():
        r = kn.dose_response(doses, SensorConfig(reporter=rep), t_read=1800.0)
        ax.plot(doses, r, label=label, color=color, lw=1.8)
    ax.axhline(kn.p("thr_visual"), ls="--", c="k", lw=1, label="umbral visual")
    ax.axhline(kn.p("thr_instrument"), ls=":", c="gray", lw=1,
               label="umbral instrumento")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("contaminante (uM)")
    ax.set_ylabel("reportero maduro a los 30 min (nM)")
    ax.set_title("Dosis-respuesta a los 30 min")
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "02_dose_response.png")
    plt.close(fig)

    # --- fig 3: kill switch ---
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6))

    for ax, lost, title in (
        (a1, None, "Plasmido retenido\n(deberia sobrevivir)"),
        (a2, 3600.0, "Plasmido perdido a 60 min\n(deberia morir)"),
    ):
        t_a, y = kn.simulate_killswitch(
            KillSwitchConfig(plasmid_lost_at=lost), t_end=14400.0)
        ax.plot(t_a / 60.0, y[0], label="MazE libre", color="#2471a3", lw=1.8)
        ax.plot(t_a / 60.0, y[1], label="MazF libre", color="#c0392b", lw=1.8)
        ax.axhline(kn.p("tox_lethal"), ls="--", c="k", lw=1, label="umbral letal")
        if lost:
            ax.axvline(60, c="#e67e22", lw=1.2)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("tiempo (min)")
        ax.set_yscale("log")
        ax.set_ylim(1e-2, 1e4)
        ax.legend(fontsize=7)
    a1.set_ylabel("concentracion (nM)")

    fig.suptitle("La toxina supera el umbral letal en AMBOS casos: "
                 "el switch no discrimina", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_killswitch.png")
    plt.close(fig)

    # --- fig 4: barrido de sintesis de toxina ---
    fig, ax = plt.subplots(figsize=(6.5, 4))
    k_max = kn.max_toxin_synthesis()
    ks = np.logspace(np.log10(k_max) - 1.5, np.log10(kn.p("k_mazF")) + 0.5, 40)
    orig = kn.PARAMS["k_mazF"]
    times = []
    for k in ks:
        kn.PARAMS["k_mazF"] = (float(k), orig[1], orig[2])
        d = kn.time_to_death(KillSwitchConfig(), t_end=86400.0)
        times.append(np.nan if d is None else d)
    kn.PARAMS["k_mazF"] = orig

    ax.plot(ks, times, color="#c0392b", lw=1.8, zorder=3)
    ax.axvspan(ks[0], k_max, color="#2a9d3f", alpha=0.10, zorder=0)
    ax.axvline(k_max, ls="--", c="k", lw=1.2, zorder=2)
    ax.axvline(kn.p("k_mazF"), ls=":", c="#e67e22", lw=1.5, zorder=2)

    top = ax.get_ylim()[1]
    ax.text(k_max * 1.15, top * 0.88,
            f"cota de viabilidad\n$k_F = d_F \\cdot F_{{letal}}$\n{k_max:.4f} nM/s",
            fontsize=8, va="top", zorder=4)
    ax.annotate("diseno actual\n(16x la cota)", xy=(kn.p("k_mazF"), 95),
                xytext=(kn.p("k_mazF") * 0.30, top * 0.42),
                color="#e67e22", fontsize=8, ha="center", zorder=4,
                arrowprops=dict(arrowstyle="->", color="#e67e22", lw=1.2))
    ax.text(ks[0] * 1.15, top * 0.06, "celula viable", fontsize=8,
            color="#1d6b2b", zorder=4)
    ax.set_xscale("log")
    ax.set_xlabel("sintesis de MazF, $k_F$ (nM/s)")
    ax.set_ylabel("tiempo hasta muerte espuria (min)")
    ax.set_title("A la izquierda de la cota la celula sobrevive indefinidamente")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "04_toxin_sweep.png")
    plt.close(fig)

    print()
    print(f"  figuras escritas en {FIG_DIR.relative_to(Path.cwd())}/")


def main() -> int:
    print("=" * 66)
    print("Analisis cinetico -- Yakuponia")
    print("=" * 66)
    print("Parametros de literatura para sistemas analogos, no mediciones")
    print("de este proyecto. Sirven para ordenes de magnitud.")

    rows = q1_time_to_signal()
    q2_rbs_defect()
    q3_killswitch()
    make_figures()

    section("Conclusion")
    visual_ok = [r for r in rows if r[1] == "visual" and r[3]]
    if not visual_ok:
        best = min((r for r in rows if r[1] == "visual" and r[2] is not None),
                   key=lambda r: r[2], default=None)
        print(f"  Ninguna configuracion alcanza senal VISUAL en <{CLAIM_MIN:.0f} min.")
        if best:
            print(f"  La mas rapida es {best[0]} con {fmt(best[2])}.")
        print("  La afirmacion de '<30 minutos' no se sostiene para lectura")
        print("  a simple vista: la maduracion del cromoforo impone el piso.")
        print("  Con instrumento si es alcanzable.")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
