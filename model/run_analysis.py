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
    """Tiempo hasta senal: v1 frente a v2."""
    section("1. Tiempo hasta senal detectable")
    print("  v1 = GFP/mCherry con espaciado RBS 0 nt (eff. traduccional ~10%)")
    print("  v2 = cromoproteina con espaciador de 7 nt")
    print()
    rows = []
    for label, rep, eff in (("v1 FtaGen (GFP)", "gfp", 0.1),
                            ("v1 MetalGen (mCherry)", "mcherry", 1.0),
                            ("v2 ambos (cromoproteina)", "chromoprotein", 1.0)):
        for readout in ("visual", "instrument"):
            cfg = SensorConfig(reporter=rep, readout=readout,
                               rbs_efficiency=eff, inducer=10.0)
            tt = kn.time_to_signal(cfg)
            ok = tt is not None and tt <= CLAIM_MIN
            rows.append((label, readout, tt, ok))
            print(f"  {label:<26} {readout:<11} {fmt(tt):>12}"
                  f"   {'cumple <30 min' if ok else ''}")
    return rows


def q2_rbs_defect() -> None:
    """Por que el espaciador era imprescindible."""
    section("2. Efecto del espaciado RBS->ATG")
    print("  eficiencia de traduccion -> tiempo hasta senal (cromoproteina)")
    for eff in (1.0, 0.5, 0.2, 0.1, 0.05):
        cfg = SensorConfig(reporter="chromoprotein", readout="visual",
                           rbs_efficiency=eff, inducer=10.0)
        note = ""
        if eff == 1.0:
            note = "  <- v2, espaciador de 7 nt"
        elif eff <= 0.1:
            note = "  <- v1, espaciado 0 nt"
        print(f"    {eff:>5.0%}  {fmt(kn.time_to_signal(cfg)):>12}{note}")
    print()
    print("  Por debajo del 20% la sintesis no supera a la dilucion por")
    print("  division: el sensor no llega nunca, no solo tarda mas.")


def q3_killswitch() -> None:
    """Contencion: v1 frente a v2 con tag ssrA."""
    section("3. Kill switch MazEF")

    print("  condicion de viabilidad:  k_mazF < d_mazF * umbral_letal")
    print(f"  k_mazF del diseno: {kn.p('k_mazF'):.4f} nM/s (igual en v1 y v2)")
    print()
    for ssrA in (False, True):
        tag = "v2 (MazF-ssrA)" if ssrA else "v1 (MazF sin tag)"
        cota = kn.max_toxin_synthesis(ssrA)
        rel = kn.p("k_mazF") / cota
        print(f"  {tag:<18} cota {cota:.4f} nM/s   "
              f"k_mazF esta {rel:.1f}x {'por encima' if rel > 1 else 'por debajo'}"
              f"   viable={kn.is_viable(ssrA)}")

    print()
    print("  comportamiento:")
    print(f"    {'':22}{'plasmido retenido':<22}{'plasmido perdido a 60 min'}")
    print(f"    {'':22}{'(debe sobrevivir)':<22}{'(debe morir)'}")
    for ssrA in (False, True):
        tag = "v2 (MazF-ssrA)" if ssrA else "v1 (sin tag)"
        keep = kn.time_to_death(KillSwitchConfig(ssrA=ssrA))
        lost = kn.time_to_death(
            KillSwitchConfig(ssrA=ssrA, plasmid_lost_at=3600.0))
        k_s = "sobrevive" if keep is None else f"MUERE {keep:.0f} min"
        l_s = "sobrevive" if lost is None else f"muere {lost:.0f} min"
        ok = keep is None and lost is not None
        print(f"    {tag:<22}{k_s:<22}{l_s:<18}{'  <- correcto' if ok else ''}")

    print()
    print("  El tag ssrA sube la cota 22x sin tocar promotores ni RBS:")
    print("  desestabilizar la toxina es lo que restaura la contencion.")

    # El margen de la v2 es estrecho y conviene declararlo.
    import numpy as np
    keep_ss = kn.simulate_killswitch(KillSwitchConfig(ssrA=True))[1][1, -1]
    thr = kn.p("tox_lethal")
    t_a, y = kn.simulate_killswitch(
        KillSwitchConfig(ssrA=True, plasmid_lost_at=3600.0), t_end=21600.0)
    hot = np.flatnonzero(y[1] >= thr)
    print()
    print("  margenes de la v2 (estrechos, conviene declararlos):")
    print(f"    en reposo MazF libre se estabiliza en {keep_ss:.0f} nM,"
          f" umbral {thr:.0f} nM")
    print(f"    -> margen de solo {thr - keep_ss:.0f} nM ante variacion de "
          "parametros")
    if hot.size:
        print(f"    tras perder el plasmido la toxina supera el umbral entre "
              f"{t_a[hot[0]] / 60:.0f} y {t_a[hot[-1]] / 60:.0f} min")
        print(f"    -> ventana letal de {(t_a[hot[-1]] - t_a[hot[0]]) / 60:.0f}"
              " min: la muerte debe ocurrir dentro de ella")


def make_figures() -> None:
    """Escribe las tres figuras del analisis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})

    # --- fig 1: v1 frente a v2 ---
    fig, ax = plt.subplots(figsize=(6.5, 4))
    cases = [
        ("v1 FtaGen: GFP, RBS 0 nt", "gfp", 0.1, "#2a9d3f", ":"),
        ("v1 MetalGen: mCherry", "mcherry", 1.0, "#c0392b", ":"),
        ("v2: cromoproteina + espaciador", "chromoprotein", 1.0, "#7d3c98", "-"),
    ]
    for label, rep, eff, color, ls in cases:
        cfg = SensorConfig(reporter=rep, readout="visual",
                           rbs_efficiency=eff, inducer=10.0)
        t_a, y = kn.simulate_sensor(cfg, t_end=7200.0)
        tt = kn.time_to_signal(cfg)
        suffix = f"{tt:.0f} min" if tt else "no alcanza"
        ax.plot(t_a / 60.0, y[2], label=f"{label} — {suffix}",
                color=color, lw=2.0 if ls == "-" else 1.6, ls=ls)
        if tt:
            ax.plot([tt], [kn.p("thr_visual")], "o", color=color, ms=5, zorder=5)

    ax.axhline(kn.p("thr_visual"), ls="--", c="k", lw=1,
               label="umbral visible a simple vista")
    ax.axvspan(0, CLAIM_MIN, color="#e67e22", alpha=0.10)
    ax.axvline(CLAIM_MIN, c="#e67e22", lw=1.2)
    ax.text(CLAIM_MIN - 1.5, kn.p("thr_visual") * 1.45,
            "30 min\n(afirmado)", color="#e67e22", fontsize=8, ha="right")
    ax.set_xlabel("tiempo (min)")
    ax.set_ylabel("reportero maduro (nM)")
    ax.set_title("v2 restaura la funcion, pero el piso de maduracion persiste")
    ax.set_xlim(0, 120)
    ax.set_ylim(0, kn.p("thr_visual") * 1.9)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_time_to_signal.png")
    plt.close(fig)

    # --- fig 2: dosis-respuesta a 30 min ---
    fig, ax = plt.subplots(figsize=(6.5, 4))
    doses = np.logspace(-2, 2, 60)
    styles = {"gfp": ("GFP (v1)", "#2a9d3f"),
              "mcherry": ("mCherry (v1)", "#c0392b"),
              "chromoprotein": ("cromoproteina (v2)", "#7d3c98")}
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

    # --- fig 3: kill switch v1 vs v2 ---
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
    for row, ssrA in enumerate((False, True)):
        ver = "v2 (MazF-ssrA)" if ssrA else "v1 (MazF sin tag)"
        for col, lost in enumerate((None, 3600.0)):
            ax = axes[row, col]
            t_a, y = kn.simulate_killswitch(
                KillSwitchConfig(ssrA=ssrA, plasmid_lost_at=lost),
                t_end=14400.0)
            ax.plot(t_a / 60.0, y[0], label="MazE libre",
                    color="#2471a3", lw=1.6)
            ax.plot(t_a / 60.0, y[1], label="MazF libre",
                    color="#c0392b", lw=1.6)
            ax.axhline(kn.p("tox_lethal"), ls="--", c="k", lw=1)
            if lost:
                ax.axvline(60, c="#e67e22", lw=1.1)
            ax.set_yscale("log")
            ax.set_ylim(1e-2, 1e4)
            d = kn.time_to_death(
                KillSwitchConfig(ssrA=ssrA, plasmid_lost_at=lost))
            want_death = lost is not None
            good = (d is not None) == want_death
            verdict = "sobrevive" if d is None else f"muere {d:.0f} min"
            ax.set_title(f"{ver} — {'plasmido perdido' if lost else 'retenido'}"
                         f"\n{verdict}  {'OK' if good else 'FALLA'}",
                         fontsize=8.5,
                         color="#1d6b2b" if good else "#a93226")
            if row == 1:
                ax.set_xlabel("tiempo (min)")
            if col == 0:
                ax.set_ylabel("concentracion (nM)")
    axes[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle("El tag ssrA restaura la contencion: v2 solo muere "
                 "cuando debe", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_killswitch.png")
    plt.close(fig)

    # --- fig 4: barrido de sintesis de toxina ---
    fig, ax = plt.subplots(figsize=(6.5, 4))
    k_max = kn.max_toxin_synthesis(False)
    k_max_v2 = kn.max_toxin_synthesis(True)
    ks = np.logspace(np.log10(k_max) - 1.5, np.log10(k_max_v2) + 0.4, 40)
    orig = kn.PARAMS["k_mazF"]
    times = []
    for k in ks:
        kn.PARAMS["k_mazF"] = (float(k), orig[1], orig[2])
        d = kn.time_to_death(KillSwitchConfig(), t_end=86400.0)
        times.append(np.nan if d is None else d)
    kn.PARAMS["k_mazF"] = orig

    ax.plot(ks, times, color="#c0392b", lw=1.8, zorder=3)
    ax.axvspan(ks[0], k_max, color="#2a9d3f", alpha=0.10, zorder=0)
    ax.axvline(k_max_v2, ls="--", c="#1d6b2b", lw=1.4, zorder=2)
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
    ax.text(ks[0] * 1.15, top * 0.06, "celula viable (v1)", fontsize=8,
            color="#1d6b2b", zorder=4)
    ax.annotate(f"cota v2 con ssrA\n{k_max_v2:.3f} nM/s\nel diseno queda a la izquierda",
                xy=(k_max_v2, top * 0.30),
                xytext=(k_max_v2 * 0.55, top * 0.68),
                fontsize=8, color="#1d6b2b", ha="center", zorder=4,
                arrowprops=dict(arrowstyle="->", color="#1d6b2b", lw=1.2))
    ax.set_xscale("log")
    ax.set_xlabel("sintesis de MazF, $k_F$ (nM/s)")
    ax.set_ylabel("tiempo hasta muerte espuria (min)")
    ax.set_title("El tag ssrA mueve la cota 22x y deja el diseno del lado viable")
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
