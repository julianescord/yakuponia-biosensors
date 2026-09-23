#!/usr/bin/env python3
"""Figuras suplementarias, una por panel.

Las figuras principales fusionan varios paneles en una sola afirmacion. Estas
conservan cada pieza por separado, para el material suplementario:

  S1  curso temporal del sensor          (panel a de la Fig. 1)
  S2  dosis-respuesta a los 30 min       (panel b de la Fig. 1)
  S3  barrido de sintesis de toxina      (antes figura 04, Extended Data)
  S4  reparto por cuadrantes             (barras horizontales)

La S4 vuelve al formato de barras de la version original: en la Fig. 3 el
reparto ya lo muestra el panel b, asi que aqui no compite con nada y las
cuatro categorias se leen mejor en barras que en una franja apilada.

Uso:  .venv/bin/python model/figS_supplementary.py [--fast]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import kinetics as kn  # noqa: E402
from model.kinetics import KillSwitchConfig  # noqa: E402

CENSOR_MIN = 360.0


def panel_quadrants(Y):
    """S4 — fraccion de muestras en cada cuadrante, en barras horizontales."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

    safe, kills = Y[:, 1] > 0, Y[:, 2] < CENSOR_MIN
    cats = [("Cumple ambas", safe & kills, "v3"),
            ("Segura, no contiene", safe & ~kills, "mazF"),
            ("Contiene, dispara sola", ~safe & kills, "v1"),
            ("Ninguna", ~safe & ~kills, "muted")]

    fig = plt.figure(figsize=(4.40, 1.45))
    ax = fig.add_axes([0.335, 0.215, 0.615, 0.735])

    ypos = np.arange(len(cats))[::-1]
    for yp, (lab, mask, role) in zip(ypos, cats):
        frac = mask.mean()
        ax.barh(yp, frac * 100, height=0.62, color=fs.c(role), zorder=2)
        # Una categoria vacia no dibuja barra: se marca el cero para que no
        # se confunda con un dato que falta.
        ax.text(frac * 100 + 1.8, yp, f"{frac:.0%}", va="center",
                color="black", zorder=3)

    ax.set_yticks(ypos)
    ax.set_yticklabels([c[0] for c in cats])
    ax.tick_params(axis="y", length=0, pad=2.0)
    ax.set_ylim(-0.7, len(cats) - 0.3)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Muestras (%)")
    return fig


def panel_toxin_sweep():
    """S3 — barrido de sintesis de MazF frente a la cota de viabilidad."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

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

    fig = plt.figure(figsize=(4.40, 1.75))
    ax = fig.add_axes([0.145, 0.195, 0.825, 0.640])

    ax.plot(ks, times, color=fs.c("mazF"), lw=1.4, zorder=4)
    ax.axvspan(ks[0], k_max, color=fs.c("ok"), alpha=0.075, zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(ks[0], ks[-1])
    # Sin exponentes: mathtext dibuja el superindice al 70% del cuerpo y cae
    # bajo el piso de 5 pt de Nature. La notacion decimal evita el problema.
    ax.set_xticks([1e-3, 1e-2, 1e-1])
    ax.set_xticklabels(["0.001", "0.01", "0.1"])
    top = max(t for t in times if not np.isnan(t)) * 1.06
    ax.set_ylim(0, top)

    # Las cotas se rotulan por encima del area de datos: dentro cruzaban la
    # curva y competian entre si.
    marks = [(k_max, fs.c("threshold"), (0, (4, 2)), f"Cota v1\n{k_max:.4f}",
              1.04, "center"),
             (kn.p("k_mazF"), fs.c("event"), (0, (1.5, 1.5)),
              f"Diseno\n{kn.p('k_mazF'):.2f}", 1.04, "right"),
             (k_max_v2, fs.c("ok"), (0, (4, 2)),
              f"Cota v2\n{k_max_v2:.3f}", 1.30, "left")]
    for xv, color, ls, txt, yf, ha in marks:
        ax.axvline(xv, ls=ls, color=color, lw=0.8, zorder=3)
        ax.text(xv, top * yf, txt, color=color, ha=ha, va="bottom",
                linespacing=1.3, zorder=6, clip_on=False)

    ax.text(ks[0] * 1.2, top * 0.09, "Celula viable", color=fs.c("ok"),
            zorder=5)
    ax.set_xlabel("Sintesis de MazF (nM/s)")
    ax.set_ylabel("Tiempo hasta\nmuerte espuria (min)")
    return fig


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="usa la cache del muestreo reducido")
    args = ap.parse_args()

    from model import figstyle as fs
    from model.fig1_sensor import panel_dose, panel_timecourse
    from model.fig3_sensitivity import load

    fs.apply_style()
    written = []

    # S1 y S2 reutilizan los generadores de la Fig. 1: son el mismo panel,
    # no una version distinta del mismo grafico.
    fig_a, _ = panel_timecourse()
    written.append(fs.save_panel(fig_a, "figS1_timecourse"))
    written.append(fs.save_panel(panel_dose(), "figS2_dose_response"))
    written.append(fs.save_panel(panel_toxin_sweep(), "figS3_toxin_sweep"))

    _, _, Y = load(args.fast)
    written.append(fs.save_panel(panel_quadrants(Y), "figS4_quadrants"))

    print(f"  {len(written)} paneles suplementarios en figures/panels/")
    for p in written:
        print(f"    {p.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
