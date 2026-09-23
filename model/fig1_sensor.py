#!/usr/bin/env python3
"""Fig. 1 — El sensor no alcanza senal visible en 30 minutos.

Afirmacion de la figura (una sola, y los dos paneles la sirven):
  ninguna version del biosensor produce color visible a simple vista dentro
  de los 30 min que afirma el proyecto, y el limite no es de concentracion.

Cadena de evidencia:
  a  evidencia primaria -- curso temporal de las tres versiones frente al
     umbral visual; solo la v3 lo cruza, y lo hace a los 50 min
  b  control de la explicacion alternativa -- dosis-respuesta a los 30 min:
     ningun reportero cruza el umbral a ninguna concentracion, asi que el
     techo lo impone la maduracion del cromoforo, no la dosis

Las dos eran figuras separadas (01 y 02) que afirmaban lo mismo. Fusionarlas
deja una sola afirmacion con dos roles inferenciales distintos.

Uso:  .venv/bin/python model/fig1_sensor.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import kinetics as kn  # noqa: E402
from model.kinetics import SensorConfig  # noqa: E402

CLAIM_MIN = 30.0     # la afirmacion de "<30 minutos" del proyecto
T_MAX = 140.0        # ventana visible del panel a

# Cada version con su reportero, espaciado y apareamiento SD del RBS.
VERSIONS = {
    "v1": dict(reporter="gfp", spacer_ok=False, sd=4),
    "v2": dict(reporter="chromoprotein", spacer_ok=True, sd=4),
    "v3": dict(reporter="chromoprotein", spacer_ok=True, sd=6),
}


def eff_of(v: dict) -> float:
    """Eficiencia traduccional: apareamiento SD, penalizado si falta espaciador."""
    e = kn.rbs_efficiency_from_sd(v["sd"])
    return e if v["spacer_ok"] else e * 0.67


def panel_timecourse():
    """a — reportero maduro frente al tiempo, por version."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

    thr = kn.p("thr_visual")
    # El area de ejes deja margen a la derecha: los rotulos directos
    # viven dentro del lienzo del panel, que no se recorta con `tight`.
    fig = plt.figure(figsize=(2.55, 1.70))
    ax = fig.add_axes([0.165, 0.185, 0.645, 0.775])

    ends = {}
    for ver, v in VERSIONS.items():
        cfg = SensorConfig(reporter=v["reporter"], readout="visual",
                           rbs_efficiency=eff_of(v), inducer=10.0)
        t, y = kn.simulate_sensor(cfg, t_end=9000.0)
        mins = t / 60.0
        # Se recorta la serie al rango visible en vez de dejar que matplotlib
        # la recorte: una curva recortada conserva una caja que abarca el
        # panel entero y colisiona con cualquier texto.
        vis = mins <= T_MAX
        hero = ver == "v3"
        ax.plot(mins[vis], y[2][vis], color=fs.c(ver),
                lw=1.4 if hero else 1.0,
                ls="-" if hero else (0, (3, 1.6)),
                zorder=4 if hero else 3)
        tt = kn.time_to_signal(cfg)
        ends[ver] = (mins[vis][-1], y[2][vis][-1], tt)
        if tt:
            ax.plot([tt], [thr], "o", color=fs.c(ver), ms=2.8, zorder=6,
                    markeredgecolor="white", markeredgewidth=0.4)

    ymax = max(e[1] for e in ends.values()) * 1.06
    ax.axhline(thr, ls=(0, (4, 2)), color=fs.c("threshold"), lw=0.6, zorder=2)
    # La ventana afirmada se marca con una vertical punteada, no con un
    # sombreado: el relleno competia con las curvas.
    ax.axvline(CLAIM_MIN, ls=(0, (1.5, 1.5)), color=fs.c("threshold"),
               lw=0.7, zorder=2)

    # Rotulos directos al final de cada curva. La v1 y la v2 terminan a menos
    # de una altura de texto entre si, asi que se separan un minimo en vez de
    # anclarse al valor exacto.
    order = sorted(ends, key=lambda v: ends[v][1])
    gap = ymax * 0.085
    ys, prev = {}, None
    for ver in order:
        y = ends[ver][1] if prev is None else max(ends[ver][1], prev + gap)
        ys[ver], prev = y, y
    for ver, y in ys.items():
        ax.text(T_MAX * 1.02, y, ver, color=fs.c(ver),
                va="center", ha="left", clip_on=False)

    ax.set_xlim(0, T_MAX)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 30, 60, 90, 120])
    ax.set_yticks([0, 5000, 10000])
    ax.set_yticklabels(["0", "5000", "10 000"])
    ax.set_xlabel("Tiempo (min)")
    ax.set_ylabel("Reportero maduro (nM)")
    return fig, ends


def panel_dose():
    """b — dosis-respuesta a los 30 min, por reportero."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

    fig = plt.figure(figsize=(2.55, 1.70))
    ax = fig.add_axes([0.175, 0.185, 0.56, 0.775])

    doses = np.logspace(-2, 2, 60)
    series = [("gfp", "GFP", "v1"),
              ("mcherry", "mCherry", "v2"),
              ("chromoprotein", "Cromoproteina", "v3")]
    # Las tres mesetas quedan a menos de media decada, demasiado juntas para
    # rotular a la altura exacta de cada curva. Los rotulos se reparten en
    # posiciones fijas del eje, en el mismo orden vertical que las curvas.
    slots = {"chromoprotein": 0.90, "gfp": 0.79, "mcherry": 0.68}
    for rep, label, role in series:
        r = kn.dose_response(doses, SensorConfig(reporter=rep), t_read=1800.0)
        hero = rep == "chromoprotein"
        ax.plot(doses, r, color=fs.c(role), lw=1.4 if hero else 1.0,
                ls="-" if hero else (0, (3, 1.6)), zorder=4 if hero else 3)
        ax.text(1.06, slots[rep], label, transform=ax.transAxes,
                color=fs.c(role), va="center", ha="left", clip_on=False)

    for val in (kn.p("thr_visual"), kn.p("thr_instrument")):
        ax.axhline(val, ls=(0, (4, 2)), color=fs.c("threshold"), lw=0.6,
                   zorder=2)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(doses[0], doses[-1])
    # Sin exponentes: mathtext los dibuja al 70% del cuerpo y caen bajo el
    # piso de 5 pt; los superindices Unicode faltan en Liberation Sans.
    ax.set_xticks([0.01, 1, 100])
    ax.set_xticklabels(["0.01", "1", "100"])
    ax.set_yticks([1, 100, 10000])
    ax.set_yticklabels(["1", "100", "10 000"])
    ax.set_xlabel("Contaminante (µM)")
    ax.set_ylabel("Reportero a los 30 min (nM)")
    return fig


def main() -> int:
    from model import figstyle as fs

    fs.apply_style()
    written = []

    fig_a, ends = panel_timecourse()
    written.append(fs.save_panel(fig_a, "fig1_timecourse"))
    written.append(fs.save_panel(panel_dose(), "fig1_dose"))

    print(f"  {len(written)} paneles escritos en figures/panels/")
    for p in written:
        print(f"    {p.name}")
    print()
    for ver, (_, _, tt) in ends.items():
        print(f"    {ver}: {'no alcanza' if tt is None else f'{tt:.0f} min'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
