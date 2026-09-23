#!/usr/bin/env python3
"""Fig. 2 — Solo la v3 contiene la celula.

Afirmacion de la figura (una sola, y todos los paneles la sirven):
  el kill switch solo discrimina entre retener y perder el plasmido en la v3.

Cadena de evidencia, un rol inferencial por fila:
  a, b  v1 -- evidencia primaria: sin ssrA la toxina mata con el plasmido puesto
  c, d  v2 -- control: el ssrA evita el disparo espurio pero pierde la contencion
  e, f  v3 -- caso decisivo: sobrevive con el plasmido y muere al perderlo

El barrido de sintesis de toxina (antes figura 04) sale de aqui: repetia esta
misma afirmacion por via analitica. Va a Extended Data como `ED1`.

Cada panel se dibuja por separado y se exporta a SVG; el tamano fisico final
lo fija el compositor en mm. Aqui no hay titulos ni texto explicativo: eso
vive en `figures/legends.md`.

Tampoco hay leyenda de color: el pie ya dice que MazE va en azul y MazF en
vino, y el bloque de leyenda dejaba una franja vacia bajo la rejilla.

Uso:  .venv/bin/python model/fig2_killswitch.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import kinetics as kn  # noqa: E402
from model.kinetics import KillSwitchConfig  # noqa: E402

# (version, SD de mazE, SD de mazF, tag ssrA)
VERSIONS = [("v1", 5, 4, False), ("v2", 5, 4, True), ("v3", 7, 6, True)]
# (momento de perdida del plasmido, si aplica)
FATES = [None, 3600.0]

T_END = 21600.0
FLOOR = 6e-3   # piso visible; por debajo se enmascara en vez de recortar


def draw_panel(sd_e: int, sd_f: int, ssrA: bool, lost: float | None,
               *, show_y: bool, show_x: bool):
    """Un panel: MazE y MazF libres frente al tiempo, en escala log."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

    lethal = kn.p("tox_lethal")
    cfg = KillSwitchConfig(rbs_mazE=kn.rbs_efficiency_from_sd(sd_e),
                           rbs_mazF=kn.rbs_efficiency_from_sd(sd_f),
                           ssrA=ssrA, plasmid_lost_at=lost)
    t, y = kn.simulate_killswitch(cfg, t_end=T_END)
    mins = t / 60.0

    # El area de ejes se fija en coordenadas absolutas y el lienzo se recorta
    # por igual en los seis. Con `tight` y margenes automaticos, un panel con
    # eje Y sale mas ancho que uno sin el y las columnas dejan de alinear.
    fig = plt.figure(figsize=(2.55, 1.38))
    ax = fig.add_axes([0.175, 0.225, 0.80, 0.735])

    # Se enmascara lo que cae bajo el piso visible en vez de dejar que
    # matplotlib lo recorte: una curva recortada conserva una caja que llega
    # al borde del panel y pisa los rotulos del eje.
    for series, role in ((y[0], "mazE"), (y[1], "mazF")):
        ax.plot(mins, np.where(series >= FLOOR, series, np.nan),
                color=fs.c(role), lw=1.1, zorder=3)

    ax.axhline(lethal, ls=(0, (3, 2)), color=fs.c("threshold"), lw=0.6,
               zorder=2)
    if lost is not None:
        ax.axvline(lost / 60.0, color=fs.c("event"), lw=0.7,
                   ls=(0, (1.5, 1.5)), zorder=2)

    ax.set_yscale("log")
    ax.set_ylim(4e-3, 1e4)
    ax.set_xlim(0, 360)
    # Sin exponentes: mathtext dibuja el superindice al 70% del cuerpo y
    # genera digitos de 3.7 pt, bajo el piso de 5 pt; y los superindices
    # Unicode (U+2070, U+207B) no existen en Liberation Sans, asi que salen
    # como cajas vacias. La notacion decimal evita ambos problemas y en este
    # rango se lee igual de bien.
    ax.set_yticks([1e-2, 1e0, 1e2, 1e4])
    ax.set_yticklabels(["0.01", "1", "100", "10 000"])
    ax.set_xticks([0, 120, 240, 360])

    # El veredicto es un dato, no un juicio: va en negro y solo dice que pasa.
    # Codificarlo en verde/bermellon tenia dos defectos: el bermellon era el
    # mismo de la curva de MazF, asi que un color significaba dos cosas en el
    # mismo panel, y el color quedaba como unica senal de acierto o fallo.
    # Que comportamiento es el buscado se explica en el pie de figura.
    d = kn.time_to_death(cfg)
    verdict = "Sobrevive" if d is None else f"Muere {d:.0f} min"
    ax.text(0.97, 0.90, verdict, transform=ax.transAxes, ha="right",
            va="top", color="black", zorder=6)

    if show_y:
        ax.set_ylabel("Proteina libre (nM)")
    else:
        ax.tick_params(labelleft=False)
    if show_x:
        ax.set_xlabel("Tiempo (min)")
    else:
        ax.tick_params(labelbottom=False)

    return fig


def main() -> int:
    from model import figstyle as fs

    fs.apply_style()
    written = []
    for row, (ver, sd_e, sd_f, ssrA) in enumerate(VERSIONS):
        for col, lost in enumerate(FATES):
            fig = draw_panel(sd_e, sd_f, ssrA, lost,
                             show_y=(col == 0), show_x=(row == 2))
            name = f"fig2_{ver}_{'lost' if lost else 'kept'}"
            written.append(fs.save_panel(fig, name))

    print(f"  {len(written)} paneles escritos en figures/panels/")
    for p in written:
        print(f"    {p.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
