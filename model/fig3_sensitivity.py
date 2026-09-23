#!/usr/bin/env python3
"""Fig. 3 — Las conclusiones bajo incertidumbre parametrica.

Afirmacion de la figura (una sola, y los tres paneles la sirven):
  las conclusiones del modelo no dependen por igual de los parametros, y la
  contencion falla sobre todo por no matar, no por disparar sola.

Cadena de evidencia:
  a  descomposicion -- que parametros gobiernan cada conclusion (Sobol)
  b  evidencia primaria -- la tension entre las dos condiciones del kill
     switch, sobre TODAS las muestras

El reparto por cuadrantes era un tercer panel y salio: cuantificaba las mismas
categorias que ya muestra `b`, asi que era un subconjunto suyo en otra forma
visual. Las cifras viven ahora en el pie de figura.

Reutiliza el barrido de Sobol cacheado en `model/.cache/`: el muestreo cuesta
minutos y las figuras se reajustan en segundos.

Uso:  .venv/bin/python model/fig3_sensitivity.py [--fast]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import sensitivity as sn  # noqa: E402

CENSOR_MIN = 360.0   # 6 h: las muestras que no matan antes quedan censuradas

# Nombre de manuscrito de cada parametro. Los identificadores del codigo
# (`k_txn_max`, `d_mazF_ssrA`) no van en una figura: las barras bajas y las
# abreviaturas son convenciones de programacion, no de publicacion.
PARAM_LABELS = {
    "k_txn_max": "Transcripcion maxima",
    "k_tln": "Traduccion por mRNA",
    "leak": "Expresion basal",
    "d_mrna": "Degradacion de mRNA",
    "d_prot": "Dilucion por division",
    "t_mat_chromo": "Maduracion del cromoforo",
    "thr_visual": "Umbral visual",
    "k_mazF": "Sintesis de MazF",
    "k_mazE": "Sintesis de MazE",
    "k_on": "Union MazE-MazF",
    "d_mazE": "Degradacion de MazE",
    "d_mazF_ssrA": "Degradacion de MazF-ssrA",
    "tox_lethal": "Umbral letal de MazF",
}

# Titulo de manuscrito de cada conclusion, mas corto que el del eje.
OUTPUT_TITLES = {
    "t_signal": "Tiempo hasta senal visible",
    "margin": "Margen bajo el umbral letal",
    "kill_delay": "Retardo hasta la muerte",
}


def load(fast: bool):
    """Indices de Sobol y salidas crudas, desde la cache."""
    from SALib.analyze import sobol as sobol_analyze

    prob = sn.problem()
    n = 128 if fast else 1024
    if not sn._cache_key(prob, n).exists():
        raise SystemExit(
            "no hay cache para esta configuracion; corre primero\n"
            "  .venv/bin/python model/sensitivity.py"
            + (" --fast" if fast else ""))
    _, Y = sn.sample_model(prob, n)
    res = {k: sobol_analyze.analyze(prob, Y[:, j], calc_second_order=False,
                                    print_to_console=False)
           for j, (k, _) in enumerate(sn.OUTPUTS)}
    return prob, res, Y


def panel_sobol(prob, results):
    """a — indices de Sobol por conclusion, en un solo panel agrupado."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

    # Un solo panel con los tres bloques apilados, en vez de tres subpaneles:
    # asi el grosor de barra es identico en todos y no se lee como
    # importancia.
    blocks = []
    for key, _ in sn.OUTPUTS:
        si = results[key]
        keep = [i for i in np.argsort(si["ST"]) if si["ST"][i] >= 0.02]
        blocks.append((OUTPUT_TITLES[key], keep, si))

    # Tres subpaneles en horizontal, uno por conclusion: cada titulo va sobre
    # su propio eje, que es donde calza sin pelear con la columna de nombres.
    # El eje y se fija al bloque mas poblado en todos, asi que el grosor de
    # barra es identico y no se lee como importancia.
    n_max = max(len(k) for _, k, _ in blocks)
    fig = plt.figure(figsize=(7.20, 0.125 * n_max + 0.86))
    gs = fig.add_gridspec(1, 3, wspace=0.60, left=0.145, right=0.995,
                          top=0.815, bottom=0.215)
    axes = [fig.add_subplot(gs[0, j]) for j in range(3)]

    for ax, (title, keep, si) in zip(axes, blocks):
        y = np.arange(len(keep))[::-1]
        # Escala de grises: estas barras miden magnitud, no categorias.
        ax.barh(y, [si["ST"][i] for i in keep], height=0.66,
                color="#D6D6D6", zorder=2)
        ax.barh(y, [max(si["S1"][i], 0) for i in keep], height=0.66,
                color="#4D4D4D", zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([PARAM_LABELS.get(prob["names"][i],
                                             prob["names"][i]) for i in keep])
        ax.tick_params(axis="y", length=0, pad=2.0)
        # tope comun: iguala el grosor de barra entre subpaneles
        # El eje se invierte para que el ST mayor quede arriba y las barras
        # se apoyen en el borde superior, no en el inferior.
        ax.set_ylim(n_max - 0.3, -0.7)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_title(title, fontweight="bold", color="black", pad=4)

    axes[1].set_xlabel("Indice de Sobol")

    # Clave bajo el subpanel derecho, donde el eje deja sitio. Nature pide la
    # clave dentro de la figura antes que describir los colores en el pie.
    ax = axes[2]
    ky = -0.30
    ax.add_patch(plt.Rectangle((0.10, ky - 0.030), 0.14, 0.060,
                               transform=ax.transAxes, clip_on=False,
                               facecolor="#4D4D4D", zorder=5))
    ax.text(0.28, ky, "S1", transform=ax.transAxes, ha="left", va="center",
            color="black", clip_on=False)
    ax.add_patch(plt.Rectangle((0.50, ky - 0.030), 0.14, 0.060,
                               transform=ax.transAxes, clip_on=False,
                               facecolor="#D6D6D6", zorder=5))
    ax.text(0.68, ky, "ST", transform=ax.transAxes, ha="left", va="center",
            color="black", clip_on=False)
    return fig


def panel_tradeoff(Y):
    """b — margen en reposo frente a retardo, con TODAS las muestras.

    Las censuradas (nunca alcanzan nivel letal en 6 h) van en una banda
    declarada sobre el eje, separada por un corte. Excluirlas, como se hacia
    antes, escondia justamente el grupo que demuestra la afirmacion de la
    figura: a medida que el margen en reposo crece, la toxina deja de matar.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import figstyle as fs

    margin, delay = Y[:, 1], Y[:, 2]
    safe, kills = margin > 0, delay < CENSOR_MIN

    # Ancho completo, igual que el panel a: la nube gana resolucion
    # horizontal y desaparece la franja vacia de la derecha.
    fig = plt.figure(figsize=(7.20, 2.05))
    ax = fig.add_axes([0.088, 0.185, 0.905, 0.600])   # nube
    axc = fig.add_axes([0.088, 0.845, 0.905, 0.075])  # banda de censura

    for a in (ax, axc):
        a.set_xlim(-120, 205)
        a.axvline(0, color=fs.c("threshold"), lw=0.6, ls=(0, (4, 2)),
                  zorder=4)

    # nube: solo las que si mueren, con su eje real
    for mask, role in ((kills & safe, "v3"), (kills & ~safe, "mazF")):
        ax.scatter(margin[mask], delay[mask], s=1.8, alpha=0.30,
                   color=fs.c(role), linewidths=0, rasterized=True)
    ax.set_ylim(0, 115)
    ax.set_yticks([0, 50, 100])
    ax.set_xlabel("Margen bajo el umbral letal (nM)")
    ax.set_ylabel("Retardo hasta\nla muerte (min)")

    # banda de censura: mismo eje x, sin eje y propio. Los puntos se dispersan
    # en vertical solo para que no se apilen en una linea.
    rng = np.random.default_rng(0)
    cen = ~kills
    axc.scatter(margin[cen], rng.uniform(0, 1, cen.sum()), s=1.8, alpha=0.30,
                color=fs.c("muted"), linewidths=0, rasterized=True)
    axc.set_ylim(-0.2, 1.2)
    axc.set_yticks([])
    axc.tick_params(labelbottom=False, length=0)
    for sp in ("left", "right", "bottom"):
        axc.spines[sp].set_visible(False)
    axc.spines["top"].set_visible(False)
    # Rotulo directo DENTRO de la banda, sobre su propio grupo: colgado a la
    # derecha obligaba a reservar margen y dejaba una franja vacia.
    axc.text(0.015, 0.5, "Nunca alcanza nivel letal (6 h)",
             transform=axc.transAxes, va="center", ha="left", zorder=6)

    # marca de corte: el eje y esta partido entre la nube y la banda
    for yf in (0.800, 0.818):
        fig.add_artist(plt.Line2D([0.081, 0.095], [yf, yf + 0.008],
                                  color="black", lw=0.6))
    return fig


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="usa la cache del muestreo reducido")
    args = ap.parse_args()

    from model import figstyle as fs

    prob, results, Y = load(args.fast)
    fs.apply_style()

    written = [fs.save_panel(panel_sobol(prob, results), "fig3_sobol"),
               fs.save_panel(panel_tradeoff(Y), "fig3_tradeoff")]

    print(f"  {len(written)} paneles escritos en figures/panels/")
    for p in written:
        print(f"    {p.name}")

    safe, kills = Y[:, 1] > 0, Y[:, 2] < CENSOR_MIN
    print()
    print(f"    cumple ambas          {(safe & kills).mean():>6.1%}")
    print(f"    segura, no contiene   {(safe & ~kills).mean():>6.1%}")
    print(f"    contiene, dispara     {(~safe & kills).mean():>6.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
