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


# Cada version con su reportero, espaciado y apareamiento SD del RBS.
# eff = eficiencia traduccional resultante.
VERSIONS = {
    "v1": dict(reporter="gfp", spacer_ok=False, sd=4,
               note="GFP, sin espaciador, B0032"),
    "v2": dict(reporter="chromoprotein", spacer_ok=True, sd=4,
               note="eforRed, espaciador 7 nt, B0032"),
    "v3": dict(reporter="chromoprotein", spacer_ok=True, sd=6,
               note="eforRed, espaciador 7 nt, RBS B. subtilis"),
}


def eff_of(v: dict) -> float:
    """Eficiencia traduccional: apareamiento SD, penalizado si falta espaciador."""
    e = kn.rbs_efficiency_from_sd(v["sd"])
    return e if v["spacer_ok"] else e * 0.67


def q1_time_to_signal() -> list[tuple]:
    """Tiempo hasta senal para cada version, en chasis B. subtilis."""
    section("1. Tiempo hasta senal detectable (FtaGen, RBS medio)")
    rows = []
    for ver, v in VERSIONS.items():
        eff = eff_of(v)
        for readout in ("visual", "instrument"):
            cfg = SensorConfig(reporter=v["reporter"], readout=readout,
                               rbs_efficiency=eff, inducer=10.0)
            tt = kn.time_to_signal(cfg)
            ok = tt is not None and tt <= CLAIM_MIN
            rows.append((ver, readout, tt, ok))
            print(f"  {ver}  {v['note']:<38} {readout:<11} "
                  f"{fmt(tt):>12}   {'cumple <30 min' if ok else ''}")
    print()
    print("  La v2 corrigio el espaciador pero conservo un SD de 4 nt, que en")
    print("  B. subtilis no basta: el sensor seguia sin dar senal visible.")
    return rows


def q2_rbs_defect() -> None:
    """Fuerza del SD y su efecto sobre la traduccion."""
    section("2. Shine-Dalgarno en B. subtilis")
    print("  El 16S de B. subtilis exige mas apareamiento contiguo que el de")
    print("  E. coli. Los RBS del registro iGEM estan caracterizados en E. coli.")
    print()
    print("  RBS                      SD   apareamiento  eficiencia")
    for name, sd in (("BBa_B0032 (v1, v2)", 4), ("BBa_B0034 (v1, v2)", 5),
                     ("Bsub medium (v3)", 6), ("Bsub strong (v3)", 7)):
        e = kn.rbs_efficiency_from_sd(sd)
        print(f"  {name:<24} {sd:>2} nt        {e:>5.0%}")
    print()
    print("  efecto sobre el tiempo hasta senal (cromoproteina, visual):")
    for sd in (4, 5, 6, 7):
        cfg = SensorConfig(reporter="chromoprotein", readout="visual",
                           rbs_efficiency=kn.rbs_efficiency_from_sd(sd),
                           inducer=10.0)
        print(f"    SD {sd} nt  {fmt(kn.time_to_signal(cfg)):>12}")


def q3_killswitch() -> None:
    """Contencion en el chasis real, version por version."""
    section("3. Kill switch MazEF en B. subtilis")

    print("  condicion de viabilidad:  k_mazF < d_mazF * umbral_letal")
    for ssrA in (False, True):
        tag = "con tag ssrA" if ssrA else "sin tag"
        cota = kn.max_toxin_synthesis(ssrA)
        print(f"    MazF {tag:<14} cota {cota:.4f} nM/s   "
              f"viable={kn.is_viable(ssrA)}")

    print()
    print(f"  {'':34}{'retenido':<14}{'perdido a 60 min'}")
    print(f"  {'':34}{'(sobrevivir)':<14}{'(morir)'}")
    configs = [
        ("v1  B0032/B0034, sin ssrA", 5, 4, False),
        ("v2  B0032/B0034, ssrA", 5, 4, True),
        ("v3  RBS B. subtilis, ssrA", 7, 6, True),
    ]
    for label, sd_e, sd_f, ssrA in configs:
        e_e = kn.rbs_efficiency_from_sd(sd_e)
        e_f = kn.rbs_efficiency_from_sd(sd_f)
        keep = kn.time_to_death(
            KillSwitchConfig(rbs_mazE=e_e, rbs_mazF=e_f, ssrA=ssrA))
        lost = kn.time_to_death(
            KillSwitchConfig(rbs_mazE=e_e, rbs_mazF=e_f, ssrA=ssrA,
                             plasmid_lost_at=3600.0))
        k_s = "sobrevive" if keep is None else f"MUERE {keep:.0f}m"
        l_s = "sobrevive" if lost is None else f"muere {lost:.0f}m"
        ok = keep is None and lost is not None
        print(f"  {label:<34}{k_s:<14}{l_s:<18}"
              f"{'  <- correcto' if ok else '  <- FALLA'}")

    print()
    print("  La v2 no mata al perder el plasmido: con SD de 4-5 nt la toxina")
    print("  nunca alcanza concentracion letal. La contencion parecia resuelta")
    print("  solo porque el modelo asumia traduccion plena.")

    import numpy as np
    e_e = kn.rbs_efficiency_from_sd(7)
    e_f = kn.rbs_efficiency_from_sd(6)
    keep_ss = kn.simulate_killswitch(
        KillSwitchConfig(rbs_mazE=e_e, rbs_mazF=e_f, ssrA=True))[1][1, -1]
    thr = kn.p("tox_lethal")
    t_a, y = kn.simulate_killswitch(
        KillSwitchConfig(rbs_mazE=e_e, rbs_mazF=e_f, ssrA=True,
                         plasmid_lost_at=3600.0), t_end=21600.0)
    hot = np.flatnonzero(y[1] >= thr)
    print()
    print("  margenes de la v3:")
    print(f"    en reposo MazF libre {keep_ss:.0f} nM contra umbral {thr:.0f} nM")
    if hot.size:
        print(f"    ventana letal tras perder el plasmido: "
              f"{(t_a[hot[-1]] - t_a[hot[0]]) / 60:.0f} min")


def make_figures() -> None:
    """Escribe las cuatro figuras del analisis cinetico."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    from model import figstyle as fs

    fs.apply_style()

    # ----------------------------------------------------------------
    # Fig 1 — tiempo hasta senal visible, por version
    # ----------------------------------------------------------------
    # La leyenda antes caia encima de las curvas v1/v2. Ahora cada curva se
    # rotula donde termina, que ademas elimina el viaje de ojo a la leyenda.
    fig, ax = plt.subplots(figsize=(fs.COL_SINGLE, 2.6))
    # `top` deja aire suficiente para que el titulo no toque la curva v3,
    # que llega hasta el borde superior del area de datos.
    fig.subplots_adjust(left=0.165, right=0.70, top=0.845, bottom=0.175)

    thr = kn.p("thr_visual")
    colors = {"v1": fs.C["v1"], "v2": fs.C["v2"], "v3": fs.C["v3"]}
    ends = []
    for ver, v in VERSIONS.items():
        cfg = SensorConfig(reporter=v["reporter"], readout="visual",
                           rbs_efficiency=eff_of(v), inducer=10.0)
        t_a, y = kn.simulate_sensor(cfg, t_end=9000.0)
        tt = kn.time_to_signal(cfg)
        hero = ver == "v3"
        # Se recorta la serie al rango visible en vez de dejar que matplotlib
        # la recorte: una curva que sale por el borde superior conserva una
        # caja que abarca toda la pagina y colisiona con cualquier titulo.
        mins = t_a / 60.0
        vis = mins <= 140
        ax.plot(mins[vis], y[2][vis], color=colors[ver],
                lw=1.5 if hero else 1.0, ls="-" if hero else (0, (3, 1.6)),
                zorder=4 if hero else 3)
        keep = vis
        ends.append((ver, (t_a / 60.0)[keep][-1], y[2][keep][-1], tt))
        if tt:
            ax.plot([tt], [thr], "o", color=colors[ver], ms=3.4, zorder=6,
                    markeredgecolor="white", markeredgewidth=0.5)

    # Rotulos directos a la derecha, anclados al punto real donde termina cada
    # curva dentro del eje. El ancla se limita al recorte visible para que la
    # guia no salga por arriba del panel.
    # el tope cubre el maximo visible: ninguna curva se recorta contra el
    # borde superior, que es lo que inflaba su caja a toda la pagina
    ymax = max(ye for _, _, ye, _ in ends) * 1.06
    # el hueco de 0.53 queda libre para el rotulo del umbral, que va al mismo
    # margen derecho a la altura de su propia linea
    # el rotulo del umbral ocupa ~0.374 (su propia altura), asi que v2 y v1
    # se reparten por debajo sin tocarlo
    slots = {"v3": 0.80, "v2": 0.225, "v1": 0.085}
    for ver, xe, ye, tt in ends:
        note_txt = f"{tt:.0f} min" if tt else "no alcanza"
        ax.annotate(f"{ver} · {note_txt}",
                    xy=(min(xe, 139), min(ye, ymax * 0.97)),
                    xytext=(146, ymax * slots[ver]),
                    fontsize=6.3, color=colors[ver], va="center", ha="left",
                    annotation_clip=False, zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.45,
                                    color=colors[ver], shrinkA=1.0,
                                    shrinkB=1.0))
    # pie del rotulo v3, anclado justo debajo de su slot
    ax.text(146, ymax * (slots["v3"] - 0.045),
            "eforRed, espaciador 7 nt\nRBS $\\it{B.\\ subtilis}$",
            fontsize=5.4, color=fs.C["mid"], va="top", ha="left",
            linespacing=1.4, clip_on=False)

    # El umbral se rotula dentro del panel, sobre la zona vacia de la
    # izquierda: el margen derecho ya lo ocupan los tres rotulos de version.
    ax.axhline(thr, ls="--", color=fs.C["ink"], lw=0.7, zorder=2)
    # a la izquierda de donde la v3 cruza el umbral, con fondo por si alguna
    # raya de la propia linea queda debajo
    # todos los rotulos viven en el margen derecho, fuera del area de datos:
    # dentro, la caja de cualquier curva abarca el panel entero y colisiona
    ax.text(146, thr, "umbral visible\na simple vista", fontsize=5.6,
            color=fs.C["ink"], va="center", ha="left", linespacing=1.35,
            zorder=6, clip_on=False)
    ax.axvspan(0, CLAIM_MIN, color=fs.C["warn"], alpha=0.085, zorder=0)
    ax.axvline(CLAIM_MIN, color=fs.C["warn"], lw=0.8, zorder=2)
    ax.text(CLAIM_MIN - 2.5, ymax * 0.975, "30 min\nafirmado", color=fs.C["warn"],
            fontsize=5.8, ha="right", va="top", linespacing=1.25)
    ax.set_xlabel("tiempo (min)")
    ax.set_ylabel("reportero maduro (nM)")
    ax.set_xlim(0, 140)
    ax.set_ylim(0, ymax)
    # El titulo va como texto de figura, por encima del area de ejes: como
    # titulo de eje su caja se solapaba con la de la curva v3, que llega al
    # borde superior del panel.
    fig.text(0.035, 0.945, "Solo la v3 da senal visible en "
             "$\\it{B.\\ subtilis}$", fontsize=7.5, va="top", ha="left",
             color=fs.C["ink"])
    fs.save(fig, "01_time_to_signal")

    # ----------------------------------------------------------------
    # Fig 2 — dosis-respuesta a los 30 min
    # ----------------------------------------------------------------
    # Los dos umbrales se rotulan sobre su propia linea en vez de competir
    # por un hueco en la leyenda con los tres reporteros.
    fig, ax = plt.subplots(figsize=(fs.COL_SINGLE, 2.6))
    fig.subplots_adjust(left=0.175, right=0.635, top=0.845, bottom=0.175)

    doses = np.logspace(-2, 2, 60)
    styles = {"gfp": ("GFP (v1)", fs.C_REPORTER["gfp"]),
              "mcherry": ("mCherry (v1)", fs.C_REPORTER["mcherry"]),
              "chromoprotein": ("cromoproteina (v2, v3)",
                                fs.C_REPORTER["chromoprotein"])}
    last = {}
    for rep, (label, color) in styles.items():
        r = kn.dose_response(doses, SensorConfig(reporter=rep), t_read=1800.0)
        hero = rep == "chromoprotein"
        ax.plot(doses, r, color=color, lw=1.5 if hero else 1.0,
                zorder=4 if hero else 3)
        last[rep] = (label, color, r[-1])

    # Las tres curvas y los dos umbrales llegan planos al borde derecho, asi
    # que cada rotulo se escribe a la altura de su propia serie y no hace
    # falta ninguna guia: sin lineas, no hay cruces que resolver.
    for val, ls in ((kn.p("thr_visual"), "--"),
                    (kn.p("thr_instrument"), ":")):
        ax.axhline(val, ls=ls, color=fs.C["ink"], lw=0.7, zorder=2)

    right = doses[-1] * 1.35
    marks = [(last["chromoprotein"][2], last["chromoprotein"][0],
              fs.C_REPORTER["chromoprotein"], 6.3),
             (last["gfp"][2], last["gfp"][0], fs.C_REPORTER["gfp"], 6.3),
             (last["mcherry"][2], last["mcherry"][0],
              fs.C_REPORTER["mcherry"], 6.3),
             (kn.p("thr_visual"), "umbral visual", fs.C["ink"], 5.8),
             (kn.p("thr_instrument"), "umbral instrumento", fs.C["ink"], 5.8)]
    for yv, txt, color, size in marks:
        ax.text(right, yv, txt, fontsize=size, color=color, va="center",
                ha="left", zorder=6, clip_on=False)

    ax.set_xscale("log")
    ax.set_yscale("log")
    fs.log_decades(ax, "y")
    # los exponentes de la notacion cientifica se dibujan al 70% del tamano
    # del tick: con 6.5 pt caerian a 4.55 pt, por debajo del piso de 5 pt
    ax.tick_params(axis="both", labelsize=7.4)
    ax.set_xlabel("contaminante (µM)")
    ax.set_ylabel("reportero maduro a los 30 min (nM)")
    ax.set_xlim(doses[0], doses[-1])
    # titulo como texto de figura: ver nota en la fig 1
    fig.text(0.035, 0.945, "Ningun reportero cruza el umbral visual "
             "a los 30 min", fontsize=7.5, va="top", ha="left",
             color=fs.C["ink"])
    fs.save(fig, "02_dose_response")

    # ----------------------------------------------------------------
    # Fig 3 — kill switch: version x destino del plasmido
    # ----------------------------------------------------------------
    # Matriz 3x2. La condicion (columna) y la version (fila) se rotulan una
    # sola vez en los margenes; cada panel conserva solo su veredicto, como
    # insignia en vez de titulo largo. Asi el titulo deja de competir con los
    # datos y las seis celdas se leen como una tabla.
    ks_cfgs = [("v1", "B0032/B0034\nsin ssrA", 5, 4, False),
               ("v2", "B0032/B0034\nssrA", 5, 4, True),
               ("v3", "RBS $\\it{B.\\ subtilis}$\nssrA", 7, 6, True)]
    cols = [(None, "plasmido retenido"), (3600.0, "plasmido perdido a 60 min")]

    fig = plt.figure(figsize=(fs.COL_DOUBLE, 4.6))
    gs = GridSpec(3, 2, figure=fig, hspace=0.34, wspace=0.10,
                  left=0.165, right=0.795, top=0.825, bottom=0.105)
    axes = np.empty((3, 2), dtype=object)
    lethal = kn.p("tox_lethal")

    for row, (ver, note_txt, sd_e, sd_f, ssrA) in enumerate(ks_cfgs):
        e_e = kn.rbs_efficiency_from_sd(sd_e)
        e_f = kn.rbs_efficiency_from_sd(sd_f)
        for col, (lost, _) in enumerate(cols):
            ax = fig.add_subplot(gs[row, col])
            axes[row, col] = ax
            cfg = KillSwitchConfig(rbs_mazE=e_e, rbs_mazF=e_f, ssrA=ssrA,
                                   plasmid_lost_at=lost)
            t_a, y = kn.simulate_killswitch(cfg, t_end=21600.0)
            # Se enmascara lo que cae por debajo del piso visible en vez de
            # dejar que matplotlib lo recorte: una curva recortada conserva
            # una caja que llega al borde y pisa los rotulos del eje.
            mins = t_a / 60.0
            floor = 6e-3
            for series, color, lab, z in ((y[0], fs.C["mazE"], "MazE libre", 3),
                                          (y[1], fs.C["mazF"], "MazF libre", 4)):
                vis = np.where(series >= floor, series, np.nan)
                ax.plot(mins, vis, color=color, lw=1.0, label=lab, zorder=z)
            ax.axhline(lethal, ls="--", color=fs.C["ink"], lw=0.6, zorder=2)
            if lost:
                ax.axvline(60, color=fs.C["warn"], lw=0.8, zorder=2)

            ax.set_yscale("log")
            # el piso baja un poco por debajo del ultimo tick: asi ninguna
            # curva toca el borde inferior ni invade los rotulos del eje x
            ax.set_ylim(4e-3, 1e4)
            ax.set_xlim(0, 360)
            ax.set_yticks([1e-2, 1e0, 1e2, 1e4])
            # sin el 0 de la columna derecha: chocaba con el 360 de la
            # izquierda al compartir el canal entre paneles
            ax.set_xticks([60, 120, 180, 240, 300, 360] if col
                          else [0, 60, 120, 180, 240, 300, 360])
            # ver nota de la fig 2: el exponente baja al 70% del tick
            ax.tick_params(axis="both", labelsize=7.4)

            d = kn.time_to_death(cfg)
            good = (d is not None) == (lost is not None)
            verdict = "sobrevive" if d is None else f"muere {d:.0f} min"
            # el veredicto va encima del panel: dentro lo cruzan las curvas
            ax.text(1.0, 1.015, verdict, transform=ax.transAxes,
                    fontsize=6.2, ha="right", va="bottom", zorder=7,
                    color=fs.C["ok"] if good else fs.C["fail"])

            if col != 0:
                ax.tick_params(labelleft=False)
            if row < 2:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("tiempo (min)")

        # Rotulo de fila en el margen izquierdo de la figura, no del eje: a la
        # altura del panel pero fuera de el, para no pisar los ticks.
        pos = axes[row, 0].get_position()
        y_mid = pos.y0 + pos.height / 2
        fig.text(0.008, y_mid + 0.052, ver, fontsize=8, fontweight="bold",
                 ha="left", va="center", color=fs.C["ink"])
        fig.text(0.008, y_mid - 0.038, note_txt, fontsize=5.3, ha="left",
                 va="center", color=fs.C["mid"], linespacing=1.45)

    # encabezado de columna, una sola vez, por encima de los veredictos
    for col, (_, title) in enumerate(cols):
        axes[0, col].set_title(title, fontsize=7.0, pad=21, color=fs.C["ink"])

    # un solo rotulo de eje y para las tres filas
    left_edge = axes[0, 0].get_position().x0
    fig.text(left_edge - 0.065, 0.475, "MazE / MazF libre (nM)", fontsize=7,
             rotation=90, va="center", ha="center", color=fs.C["ink"])

    # Etiqueta por encima del area de datos: dentro del panel, la caja de
    # cualquier curva la cruzaria, y pegada al eje choca con el tick 10^4.
    # El desplazamiento negativo solo cabe en la columna izquierda; en la
    # derecha la etiqueta se apoya en el borde del propio panel.
    for idx, ax in enumerate(axes.ravel()):
        dx = -20.0 if idx % 2 == 0 else -6.0
        fs.panel_label(ax, "abcdef"[idx], dx_pt=dx, dy_pt=9.0)

    handles = [
        plt.Line2D([], [], color=fs.C["mazE"], lw=1.0, label="MazE libre"),
        plt.Line2D([], [], color=fs.C["mazF"], lw=1.0, label="MazF libre"),
        plt.Line2D([], [], color=fs.C["ink"], lw=0.6, ls="--",
                   label=f"umbral letal ({lethal:.0f} nM)"),
        plt.Line2D([], [], color=fs.C["warn"], lw=0.8,
                   label="perdida del plasmido"),
    ]
    fig.legend(handles=handles, loc="center left",
               bbox_to_anchor=(0.815, 0.62), fontsize=6.3,
               handlelength=1.5, labelspacing=0.55)
    fig.text(0.815, 0.35, "Solo la v3 separa los dos\ndestinos: sobrevive con "
             "el\nplasmido y muere al perderlo.", fontsize=6.0,
             color=fs.C["mid"], va="top", linespacing=1.5)

    fig.suptitle("Solo la v3 discrimina: la v2 no mata ni al perder el "
                 "plasmido", fontsize=8.5, y=0.945)
    fs.save(fig, "03_killswitch")

    # ----------------------------------------------------------------
    # Fig 4 — barrido de sintesis de toxina
    # ----------------------------------------------------------------
    # Tres verticales competian sin rotulo propio y el texto verde cruzaba su
    # propia linea. Ahora cada cota se rotula arriba, fuera del area de datos,
    # y las anotaciones salen del eje hacia el margen derecho.
    fig, ax = plt.subplots(figsize=(fs.COL_SINGLE, 2.75))
    fig.subplots_adjust(left=0.165, right=0.985, top=0.695, bottom=0.165)

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

    ax.plot(ks, times, color=fs.C["mazF"], lw=1.5, zorder=4)
    ax.axvspan(ks[0], k_max, color=fs.C["ok"], alpha=0.075, zorder=0)
    ax.set_xscale("log")
    # ver nota de la fig 2: el exponente se dibuja al 70% del tick
    ax.tick_params(axis="both", labelsize=7.4)
    ax.set_xlim(ks[0], ks[-1])
    top = max(t for t in times if not np.isnan(t)) * 1.06
    ax.set_ylim(0, top)

    # Las cotas v2 y el diseno caen a 0.144 y 0.10: demasiado juntas para dos
    # rotulos centrados. Se escalonan en altura y se alinean por su lado.
    marks = [
        (k_max, fs.C["ink"], "--", f"cota v1\n{k_max:.4f}", 1.035, "center"),
        (kn.p("k_mazF"), fs.C["warn"], ":", f"diseno\n{kn.p('k_mazF'):.2f}",
         1.035, "center"),
        (k_max_v2, fs.C["ok"], "--", f"cota v2 (ssrA)\n{k_max_v2:.3f}",
         1.175, "right"),
    ]
    for xv, color, ls, txt, yf, ha in marks:
        ax.axvline(xv, ls=ls, color=color, lw=0.9, zorder=3)
        # rotulo por encima del area de datos: no cruza ni curva ni linea
        ax.text(xv, top * yf, txt, fontsize=5.8, color=color, ha=ha,
                va="bottom", linespacing=1.3, zorder=6, clip_on=False)

    # Ambos textos van en la banda vacia de la izquierda, dentro de la zona
    # viable: a la derecha cruzaban la curva y las tres verticales.
    ax.text(ks[0] * 1.2, top * 0.085, "celula viable", fontsize=6.0,
            color=fs.C["ok"], zorder=5)
    ax.text(ks[0] * 1.2, top * 0.30,
            "con ssrA la cota sube 22x\ny el diseno queda\ndel lado viable",
            fontsize=5.9, color=fs.C["ok"], ha="left", va="center",
            linespacing=1.45, zorder=6)
    # el subindice de $k_F$ se dibuja al 70%: a 7 pt cae a 4.9 pt, bajo el
    # piso de 5 pt, asi que la etiqueta sube a 7.5 pt
    ax.set_xlabel("sintesis de MazF, $k_F$ (nM/s)", fontsize=7.5)
    ax.set_ylabel("tiempo hasta muerte espuria (min)")
    # titulo como texto de figura: ver nota en la fig 1
    fig.text(0.035, 0.955, "El tag ssrA mueve la cota 22x y deja el diseno "
             "del lado viable", fontsize=7.5, va="top", ha="left",
             color=fs.C["ink"])
    fs.save(fig, "04_toxin_sweep")

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
