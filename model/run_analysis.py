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
    """Escribe las tres figuras del analisis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "figure.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})

    # --- fig 1: las tres versiones ---
    fig, ax = plt.subplots(figsize=(6.5, 4))
    colors = {"v1": "#2a9d3f", "v2": "#c0392b", "v3": "#7d3c98"}
    for ver, v in VERSIONS.items():
        cfg = SensorConfig(reporter=v["reporter"], readout="visual",
                           rbs_efficiency=eff_of(v), inducer=10.0)
        t_a, y = kn.simulate_sensor(cfg, t_end=9000.0)
        tt = kn.time_to_signal(cfg)
        suffix = f"{tt:.0f} min" if tt else "no alcanza"
        ax.plot(t_a / 60.0, y[2], color=colors[ver], lw=2.0 if ver == "v3" else 1.6,
                ls="-" if ver == "v3" else ":",
                label=f"{ver}: {v['note']} — {suffix}")
        if tt:
            ax.plot([tt], [kn.p("thr_visual")], "o", color=colors[ver],
                    ms=5, zorder=5)
    ax.axhline(kn.p("thr_visual"), ls="--", c="k", lw=1,
               label="umbral visible a simple vista")
    ax.axvspan(0, CLAIM_MIN, color="#e67e22", alpha=0.10)
    ax.axvline(CLAIM_MIN, c="#e67e22", lw=1.2)
    ax.text(CLAIM_MIN - 2, kn.p("thr_visual") * 1.45, "30 min\n(afirmado)",
            color="#e67e22", fontsize=8, ha="right")
    ax.set_xlabel("tiempo (min)")
    ax.set_ylabel("reportero maduro (nM)")
    ax.set_title("Solo la v3 da senal visible en B. subtilis")
    ax.set_xlim(0, 140)
    ax.set_ylim(0, kn.p("thr_visual") * 1.9)
    ax.legend(fontsize=6.5, loc="lower right")
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

    # --- fig 3: kill switch por version ---
    ks_cfgs = [("v1 sin ssrA", 5, 4, False), ("v2 con ssrA", 5, 4, True),
               ("v3 RBS B. subtilis", 7, 6, True)]
    fig, axes = plt.subplots(3, 2, figsize=(9, 8), sharex=True)
    for row, (label, sd_e, sd_f, ssrA) in enumerate(ks_cfgs):
        e_e = kn.rbs_efficiency_from_sd(sd_e)
        e_f = kn.rbs_efficiency_from_sd(sd_f)
        for col, lost in enumerate((None, 3600.0)):
            ax = axes[row, col]
            cfg = KillSwitchConfig(rbs_mazE=e_e, rbs_mazF=e_f, ssrA=ssrA,
                                   plasmid_lost_at=lost)
            t_a, y = kn.simulate_killswitch(cfg, t_end=21600.0)
            ax.plot(t_a / 60.0, y[0], color="#2471a3", lw=1.5,
                    label="MazE libre")
            ax.plot(t_a / 60.0, y[1], color="#c0392b", lw=1.5,
                    label="MazF libre")
            ax.axhline(kn.p("tox_lethal"), ls="--", c="k", lw=1)
            if lost:
                ax.axvline(60, c="#e67e22", lw=1.1)
            ax.set_yscale("log")
            ax.set_ylim(1e-2, 1e4)
            d = kn.time_to_death(cfg)
            good = (d is not None) == (lost is not None)
            verdict = "sobrevive" if d is None else f"muere {d:.0f} min"
            ax.set_title(f"{label} — {'perdido' if lost else 'retenido'}: "
                         f"{verdict}  {'OK' if good else 'FALLA'}",
                         fontsize=8, color="#1d6b2b" if good else "#a93226")
            if col == 0:
                ax.set_ylabel("nM")
    for ax in axes[2]:
        ax.set_xlabel("tiempo (min)")
    axes[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Solo la v3 discrimina: la v2 ni siquiera mata al perder "
                 "el plasmido", fontsize=10)
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
