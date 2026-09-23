#!/usr/bin/env python3
"""Corre el analisis cinetico y escribe la tabla de resultados.

Responde tres preguntas que el diseno de secuencia no puede contestar:

  1. Es alcanzable la afirmacion de "<30 min" del proyecto?
  2. Cuanta diferencia hace la cromoproteina frente a la fluorescente?
  3. El kill switch contiene la celula, y sobrevive al defecto de RBS?

Las figuras del manuscrito se generan aparte, con `model/fig1_sensor.py`,
`model/fig2_killswitch.py` y `model/compose_figures.py`.

Uso:  .venv/bin/python model/run_analysis.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import kinetics as kn  # noqa: E402
from model.kinetics import KillSwitchConfig, SensorConfig  # noqa: E402

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


def main() -> int:
    print("=" * 66)
    print("Analisis cinetico -- Yakuponia")
    print("=" * 66)
    print("Parametros de literatura para sistemas analogos, no mediciones")
    print("de este proyecto. Sirven para ordenes de magnitud.")

    rows = q1_time_to_signal()
    q2_rbs_defect()
    q3_killswitch()

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
