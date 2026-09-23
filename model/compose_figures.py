#!/usr/bin/env python3
"""Compone los paneles SVG en figuras de manuscrito y las exporta.

Une las dos mitades del flujo: los generadores (`fig2_killswitch.py` y
companeros) dibujan paneles sueltos a SVG, y este script los monta en un
lienzo del ancho exacto de la revista, les pone las etiquetas de panel y
exporta PDF y PNG.

Dos detalles que no son obvios y que estaban detras de casi todos los fallos:

  1. Un panel de matplotlib viene dimensionado en **pt**; el lienzo del
     compositor lleva un viewBox en **mm**. Al insertarlo, cada numero del
     panel se lee como mm, asi que hace falta escalar por 25.4/72 ademas de
     ajustar al ancho de columna. Sin eso los paneles salen 2.83x mas grandes.
  2. El compositor dibuja su etiqueta de panel con `size=12` en unidades del
     lienzo, que en un viewBox en mm son ~34 pt. Las etiquetas se inyectan
     aqui, convertidas a pt reales.

Sobre la QA del PDF exportado: los auditores `audit_pdf_text.py` y
`audit_figure_collisions.py` de la skill `nature-figure` no sirven sobre estos
PDF. cairosvg emite `Tf 1` y aplica el tamano por matriz de texto, asi que el
primero lee 1 pt en todo; y el segundo agrupa corridas contiguas en una sola
caja, de modo que dos ticks separados del eje aparecen como solapados. La
puerta buena aqui es `validate_fonts.py` sobre el SVG compuesto, que si
acumula la escala de cada ancestro. Para comprobar el PDF final conviene medir
los spans con PyMuPDF en vez de confiar en esos dos informes.

Uso:  .venv/bin/python model/compose_figures.py
"""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
PANELS = FIG / "panels"
SKILL = Path.home() / ".claude/skills/scientific-figure/scripts"

MM_PER_PT = 25.4 / 72
COL_DOUBLE_MM = 183.0
LABEL_PT = 8.0          # minimo de Nature para etiqueta de panel


def panel_size_pt(path: Path) -> tuple[float, float]:
    """Ancho y alto declarados por el SVG, en pt."""
    head = path.read_text(encoding="utf-8")[:900]
    w = re.search(r'width="([\d.]+)pt"', head)
    h = re.search(r'height="([\d.]+)pt"', head)
    if not (w and h):
        raise ValueError(f"{path.name}: no declara tamano en pt")
    return float(w.group(1)), float(h.group(1))


def run(cmd: list[str]) -> str:
    """Corre un script de la skill con uv, que resuelve sus dependencias."""
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"fallo: {' '.join(cmd[-3:])}\n{r.stdout}{r.stderr}")
    return r.stdout


def row_layout(names: list[str], *, x0: float = 6.0, gut: float = 4.0,
               y0: float = 4.0, right: float = 3.0) -> tuple:
    """Coloca n paneles en una fila, repartidos por el ancho de columna.

    Devuelve (panels, labels, alto). Todos comparten escala, que se calcula
    del primero: la conversion pt->mm va incluida, porque el lienzo lleva un
    viewBox en mm y un numero del panel entraria como mm.
    """
    sizes = [panel_size_pt(PANELS / f"{n}.svg") for n in names]
    total_pt = sum(w for w, _ in sizes)
    avail = COL_DOUBLE_MM - x0 - right - gut * (len(names) - 1)
    scale = round(avail / total_pt, 4)

    panels, labels, x = [], [], x0
    for i, (n, (w_pt, _)) in enumerate(zip(names, sizes)):
        panels.append({"src": f"figures/panels/{n}.svg",
                       "x_mm": round(x, 2), "y_mm": round(y0, 2),
                       "scale": scale})
        labels.append({"t": "abcdefgh"[i],
                       "x": round(x - 4.2, 2), "y": round(y0 + 3.4, 2)})
        x += w_pt * scale + gut
    height = max(h for _, h in sizes) * scale
    return panels, labels, height


def build_fig1() -> dict:
    """Fig. 1 — curso temporal y dosis-respuesta, lado a lado."""
    y0 = 4.0
    panels, labels, h = row_layout(["fig1_timecourse", "fig1_dose"], y0=y0)
    return {"name": "fig1_sensor",
            "config": {"width_mm": COL_DOUBLE_MM,
                       "height_mm": round(y0 + h + 3.0, 1),
                       "journal": "nature", "panels": panels},
            "labels": labels}


def build_fig3() -> dict:
    """Fig. 3 — Sobol arriba, tension debajo; ambos al mismo ancho.

    Los dos paneles se dibujan al mismo ancho y comparten escala, asi que sus
    ejes quedan alineados sin calcularlo aparte.
    """
    x0, y0, vg, right = 6.0, 4.0, 3.5, 3.0
    w_sob, h_sob = panel_size_pt(PANELS / "fig3_sobol.svg")
    w_tra, h_tra = panel_size_pt(PANELS / "fig3_tradeoff.svg")

    scale = round((COL_DOUBLE_MM - x0 - right) / max(w_sob, w_tra), 4)
    y_tra = round(y0 + h_sob * scale + vg, 2)

    panels, labels = [], []
    for name, y in (("fig3_sobol", y0), ("fig3_tradeoff", y_tra)):
        panels.append({"src": f"figures/panels/{name}.svg",
                       "x_mm": round(x0, 2), "y_mm": round(y, 2),
                       "scale": scale})
        labels.append({"t": "ab"[len(labels)],
                       "x": round(x0 - 4.2, 2), "y": round(y + 3.4, 2)})

    return {"name": "fig3_sensitivity",
            "config": {"width_mm": COL_DOUBLE_MM,
                       "height_mm": round(y_tra + h_tra * scale + 3.0, 1),
                       "journal": "nature", "panels": panels},
            "labels": labels}


def build_fig2() -> dict:
    """Fig. 2 — rejilla 3x2 de version x destino del plasmido."""
    x0, gut, y0, vg = 6.0, 3.0, 4.0, 2.5
    pw_pt, ph_pt = panel_size_pt(PANELS / "fig2_v1_kept.svg")
    avail = (COL_DOUBLE_MM - x0 - gut - 3.0) / 2
    scale = round(avail / pw_pt, 4)     # incluye la conversion pt->mm
    pw, ph = pw_pt * scale, ph_pt * scale

    panels, labels = [], []
    for r, ver in enumerate(("v1", "v2", "v3")):
        for col, fate in enumerate(("kept", "lost")):
            x = round(x0 + col * (pw + gut), 2)
            y = round(y0 + r * (ph + vg), 2)
            panels.append({"src": f"figures/panels/fig2_{ver}_{fate}.svg",
                           "x_mm": x, "y_mm": y, "scale": scale})
            # La columna derecha no dibuja eje Y, asi que su area util
            # arranca casi en el borde del lienzo del panel.
            labels.append({"t": "abcdef"[r * 2 + col],
                           "x": round(x - 4.2 if col == 0 else x + 1.0, 2),
                           "y": round(y + 3.4, 2)})

    # Sin panel de leyenda: el pie de figura ya declara el mapeo de color, y
    # el bloque dejaba una franja vacia bajo la rejilla.
    height = round(y0 + 3 * (ph + vg) - vg + 3.0, 1)

    return {"name": "fig2_killswitch",
            "config": {"width_mm": COL_DOUBLE_MM, "height_mm": height,
                       "journal": "nature", "panels": panels},
            "labels": labels}


def inject_labels(svg: Path, labels: list[dict]) -> None:
    """Anade las etiquetas de panel a 8 pt reales.

    El viewBox esta en mm, asi que el tamano se convierte de pt a unidades de
    usuario. Hacerlo aqui, y no en el compositor, es lo que mantiene la
    etiqueta en 8 pt en vez de los ~34 pt que saldrian por defecto.
    """
    size = round(LABEL_PT * MM_PER_PT, 3)
    block = ['  <g id="panel-labels">']
    for lab in labels:
        block.append(
            f'    <text x="{lab["x"]}" y="{lab["y"]}" font-size="{size}" '
            f'font-family="Liberation Sans" font-weight="bold" '
            f'fill="#000000">{lab["t"]}</text>')
    block.append("  </g>")
    svg.write_text(svg.read_text(encoding="utf-8")
                   .replace("</svg>", "\n".join(block) + "\n</svg>"),
                   encoding="utf-8")


def compose(spec: dict) -> Path:
    name = spec["name"]
    cfg_path = FIG / f"{name}.config.json"
    cfg_path.write_text(json.dumps(spec["config"], indent=2), encoding="utf-8")

    svg = FIG / f"{name}.svg"
    run([*("uv run --with svgutils --with lxml python".split()),
         str(SKILL / "compose.py"), str(cfg_path), "-o", str(svg)])
    inject_labels(svg, spec["labels"])

    # Puerta obligatoria: ningun glifo por debajo del minimo de la revista.
    out = run([*("uv run --with lxml python".split()),
               str(SKILL / "validate_fonts.py"), str(svg),
               "--journal", "nature"])
    report = json.loads(out)
    if report["issue_count"]:
        worst = min(i["effective_pt"] for i in report["issues"])
        raise SystemExit(
            f"{name}: {report['issue_count']} glifos bajo "
            f"{report['minimum_pt']} pt (el menor, {worst:.2f} pt). "
            f"No se exporta.")
    print(f"  {name}: {report['checked_count']} textos, todos "
          f">= {report['minimum_pt']} pt")

    for ext, extra in (("pdf", []), ("png", ["--dpi", "300"])):
        run([*("uv run --with cairosvg --with lxml python".split()),
             str(SKILL / "export.py"), str(svg),
             "--out", str(FIG / f"{name}.{ext}"), *extra])
    print(f"  {name}: {spec['config']['width_mm']} x "
          f"{spec['config']['height_mm']} mm -> PDF + PNG")
    return svg


BUILDERS = [("fig1_timecourse", build_fig1, "model/fig1_sensor.py"),
            ("fig2_v1_kept", build_fig2, "model/fig2_killswitch.py"),
            ("fig3_sobol", build_fig3, "model/fig3_sensitivity.py")]


def main() -> int:
    missing = [gen for probe, _, gen in BUILDERS
               if not (PANELS / f"{probe}.svg").exists()]
    if missing:
        print("faltan paneles; corre primero:")
        for gen in missing:
            print(f"  .venv/bin/python {gen}")
        return 1
    for _, build, _ in BUILDERS:
        compose(build())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
