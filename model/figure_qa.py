"""Auditorias obligatorias sobre las figuras exportadas.

Tres controles, ninguno opcional para una figura de manuscrito:

  1. alineacion de paneles  -- bordes, anchos y canales iguales dentro de
     1.5 pt en toda figura multipanel
  2. piso tipografico       -- ningun glifo renderizado por debajo de 5 pt
  3. colisiones             -- ningun texto encima de otro texto ni de un
     trazo en el PDF final

Los auditores viven en la skill `nature-figure`; este modulo solo los localiza
y traduce su salida. Si la skill no esta instalada, el export sigue adelante y
avisa que la figura NO quedo auditada: el contrato prohibe afirmar que paso un
control que nunca corrio.
"""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

# El piso de 5 pt vale para todo glifo, incluidos super/subindices.
MIN_PT = 5.0
TOL_PT = 1.5

_SKILL_ENV = "NATURE_FIGURE_SCRIPTS"
_SKILL_DEFAULT = Path.home() / ".claude/skills/nature-figure/scripts"


def _scripts_dir() -> Path | None:
    raw = os.environ.get(_SKILL_ENV)
    cand = Path(raw).expanduser() if raw else _SKILL_DEFAULT
    return cand if (cand / "audit_panel_alignment.py").exists() else None


def _load(name: str):
    """Importa un auditor de la skill por ruta, sin tocar sys.path global."""
    d = _scripts_dir()
    if d is None:
        return None
    spec = importlib.util.spec_from_file_location(name, d / f"{name}.py")
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(name, mod)
    spec.loader.exec_module(mod)
    return mod


def _warn(msg: str) -> None:
    print(f"  [QA] NO AUDITADO: {msg}", file=sys.stderr)


def audit_alignment(fig, base: Path) -> None:
    """Mide la geometria final de los paneles. Bloquea si estan desalineados.

    Una figura de un solo panel registra NOT APPLICABLE y pasa.
    """
    mod = _load("audit_panel_alignment")
    if mod is None:
        _warn("alineacion de paneles (skill nature-figure no encontrada)")
        return
    mod.require_matplotlib_panel_alignment(
        fig,
        json_out=str(base.with_name(base.name + ".alignment.json")),
        tolerance_pt=TOL_PT,
        gutter_tolerance_pt=TOL_PT,
        strict=True,
    )


def audit_pdf(pdf: Path) -> None:
    """Corre el auditor tipografico y el de colisiones sobre el PDF final."""
    d = _scripts_dir()
    if d is None:
        _warn(f"texto y colisiones de {pdf.name}")
        return

    text = subprocess.run(
        [sys.executable, str(d / "audit_pdf_text.py"), str(pdf),
         "--min-pt", str(MIN_PT)],
        capture_output=True, text=True)
    if text.returncode:
        raise SystemExit(
            f"[QA] {pdf.name}: glifos por debajo de {MIN_PT} pt\n"
            f"{text.stdout}{text.stderr}")

    coll = subprocess.run(
        [sys.executable, str(d / "audit_figure_collisions.py"), str(pdf),
         "--json-out", str(pdf.with_suffix(".collision-audit.json"))],
        capture_output=True, text=True)
    out = coll.stdout + coll.stderr
    if coll.returncode == 1:
        raise SystemExit(f"[QA] {pdf.name}: colisiones FAIL\n{out}")
    if coll.returncode == 2:
        _warn(f"colisiones de {pdf.name} (dependencia o PDF ilegible)\n{out}")
        return
    # Se cuentan solo los hallazgos marcados, no cualquier linea que mencione
    # WARN: la nota al pie del auditor contiene esa palabra y inflaba la cuenta.
    warn = [l for l in out.splitlines() if l.lstrip().startswith("[WARN]")]
    if warn:
        print(f"  [QA] {pdf.name}: {len(warn)} WARN a revisar a tamano final")
    else:
        print(f"  [QA] {pdf.name}: sin colisiones")
