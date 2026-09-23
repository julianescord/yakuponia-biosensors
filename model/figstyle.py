"""Estilo compartido, derivado de `figures/theme.json`.

El tema es la unica fuente de color y tipografia del proyecto: una version es
el mismo color en las tres figuras, y hay exactamente dos tamanos de texto.
Este modulo no decide nada de eso, solo lo traduce a rcParams de matplotlib.

Contrato de figura (destino: manuscrito, Nature):

  - los paneles se dibujan a SVG con fondo transparente; el tamano fisico
    final lo fija el compositor en mm, no el panel
  - piso de 5 pt para cuerpo y 8 pt para etiqueta de panel
  - texto editable (`svg.fonttype='none'`, `pdf.fonttype=42`)
  - sin titulo dentro de la imagen: eso va al pie de figura

Que NO va dentro de una figura (va al pie, ver `figures/legends.md`):
titulo, frases explicativas, interpretacion. Dentro solo queda lo que no se
puede decir en el pie: unidades, categorias, umbrales y etiquetas de serie.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

ROOT = Path(__file__).resolve().parent.parent
THEME_PATH = ROOT / "figures" / "theme.json"
PANEL_DIR = ROOT / "figures" / "panels"
FIG_DIR = ROOT / "figures"

MM = 1.0 / 25.4
COL_SINGLE = 89 * MM
COL_DOUBLE = 183 * MM


@lru_cache(maxsize=1)
def theme() -> dict:
    """Carga `figures/theme.json`. Es la fuente unica de estilo."""
    with THEME_PATH.open(encoding="utf-8") as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def colors() -> dict[str, str]:
    """Roles semanticos: un color = un significado en todo el proyecto."""
    return {k: v for k, v in theme()["semantic_roles"].items()
            if not k.startswith("_")}


def c(role: str) -> str:
    """Color de un rol semantico. Falla fuerte si el rol no existe."""
    try:
        return colors()[role]
    except KeyError:
        raise KeyError(
            f"rol '{role}' no esta en semantic_roles de {THEME_PATH.name}; "
            f"disponibles: {sorted(colors())}") from None


def body_pt() -> float:
    return float(theme()["type_scale"]["body_pt"])


def label_pt() -> float:
    return float(theme()["type_scale"]["panel_label_pt"])


def apply_style() -> None:
    """rcParams derivados del tema. Llamar antes de crear figuras."""
    t = theme()
    fam = t["typography"]["family"]
    body = body_pt()

    mpl.rcParams.update({
        # texto editable: sin esto el pie no se puede corregir ni validar
        "font.family": "sans-serif",
        "font.sans-serif": [fam, "Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,

        # Dos tamanos, no mas. Todo el cuerpo comparte el mismo; la etiqueta
        # de panel es el unico texto distinto y se pone a mano.
        "font.size": body,
        "axes.titlesize": body,
        "axes.labelsize": body,
        "xtick.labelsize": body,
        "ytick.labelsize": body,
        "legend.fontsize": body,
        "figure.titlesize": body,

        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.45,
        "ytick.minor.width": 0.45,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.minor.size": 1.4,
        "ytick.minor.size": 1.4,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.1,
        "lines.markersize": 3.0,
        "patch.linewidth": 0.6,

        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.25,

        "figure.dpi": 150,
        "savefig.dpi": 600,
        "axes.labelcolor": "black",
        "text.color": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.edgecolor": "black",
    })


def panel_label(ax, label: str, dx_pt: float = -13.0,
                dy_pt: float = 3.0) -> None:
    """Etiqueta de panel con desplazamiento fisico, no fraccional.

    Un offset en fraccion de ejes da distancias distintas segun la altura del
    panel; en puntos queda igual en todos.
    """
    off = ScaledTranslation(dx_pt / 72, dy_pt / 72, ax.figure.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + off,
            fontsize=label_pt(), fontweight="bold", color="black",
            ha="left", va="bottom")


def log_decades(ax, axis: str = "y") -> None:
    """Una etiqueta por decada; evita el apilado de minor ticks."""
    from matplotlib.ticker import LogFormatterSciNotation, LogLocator
    t = ax.yaxis if axis == "y" else ax.xaxis
    t.set_major_locator(LogLocator(base=10.0))
    t.set_major_formatter(LogFormatterSciNotation(base=10.0))
    t.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10)),
                                   numticks=12))
    t.set_minor_formatter(mpl.ticker.NullFormatter())


def save_panel(fig, name: str, *, tight: bool = False) -> Path:
    """Guarda un panel como SVG para que lo componga `scientific-figure`.

    SVG conserva el texto como elementos inspeccionables y el fondo
    transparente permite montarlo sobre el lienzo.

    `tight` queda desactivado a proposito para los paneles de una rejilla:
    recorta segun el contenido, asi que un panel con eje Y sale mas ancho que
    uno sin el y las columnas dejan de alinear. Con el area de ejes fija en
    coordenadas absolutas, el lienzo sale identico en todos. Se activa solo
    para piezas sueltas, como un panel-leyenda.
    """
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    out = PANEL_DIR / f"{name}.svg"
    kw = dict(bbox_inches="tight", pad_inches=0.01) if tight else {}
    fig.savefig(out, transparent=True, **kw)
    plt.close(fig)
    return out
