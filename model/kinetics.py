"""Modelo cinetico de los biosensores y del kill switch de Yakuponia.

Implementa tres modelos de ODEs que responden a preguntas que el diseno de
secuencia por si solo no puede contestar:

  1. sensor    -- cuanto tarda un biosensor de represion en dar senal visible,
                  y si la afirmacion de "<30 min" del proyecto es alcanzable
  2. reporter  -- cuanta diferencia hace usar una cromoproteina en vez de una
                  proteina fluorescente para lectura a simple vista
  3. killswitch-- si el ratio antitoxina:toxina del diseno contiene la celula,
                  y que pasa cuando el defecto de espaciado RBS lo altera

Todos los parametros estan en PARAMS con su fuente y su incertidumbre. No son
mediciones de este proyecto: son valores de literatura para sistemas
analogos, usados para obtener ordenes de magnitud, no predicciones exactas.

Uso:
    from model import kinetics
    t, y = kinetics.simulate_sensor(inducer=1.0)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

# --------------------------------------------------------------------------
# Parametros
# --------------------------------------------------------------------------
# Cada entrada: (valor, unidad, fuente). Los ordenes de magnitud provienen de
# sistemas caracterizados en E. coli / B. subtilis; ver README del modelo.

PARAMS: dict[str, tuple[float, str, str]] = {
    # --- transcripcion / traduccion ---
    "k_txn_max": (0.5, "nM/s", "promotor sigma-A fuerte desreprimido"),
    "leak": (0.02, "-", "fraccion de expresion basal del promotor reprimido"),
    "k_tln": (0.05, "1/s", "tasa de traduccion por mRNA"),
    "d_mrna": (np.log(2) / 150, "1/s", "vida media mRNA bacteriano ~2.5 min"),
    "d_prot": (np.log(2) / 3600, "1/s", "dilucion por division, t_div ~60 min"),

    # --- union represor-inductor (ArsR/SmtB, sistema de represion) ---
    "K_ind": (1.0, "uM", "constante de disociacion represor-metal"),
    "n_hill": (2.0, "-", "cooperatividad del represor dimerico"),

    # --- maduracion del cromoforo: la diferencia clave entre reporteros ---
    # Una proteina fluorescente necesita plegarse Y oxidar su cromoforo.
    # Una cromoproteina tambien madura, pero su lectura no exige excitacion.
    "t_mat_gfp": (1800.0, "s", "maduracion GFP ~30 min"),
    "t_mat_mcherry": (2700.0, "s", "maduracion mCherry ~45 min"),
    "t_mat_chromo": (1200.0, "s", "maduracion cromoproteina ~20 min"),

    # --- umbrales de deteccion ---
    # Cuanta proteina madura hace falta para que la senal sea detectable.
    # A simple vista se necesita ~2 ordenes de magnitud mas que con un lector.
    "thr_visual": (5000.0, "nM", "umbral de color visible a simple vista"),
    "thr_instrument": (50.0, "nM", "umbral con fluorimetro / lector de placa"),

    # --- kill switch MazEF ---
    "k_mazF": (0.10, "nM/s", "sintesis de toxina, J23117 debil + RBS medio"),
    "k_mazE": (0.60, "nM/s", "sintesis de antitoxina, J23100 fuerte + RBS fuerte"),
    "k_on": (0.01, "1/(nM s)", "formacion del complejo MazE-MazF"),
    "k_off": (1e-4, "1/s", "disociacion del complejo"),
    "d_mazE": (np.log(2) / 600, "1/s", "MazE degradada por ClpAP, t1/2 ~10 min"),
    "d_mazF": (np.log(2) / 5400, "1/s", "MazF sin tag, t1/2 ~90 min (v1)"),
    "d_mazF_ssrA": (np.log(2) / 240, "1/s",
                    "MazF con tag ssrA degradada por ClpXP, t1/2 ~4 min (v2)"),
    "tox_lethal": (50.0, "nM", "MazF libre por encima del cual la celula muere"),
}


def p(name: str) -> float:
    """Devuelve el valor numerico de un parametro."""
    return PARAMS[name][0]


def d_toxin(ssrA: bool = False) -> float:
    """Degradacion de MazF, con o sin el tag ssrA de la v2."""
    return p("d_mazF_ssrA") if ssrA else p("d_mazF")


def max_toxin_synthesis(ssrA: bool = False) -> float:
    """Sintesis maxima de MazF compatible con una celula viable, en nM/s.

    En estado estacionario toda MazF sintetizada debe salir del pool libre.
    Secuestrarla en el complejo no la elimina: ClpAP degrada la MazE del
    complejo y devuelve la MazF intacta al pool. El unico sumidero real es
    la degradacion de la propia MazF, asi que la condicion de viabilidad es

        k_mazF < d_mazF * tox_lethal

    Aumentar la sintesis de antitoxina NO relaja esta cota: solo retrasa el
    momento en que se alcanza.
    """
    return d_toxin(ssrA) * p("tox_lethal")


def is_viable(ssrA: bool = False) -> bool:
    """True si el diseno permite una celula viable en equilibrio."""
    return p("k_mazF") < max_toxin_synthesis(ssrA)


# Eficiencia de iniciacion de la traduccion en B. subtilis segun cuantos
# nucleotidos contiguos aparea el SD con el 3' del 16S. La relacion es
# monotona y fuertemente no lineal: por debajo de 5 nt la iniciacion cae.
# Valores aproximados, usados para comparar disenos, no para predecir.
SD_EFFICIENCY = {4: 0.15, 5: 0.40, 6: 0.75, 7: 1.00}


def rbs_efficiency_from_sd(pairing_nt: int) -> float:
    """Eficiencia relativa de un RBS a partir de su apareamiento SD."""
    if pairing_nt in SD_EFFICIENCY:
        return SD_EFFICIENCY[pairing_nt]
    return 1.0 if pairing_nt > 7 else 0.05


@dataclass
class SensorConfig:
    """Configuracion de una simulacion de biosensor."""

    reporter: str = "chromoprotein"   # gfp | mcherry | chromoprotein
    readout: str = "visual"           # visual | instrument
    rbs_efficiency: float = 1.0       # 1.0 = espaciador correcto; <1 = defecto
    inducer: float = 10.0             # uM de contaminante

    def t_maturation(self) -> float:
        return {
            "gfp": p("t_mat_gfp"),
            "mcherry": p("t_mat_mcherry"),
            "chromoprotein": p("t_mat_chromo"),
        }[self.reporter]

    def threshold(self) -> float:
        return p("thr_visual") if self.readout == "visual" else p("thr_instrument")


# --------------------------------------------------------------------------
# Modelo 1-2: biosensor de represion
# --------------------------------------------------------------------------

def _sensor_rhs(t: float, y: np.ndarray, cfg: SensorConfig) -> list[float]:
    """mRNA -> proteina inmadura -> proteina madura (la que se ve)."""
    mrna, imm, mat = y

    # Promotor reprimido: el inductor secuestra al represor y lo libera.
    # Hill creciente porque ArsR/SmtB se suelta del DNA al unir metal.
    x = (cfg.inducer / p("K_ind")) ** p("n_hill")
    activation = p("leak") + (1.0 - p("leak")) * x / (1.0 + x)

    d_mrna = p("k_txn_max") * activation - p("d_mrna") * mrna
    k_mat = 1.0 / cfg.t_maturation()

    d_imm = (p("k_tln") * cfg.rbs_efficiency * mrna
             - k_mat * imm - p("d_prot") * imm)
    d_mat = k_mat * imm - p("d_prot") * mat
    return [d_mrna, d_imm, d_mat]


def simulate_sensor(cfg: SensorConfig | None = None, t_end: float = 14400.0):
    """Integra el biosensor. Devuelve (t, y) con y = [mRNA, inmadura, madura]."""
    cfg = cfg or SensorConfig()
    sol = solve_ivp(
        _sensor_rhs, (0.0, t_end), [0.0, 0.0, 0.0], args=(cfg,),
        method="LSODA", dense_output=True, max_step=300.0,
        rtol=1e-8, atol=1e-10,
    )
    return sol.t, sol.y


def time_to_signal(cfg: SensorConfig | None = None,
                   t_end: float = 14400.0) -> float | None:
    """Minutos hasta que la proteina madura cruza el umbral de lectura.

    Devuelve None si no lo cruza dentro de t_end.
    """
    cfg = cfg or SensorConfig()
    t, y = simulate_sensor(cfg, t_end)
    mat, thr = y[2], cfg.threshold()
    above = np.flatnonzero(mat >= thr)
    if above.size == 0:
        return None
    i = above[0]
    if i == 0:
        return 0.0
    # interpolacion lineal entre los dos puntos que cruzan el umbral
    t0, t1 = t[i - 1], t[i]
    m0, m1 = mat[i - 1], mat[i]
    t_cross = t0 + (thr - m0) * (t1 - t0) / (m1 - m0)
    return t_cross / 60.0


def dose_response(doses: np.ndarray, cfg: SensorConfig | None = None,
                  t_read: float = 1800.0) -> np.ndarray:
    """Proteina madura a t_read (default 30 min) para cada dosis de inductor."""
    cfg = cfg or SensorConfig()
    out = []
    for d in doses:
        c = SensorConfig(cfg.reporter, cfg.readout, cfg.rbs_efficiency, float(d))
        sol = solve_ivp(_sensor_rhs, (0.0, t_read), [0.0, 0.0, 0.0], args=(c,),
                        method="LSODA", max_step=300.0, rtol=1e-8, atol=1e-10)
        out.append(sol.y[2, -1])
    return np.array(out)


# --------------------------------------------------------------------------
# Modelo 3: kill switch MazEF
# --------------------------------------------------------------------------

@dataclass
class KillSwitchConfig:
    """Configuracion del kill switch.

    rbs_mazE / rbs_mazF permiten simular el defecto de espaciado RBS->ATG
    afectando a cada unidad por separado, que es el escenario peligroso.
    """

    rbs_mazE: float = 1.0
    rbs_mazF: float = 1.0
    plasmid_lost_at: float | None = None  # s; None = la celula retiene el plasmido
    ssrA: bool = False                    # True = MazF lleva el tag de la v2


def _ks_rhs(t: float, y: np.ndarray, cfg: KillSwitchConfig) -> list[float]:
    """MazE libre, MazF libre, complejo MazE-MazF."""
    e, f, c = y

    # Al perder el plasmido cesa la sintesis de ambas proteinas.
    on = 1.0 if (cfg.plasmid_lost_at is None or t < cfg.plasmid_lost_at) else 0.0

    bind = p("k_on") * e * f
    unbind = p("k_off") * c

    # ClpAP degrada la MazE del complejo y LIBERA la MazF, que es estable.
    # Modelar el complejo como si se degradara entero haria desaparecer la
    # toxina junto con la antitoxina y anularia la contencion.
    release = p("d_mazE") * c

    d_e = p("k_mazE") * cfg.rbs_mazE * on - p("d_mazE") * e - bind + unbind
    d_f = (p("k_mazF") * cfg.rbs_mazF * on - d_toxin(cfg.ssrA) * f
           - bind + unbind + release)
    d_c = bind - unbind - release
    return [d_e, d_f, d_c]


def simulate_killswitch(cfg: KillSwitchConfig | None = None,
                        t_end: float = 21600.0, dense: bool = True):
    """Integra el kill switch. Devuelve (t, y) con y = [MazE, MazF, complejo].

    `dense=False` omite la interpolacion continua y afloja la tolerancia: es
    lo que usa el muestreo de Sobol, donde importa el agregado sobre miles de
    corridas y no la trayectoria individual.
    """
    cfg = cfg or KillSwitchConfig()
    # Sistema stiff: k_on es rapido frente a las degradaciones, y d_mazE y
    # d_mazF difieren un orden de magnitud. LSODA conmuta solo entre metodos.
    # La discontinuidad al perder el plasmido se declara con t_eval para que
    # el integrador no la atraviese con un paso largo.
    kw = dict(method="LSODA", max_step=300.0, rtol=1e-8, atol=1e-10)
    if dense:
        kw["dense_output"] = True
    else:
        kw.update(rtol=1e-6, atol=1e-8, max_step=600.0,
                  t_eval=np.linspace(0.0, t_end, 361))
    sol = solve_ivp(_ks_rhs, (0.0, t_end), [0.0, 0.0, 0.0], args=(cfg,), **kw)
    return sol.t, sol.y


def time_to_death(cfg: KillSwitchConfig | None = None,
                  t_end: float = 21600.0, dense: bool = True) -> float | None:
    """Minutos hasta que MazF libre supera el umbral letal.

    Devuelve None si la celula sobrevive dentro de t_end -- que es el
    resultado deseado mientras el plasmido se mantiene.
    """
    cfg = cfg or KillSwitchConfig()
    t, y = simulate_killswitch(cfg, t_end, dense=dense)
    f = y[1]
    above = np.flatnonzero(f >= p("tox_lethal"))
    if above.size == 0:
        return None
    i = above[0]
    if i == 0:
        return 0.0
    t0, t1 = t[i - 1], t[i]
    f0, f1 = f[i - 1], f[i]
    return (t0 + (p("tox_lethal") - f0) * (t1 - t0) / (f1 - f0)) / 60.0
