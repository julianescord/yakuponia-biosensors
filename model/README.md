# Modelo cinético

Modelos de ODEs que responden tres preguntas que el diseño de secuencia por sí
solo no puede contestar.

> **Los parámetros no son mediciones de este proyecto.** Son valores de
> literatura para sistemas análogos en *E. coli* / *B. subtilis*, elegidos para
> obtener **órdenes de magnitud**, no predicciones cuantitativas. Cada uno está
> en `PARAMS` con su unidad y su justificación. Las conclusiones que siguen son
> robustas al valor exacto; lo que importa es la **estructura** del sistema.

```bash
.venv/bin/python model/run_analysis.py
```

---

## 1. ¿Es alcanzable el "<30 minutos"?

**No para lectura a simple vista, ni siquiera en la v2.** Ningún reportero
cruza el umbral visual dentro de 30 min:

| Constructo | A simple vista | Con instrumento |
|---|---|---|
| v1 FtaGen (GFP, RBS 0 nt) | **no alcanza** | 14 min |
| v1 MetalGen (mCherry) | 58 min | 6 min |
| **v2 (cromoproteína + espaciador)** | **41 min** | **5 min** |

La v2 restaura la función del sensor de ftalatos, que en v1 no llegaba nunca,
y baja el tiempo de 58 a 41 min en MetalGen. Pero no rompe el piso.

![tiempo hasta señal](figures/01_time_to_signal.png)

El cuello de botella no es la transcripción ni la traducción: es la
**maduración del cromóforo**, que impone un piso de 20–45 min según el
reportero. Ninguna optimización del promotor o del RBS lo evita.

La cromoproteína **sí era la elección correcta** —es la más rápida y no
necesita excitación— pero no basta para sostener el "<30 min". Las opciones
honestas son:

- ajustar la afirmación a **~45 minutos** para lectura visual, o
- mantener los 30 min y declarar que requiere lector.

![dosis-respuesta](figures/02_dose_response.png)

A los 30 min la respuesta ya es sigmoidea y discrimina dosis, pero se queda por
debajo del umbral visual en todo el rango: **hay señal medible antes de que
haya señal visible.**

## 2. ¿Cuánto cuesta el defecto de espaciado RBS→ATG?

Es peor que una simple demora. Por debajo del 20 % de eficiencia traduccional,
el sensor **nunca** alcanza el umbral visual:

| Eficiencia de traducción | Tiempo hasta señal | |
|---|---|---|
| **100 %** | **41 min** | v2, espaciador de 7 nt |
| 50 % | 72 min | |
| 20 % | no alcanza | |
| 10 % | no alcanza | v1, espaciado 0 nt |

La síntesis compite con la dilución por división celular. Por debajo de cierto
umbral la proteína madura se estabiliza en una meseta bajo el nivel detectable,
y esperar más no ayuda. El defecto no retrasa el sensor: **lo inutiliza**.

## 3. ¿El kill switch contiene la célula?

**En v1 no, por una razón estructural. En v2 sí.** El hallazgo más importante
del modelo, y el que motivó la corrección.

### La cota de viabilidad

Secuestrar la toxina en el complejo MazE–MazF **no la elimina**: ClpAP degrada
la MazE del complejo y libera la MazF intacta, que es estable (t½ ≈ 90 min
frente a 10 min de MazE). El único sumidero real de toxina es su propia
degradación. En estado estacionario:

$$k_F < d_F \cdot F_{letal}$$

| | cota admisible | k_mazF del diseño | |
|---|---|---|---|
| **v1** (MazF sin tag) | 0.0064 nM/s | 0.1 nM/s | **16× por encima** |
| **v2** (MazF-ssrA) | 0.1444 nM/s | 0.1 nM/s | dentro de la cota |

La corrección no toca promotores ni RBS: el tag **ssrA** (AANDENYALAA) en el
C-terminal de MazF la hace sustrato de ClpXP, bajando su t½ de ~90 a ~4 min.
Eso sube `d_mazF` y con ello la cota **22×**, sin alterar la lógica del
circuito.

**Aumentar la síntesis de antitoxina no relaja esta cota** — solo retrasa el
momento en que se alcanza. Es una propiedad del sistema, no un parámetro a
afinar.

![barrido de toxina](figures/04_toxin_sweep.png)

### Consecuencia: el switch dispara solo

| Escenario | v1 | v2 | Debería |
|---|---|---|---|
| Plásmido retenido | muere **95 min** ✗ | **sobrevive** ✓ | sobrevivir |
| Plásmido perdido a 60 min | muere 67 min ✓ | muere **68 min** ✓ | morir |

![kill switch](figures/03_killswitch.png)

En v1 el margen entre ambos escenarios era de 28 min: la célula moría tuviera o
no el plásmido. **La v2 discrimina correctamente.**

### Los márgenes de la v2 son estrechos

Conviene declararlo en vez de presentar la corrección como resuelta:

- En reposo, MazF libre se estabiliza en **35 nM contra un umbral de 50 nM**:
  quedan solo 15 nM de margen ante variación de parámetros.
- Tras perder el plásmido, la toxina supera el umbral entre los 68 y los 86
  min: una **ventana letal de 18 min**. Si la muerte no ocurre dentro de ella,
  la célula se recupera.

Ambos márgenes dependen de parámetros no medidos. Una corrección más robusta
sería **hacer la toxina condicional** en vez de constitutiva: que MazF solo se
exprese ante la señal de escape, eliminando la carrera permanente entre
síntesis y degradación. Eso es un rediseño, no un ajuste.

---

## Estructura

| Archivo | Qué hace |
|---|---|
| `kinetics.py` | ODEs, parámetros y la cota analítica de viabilidad |
| `run_analysis.py` | corre el análisis, imprime la tabla y escribe las figuras |
| `figures/` | salida, regenerable — no editar a mano |

Integración con **LSODA**, que conmuta solo entre métodos stiff y no-stiff: el
sistema lo es, porque `k_on` es rápido frente a las degradaciones y `d_mazE` y
`d_mazF` difieren un orden de magnitud.

## Qué falta

- Ajustar parámetros contra datos experimentales propios
- Análisis de sensibilidad global (Sobol) para saber qué parámetros dominan
- Modelo estocástico: a bajo número de copias el ruido puede disparar el switch
- Carga metabólica del plásmido sobre el crecimiento
