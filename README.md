# Yakuponía — Biosensores bacterianos para cultivos hidropónicos

Diseño *in silico* de tres constructos genéticos para la detección colorimétrica
de contaminantes en agua de cultivos hidropónicos, más un sistema de contención
biológica. Desarrollado por el equipo **The Operons (U2)** en el SynBio 6.0
Bootcamp.

> **Estado: diseño computacional. Ninguna construcción ha sido ensamblada ni
> validada experimentalmente.** Las limitaciones conocidas están listadas
> abajo y son detectadas automáticamente por `scripts/validate_assembly.py`.
>
> El diseño va por su **versión 3**. Cada versión corrigió defectos que el
> modelo cinético reveló, y la siguiente encontró que la anterior no bastaba.
> Las tres se conservan: `sequences/annotated/` (v1), `sequences/v2/`,
> `sequences/v3/`.

---

## El problema

Los cultivos hidropónicos recirculan agua. Si esa agua transporta **ftalatos**
(lixiviados de las tuberías y films plásticos del propio sistema) o **metales
pesados** (de fertilizantes y fuentes de agua), el contaminante se concentra en
el circuito y llega al producto final.

Detectarlos hoy requiere HPLC o espectrometría de absorción atómica: equipo de
laboratorio, personal técnico y días de espera. Para un productor pequeño, eso
equivale a no medir.

**Hipótesis de diseño:** un reportero transcripcional acoplado a un regulador
nativo de respuesta al contaminante puede dar una lectura visual en campo, sin
instrumentación.

---

## Arquitectura

Tres constructos, dos chasis de *Bacillus subtilis*:

| Constructo | Plásmido | Tamaño | Backbone | Función | Mapa |
|---|---|---|---|---|---|
| **FtaGen** | `pYAK_FtaGen` | 9 541 bp | pHT01 | detección de ftalatos | [pdf](docs/maps/ftagen_phthalate_biosensor_map.pdf) |
| **MetalGen** | `pYAK_MetalGen` | 10 920 bp | pHT01 | detección de metales pesados | [pdf](docs/maps/metalgen_heavy_metal_biosensor_map.pdf) |
| **Kill switch** | `pYAK_KillSw` | 4 720 bp | pWB980 | contención biológica | [pdf](docs/maps/killswitch_mazEF_map.pdf) |

Los tres comparten la misma lógica: `promotor inducible → RBS → reportera → terminador`.

El comportamiento dinámico de los tres está modelado en [model/](model/), que
es donde aparecieron los defectos que la secuencia por sí sola no revela.

### Versiones

| | v1 | v2 | v3 |
|---|---|---|---|
| Reportero FtaGen | GFP | eforRed | eforRed |
| Reportero MetalGen | mCherry | amilCP | amilCP |
| Espaciado RBS→ATG | 0 nt | 7 nt | 7 nt |
| RBS del reportero MetalGen | ausente | propio | propio |
| MazF | sin tag | ssrA | ssrA |
| RBS | B0032/B0034 | B0032/B0034 | **optimizado *B. subtilis*** |
| Apareamiento SD | 4–5 nt | 4–5 nt | **6–7 nt** |
| **Sensor de ftalatos** | no da señal | **no da señal** | **50 min** |
| **Kill switch con plásmido** | muere 304 min ✗ | sobrevive ✓ | sobrevive ✓ |
| **Kill switch sin plásmido** | sobrevive ✗ | **sobrevive ✗** | **muere 74 min ✓** |

Simulado en chasis *B. subtilis*. La v2 corrigió el espaciado pero conservó
RBS de *E. coli*, cuyo Shine-Dalgarno no aparea lo suficiente con el 16S de
*B. subtilis*: el sensor seguía sin dar señal y el kill switch no contenía.
**Solo la v3 funciona en el chasis real.**

Los insertos v3 están en `sequences/v3/` y pasan el validador **sin fallas**.

### FtaGen — sensor de ftalatos

```
P_pht + regulador ──▶ RBS B0032 ──▶ reportera ──▶ BBa_B1002
   1..1095            1096..1108   1109..1777    1778..1818
```

La región reguladora de 1 095 bp contiene el regulador transcripcional y su
promotor cognado. Se eligió **RBS medio (BBa_B0032)** en lugar del fuerte para
limitar la expresión basal: en un sensor, un fondo alto destruye la relación
señal/ruido más de lo que la sensibilidad gana con más reportera.

### MetalGen — sensor de metales pesados

```
P_mer + regulador ──▶ RBS B0034 ──▶ [arsR ◀──] ──▶ mCherry ──▶ BBa_B1002
   1..2109            2110..2121     2122..2445    2446..3156  3157..3197
```

**Hallazgo del análisis de marcos de lectura:** la región de 324 bp entre el RBS
y la reportera —que Benchling exportaba sin anotar— es un **CDS completo de 107
aa en la hebra reversa**, con el motivo hélice-giro-hélice `TVSHQLRFL`
característico de la familia **ArsR/SmtB** de represores metalorregulados. Es
decir, el represor que hace funcionar el sensor estaba presente pero invisible
en el diseño original. Está anotado en `sequences/annotated/`.

Aquí el RBS es **fuerte (BBa_B0034)**: a diferencia de FtaGen, el sistema
ArsR/SmtB es de represión, y el represor mantiene el fondo bajo por sí mismo.

### Kill switch — contención por MazEF

Sistema toxina–antitoxina tipo II de *E. coli*, con las dos unidades
transcripcionales en tándem:

```
J23117 ──▶ RBS B0032 ──▶ mazF (toxina)     ──▶ BBa_B0015
 débil       medio        112 aa               terminador doble

J23100 ──▶ RBS B0034 ──▶ mazE (antitoxina) ──▶ BBa_B0015
 fuerte      fuerte       83 aa                terminador doble
```

**Esta es la decisión de ingeniería central del proyecto.** MazF es una
endoribonucleasa que corta ARN en secuencias ACA y detiene la traducción; MazE
la neutraliza formando un heterocomplejo. La supervivencia celular depende de
que la antitoxina esté en exceso sobre la toxina.

Ese exceso se implementa con un **doble diferencial**, tanto transcripcional
como traduccional:

| | Toxina (MazF) | Antitoxina (MazE) |
|---|---|---|
| Promotor | J23117 (débil) | J23100 (fuerte) |
| RBS | B0032 (medio) | B0034 (fuerte) |

MazE es además intrínsecamente inestable (degradada por ClpAP), así que su
exceso debe sostenerse activamente. Si la célula escapa al ambiente y deja de
replicar el plásmido, el pool de MazE cae primero, MazF queda libre y la célula
muere. La contención no depende de un circuito adicional, sino de la asimetría
de estabilidad entre las dos proteínas.

Ambos usan **BBa_B0015** (terminador doble rrnB T1 + T7TE) para evitar
transcripción de lectura corrida entre unidades: con dos casetes en tándem,
una terminación incompleta alteraría el ratio del que depende todo el sistema.

---

## Inventario de partes

Todas verificadas por búsqueda de subcadena dentro de los constructos
(`scripts/validate_assembly.py`).

| Parte | bp | Identidad | Registro | Usada en |
|---|---|---|---|---|
| J23100 | 35 | promotor constitutivo fuerte | BBa_J23100 | KillSw |
| J23117 | 35 | promotor constitutivo débil | BBa_J23117 | KillSw |
| RBS medio | 12 | `TCACACAGGAAAG` | BBa_B0032 | FtaGen, KillSw |
| RBS fuerte | 13 | `AAAGAGGAGAAA` | BBa_B0034 | MetalGen, KillSw |
| BBa_B0015 | 129 | terminador doble | BBa_B0015 | KillSw |
| BBa_B1002 | 41 | terminador tipo rrnB | BBa_B1002 | FtaGen, MetalGen |
| mazF | 339 | toxina, endoribonucleasa ACA | — | KillSw |
| mazE | 252 | antitoxina | — | KillSw |
| P_pht + reg | 1 095 | regulador de ftalatos | — | FtaGen |
| P_mer + reg | 2 109 | regulador de metales | — | MetalGen |
| arsR-like | 324 | represor ArsR/SmtB, hebra reversa | — | MetalGen |
| mCherry | 711 | proteína fluorescente roja | — | MetalGen |
| GFP | 669 | proteína fluorescente verde | — | FtaGen |
| Prefijo BioBrick | 22 | `EcoRI-NotI-XbaI` RFC[10] | — | huérfana |

---

## Limitaciones conocidas

Declaradas explícitamente porque un diseño que no lista sus defectos no es
revisable. Las tres primeras las detecta el validador automáticamente.

### 1. Reporteros fluorescentes — *corregido en v2*

La propuesta de producto es lectura **a simple vista, sin instrumentación**.
GFP y mCherry **no cumplen eso**: requieren excitación UV/azul y un lector de
fluorescencia. A 30 min y en cultivo diluido, no hay color visible.

**Corrección pendiente:** sustituir por cromoproteínas del registro iGEM, que
producen pigmento visible sin excitación:

| Constructo | Actual | Reemplazo | Color |
|---|---|---|---|
| FtaGen | GFP | `eforRed` (BBa_K592012) | rojo |
| MetalGen | mCherry | `amilCP` (BBa_K592009) | violeta |

El diseño es modular, así que el cambio es un reemplazo de una sola parte: solo
cambia el CDS reportero, el resto del casete queda igual.

**Aplicado en la v2:** `chromoprotein_eforRed` y `chromoprotein_amilCP`.

**Pero el modelo muestra que no basta.** La cromoproteína es la más rápida
(41 min), pero **ninguna cruza el umbral visual en 30 min**: la maduración del
cromóforo impone un piso que no se evita optimizando promotor ni RBS. Queda
ajustar la afirmación a ~45 min, o declarar que los 30 min requieren lector.
Ver [model/](model/).

### 2. Iniciación de la traducción — *corregido en v3*

En FtaGen y en ambas unidades del kill switch, el codón de inicio está pegado
al RBS sin espaciador:

```
TCACACAGGAAAG|ATGAGTGTGATC     ← 0 nt entre RBS y ATG
```

El ribosoma necesita **5–12 nt** entre la secuencia Shine-Dalgarno y el ATG
para posicionarse. Con 0 nt la traducción será muy ineficiente o nula.

Es especialmente grave en el kill switch, porque afecta a toxina y antitoxina
de forma no necesariamente proporcional: el ratio del que depende la contención
podría invertirse.

**Aplicado en la v2:** espaciador de 7 nt en las tres uniones.

El modelo cuantifica el costo de no hacerlo: por debajo del 20 % de eficiencia
traduccional el sensor **nunca** alcanza el umbral, porque la síntesis no supera
a la dilución por división. No era una demora, era pérdida de función — el
FtaGen v1 no llegaba a dar señal visible en ningún tiempo.

**Defecto adicional encontrado al corregir:** el reportero de MetalGen no tenía
RBS propio. El único RBS del casete estaba a 324 nt, separado por el arsR
reverso. La v2 añade un RBS dedicado, haciendo el casete bicistrónico explícito.

**Y el espaciado era solo la mitad del problema.** Los RBS del registro iGEM
están caracterizados en *E. coli*; el 16S de *B. subtilis* exige un
Shine-Dalgarno más largo:

| RBS | Apareamiento SD | Eficiencia |
|---|---|---|
| BBa_B0032 (v1, v2) | 4 nt | 15 % |
| BBa_B0034 (v1, v2) | 5 nt | 40 % |
| Bsub medium (v3) | 6 nt | 75 % |
| Bsub strong (v3) | 7 nt | 100 % |

**Aplicado en la v3:** RBS con SD `AAAGGAGG`/`AAAGGAGGTG`. El validador mide
este apareamiento y falla por debajo de 5 nt.

### 3. El kill switch no contiene — *corregido en v3*

Hallazgo del modelo cinético, no visible en la secuencia.

Secuestrar MazF en el complejo **no la elimina**: ClpAP degrada la MazE del
complejo y libera la MazF intacta, que es estable (t½ ≈ 90 min contra 10 min de
MazE). El único sumidero real de toxina es su propia degradación, lo que impone
una cota dura sobre su síntesis:

```
k_mazF  <  d_mazF · F_letal
 0.1    <    0.0064            ← el diseño está 16x por encima
```

Consecuencia: la célula muere a los 95 min **aunque conserve el plásmido**,
contra 67 min si lo pierde. El margen de 28 min no permite discriminar — un
kill switch que mata al huésped que debía preservar no es contención.

**Aumentar la antitoxina no lo arregla**, solo retrasa el disparo.

**Aplicado en la v2:** tag ssrA (AANDENYALAA) en el C-terminal de MazF, que la
hace sustrato de ClpXP y baja su t½ de ~90 a ~4 min, subiendo la cota **22×**.

**Pero la v2 seguía sin contener.** Con los RBS de *E. coli*, la toxina nunca
alcanzaba concentración letal en *B. subtilis*: la célula sobrevivía al escape.
Solo con los RBS de la **v3** el switch discrimina — sobrevive con el plásmido,
muere a los 74 min al perderlo.

Los márgenes siguen estrechos (26 nM contra un umbral de 50; ventana letal de
10 min) y dependen de parámetros no medidos. Detalle en [model/](model/).

### 4. Estándares de ensamblaje — *resuelto al reensamblar*

En la v1 había un prefijo **BioBrick RFC[10]** huérfano, sitios **BsaI**
(Golden Gate) y varios EcoRI/XbaI/SpeI/PstI internos. Al reensamblar desde las
partes en la v2, los sitios conflictivos desaparecieron:

| | v1 | v2/v3 |
|---|---|---|
| BsaI | 1–2 por constructo | **0** |
| Prefijo BioBrick huérfano | presente | **ausente** |
| EcoRI/SpeI/PstI | varios | 1 EcoRI + 1 SpeI (solo MetalGen) |

Quedan un `EcoRI` y un `SpeI` internos en la región reguladora de MetalGen.
Solo hay que domesticarlos si se elige BioBrick RFC[10]; con Golden Gate son
irrelevantes. Hay además sitios `BsmBI`, relevantes solo para MoClo.

### 5. Partes de *E. coli* en chasis de *B. subtilis* — *parcialmente corregido en v3*

pHT01 y pWB980 son vectores de *B. subtilis*, pero los promotores Anderson, los
RBS B0032/B0034 y MazEF están caracterizados en *E. coli*.

**Los RBS ya están corregidos en la v3.** Era el punto más grave: B0032 solo
aparea 4 nt con el 16S de *B. subtilis*, insuficiente para iniciar traducción.
El validador ahora mide este apareamiento y falla por debajo de 5 nt.

**Los promotores siguen pendientes.** J23117 tiene un −35 consenso (`TTGACA`)
pero un −10 lejano al canónico (`GGATTG` frente a `TATAAT`). Las fuerzas
relativas J23100/J23117 en σA serán distintas a las publicadas, y el kill
switch depende de esa razón. Esto **no se arregla cambiando de parte: se
arregla midiendo** en el chasis. Es trabajo de banco.

### 6. Terminador BBa_B1002 sin cola poli-T

Un terminador rho-independiente necesita horquilla **más cola de U**. BBa_B0015
la tiene (`TTTT`); **BBa_B1002 no tiene ninguna**. Es el terminador de FtaGen y
MetalGen, así que la terminación será ineficiente y puede haber transcripción
de lectura corrida. Sustituible por BBa_B0015, que ya se usa en el kill switch.

### 7. Uso de codones — *evaluado, no bloqueante*

Codones raros para *B. subtilis*: 7,4 % en arsR, 4–5 % en MazEF, **0 % en las
cromoproteínas**. Por debajo del umbral que suele afectar la expresión. No se
recomienda optimizar sin datos que lo justifiquen.

### 8. Tamaño de los plásmidos

9,5 kb y 10,9 kb son grandes para *B. subtilis*: baja eficiencia de
transformación y carga metabólica que reduce la señal del reportero.

### 9. Sin validación experimental

No hay datos propios: ni curva dosis-respuesta medida, ni límite de detección,
ni tiempo real hasta señal. El modelo cinético de [model/](model/) da órdenes
de magnitud a partir de parámetros de literatura, pero **no sustituye la
validación húmeda** — sus parámetros no han sido ajustados contra mediciones de
estos constructos.

---

## Estructura del repositorio

```
sequences/
  raw/         export original de Benchling, sin modificar
  parts/       16 partes modulares individuales
  backbones/   vectores pHT01 y pWB980
  annotated/   ← constructos reanotados, generados por script
scripts/
  reannotate.py         reconstruye la anotación de la v1 desde las partes
  build_corrected.py    ensambla los insertos v2 y v3 desde las partes
  validate_assembly.py  valida ambas versiones y reporta defectos
model/
  kinetics.py           ODEs de los sensores y del kill switch
  run_analysis.py       corre el análisis y escribe las figuras
  figures/              salida, regenerable
docs/
  maps/        mapas circulares de cada constructo (Benchling)
  synbio6-bootcamp-presentation.pdf
requirements.txt
```

Nombres en `snake_case`, sin espacios, para que los scripts y `grep` funcionen
sin comillas. `raw/` conserva el export de Benchling intacto como referencia.

`sequences/annotated/` es **producto de script, no editado a mano**. Para
regenerarlo:

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
python3 scripts/reannotate.py                    # reconstruye la anotación
.venv/bin/python scripts/build_corrected.py      # ensambla v2 y v3
.venv/bin/python scripts/validate_assembly.py    # valida ambas versiones
.venv/bin/python model/run_analysis.py           # análisis cinético
```

El validador sale con código 1 mientras queden defectos abiertos. Hoy reporta
**12: de la v1 y de la v2**, ambas conservadas a propósito como referencia de
lo corregido. Los tres insertos **v3 pasan limpios**.

### Por qué reanotar

Benchling exporta casi todas las features como `misc_feature` sin metadatos. El
script recupera la anotación buscando cada parte modular dentro de cada
constructo y escribiendo features tipadas (`promoter`, `RBS`, `CDS`,
`terminator`). Esto hace que la anotación sea **reproducible y verificable**, no
un estado manual que se desincroniza. Fue así como apareció el regulador
ArsR/SmtB no anotado.

---

## Próximos pasos

- [x] Sustituir GFP/mCherry por eforRed/amilCP (limitación 1)
- [x] Insertar espaciadores RBS→ATG de 7 nt (limitación 2)
- [x] Añadir RBS propio al reportero de MetalGen
- [x] Tag ssrA en MazF para restaurar la contención (limitación 3)
- [x] RBS con Shine-Dalgarno apto para *B. subtilis* (limitación 5, parcial)
- [x] Auditar estándares de ensamblaje y uso de codones (limitaciones 4 y 7)
- [ ] Elegir un único estándar de ensamblaje y domesticar sitios (limitación 4)
- [x] Modelo cinético (ODEs): dosis-respuesta y tiempo hasta señal
- [ ] Ajustar la afirmación de tiempo de respuesta a ~45 min, o declarar lector
- [ ] Evaluar toxina condicional: los márgenes de la v2 son estrechos
- [ ] Ensamblar los insertos v3 en sus backbones (pHT01, pWB980)
- [ ] Sustituir BBa_B1002 por un terminador con cola poli-T (limitación 6)
- [ ] Medir las fuerzas reales de J23100/J23117 en σA de *B. subtilis*
- [ ] Protocolo de validación húmeda: cepas, concentraciones, controles
- [ ] Análisis de sensibilidad global sobre los parámetros del modelo

---

## Equipo

The Operons (U2) — SynBio 6.0 Bootcamp:
Valeria Quesada (U. Nacional de Costa Rica) · Pamela Elías (U. de El Salvador) ·
Naomi Macanchí (Yachay Tech) · Jovany Reyes (U. Autónoma del Estado de México) ·
Julián Escobar (Yachay Tech)
