# Yakuponía — Biosensores bacterianos para cultivos hidropónicos

Diseño *in silico* de tres constructos genéticos para la detección colorimétrica
de contaminantes en agua de cultivos hidropónicos, más un sistema de contención
biológica. Desarrollado por el equipo **The Operons (U2)** en el SynBio 6.0
Bootcamp.

> **Estado: diseño computacional. Ninguna construcción ha sido ensamblada ni
> validada experimentalmente.** Las limitaciones conocidas están listadas
> abajo y son detectadas automáticamente por `scripts/validate_assembly.py`.

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

### 1. Los reporteros son fluorescentes, no cromoproteínas — *bloqueante*

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

### 2. Espaciado RBS→ATG de 0 nt — *bloqueante*

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

**Corrección:** insertar espaciador de 6–8 nt en las tres uniones.

### 3. Estándares de ensamblaje mezclados

Hay un prefijo **BioBrick RFC[10]** que no aparece en ningún constructo
(verificado), sitios **BsaI** (Golden Gate/MoClo) y sitios EcoRI/XbaI/SpeI/PstI
internos. BioBrick y Golden Gate son estándares incompatibles; hay que elegir
uno y domesticar el resto de la secuencia en consecuencia.

### 4. Partes de *E. coli* en chasis de *B. subtilis*

pHT01 y pWB980 son vectores de *B. subtilis*, pero los promotores Anderson, los
RBS B0032/B0034 y MazEF están caracterizados en *E. coli*. Los Anderson son
σ70/σA y sí funcionan, pero **con fuerzas relativas distintas a las publicadas**
— y el kill switch depende justamente de esas fuerzas relativas. Los RBS son
subóptimos: *B. subtilis* exige Shine-Dalgarno más largo y complementario.

### 5. Tamaño de los plásmidos

9,5 kb y 10,9 kb son grandes para *B. subtilis*: baja eficiencia de
transformación y carga metabólica que reduce la señal del reportero.

### 6. Sin validación experimental

No hay curva dosis-respuesta, ni límite de detección, ni tiempo real hasta
señal. La afirmación de **"resultados en menos de 30 minutos"** que acompañó la
presentación del proyecto **no está respaldada** por ningún dato ni modelo.

---

## Estructura del repositorio

```
sequences/
  raw/         export original de Benchling, sin modificar
  parts/       16 partes modulares individuales
  backbones/   vectores pHT01 y pWB980
  annotated/   ← constructos reanotados, generados por script
scripts/
  reannotate.py         reconstruye la anotación desde las partes
  validate_assembly.py  valida el ensamblaje y reporta defectos
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
.venv/bin/python scripts/validate_assembly.py    # valida
```

El validador sale con código 1 mientras queden defectos abiertos — los 8 que
reporta hoy corresponden a las limitaciones 1 y 2.

### Por qué reanotar

Benchling exporta casi todas las features como `misc_feature` sin metadatos. El
script recupera la anotación buscando cada parte modular dentro de cada
constructo y escribiendo features tipadas (`promoter`, `RBS`, `CDS`,
`terminator`). Esto hace que la anotación sea **reproducible y verificable**, no
un estado manual que se desincroniza. Fue así como apareció el regulador
ArsR/SmtB no anotado.

---

## Próximos pasos

- [ ] Sustituir GFP/mCherry por eforRed/amilCP (limitación 1)
- [ ] Insertar espaciadores RBS→ATG de 6–8 nt (limitación 2)
- [ ] Elegir un único estándar de ensamblaje y domesticar sitios (limitación 3)
- [ ] Recalibrar el ratio toxina:antitoxina para σA de *B. subtilis* (limitación 4)
- [ ] Modelo cinético (ODEs) del circuito: dosis-respuesta y tiempo hasta señal
- [ ] Protocolo de validación húmeda: cepas, concentraciones, controles

---

## Equipo

The Operons (U2) — SynBio 6.0 Bootcamp:
Valeria Quesada (U. Nacional de Costa Rica) · Pamela Elías (U. de El Salvador) ·
Naomi Macanchí (Yachay Tech) · Jovany Reyes (U. Autónoma del Estado de México) ·
Julián Escobar (Yachay Tech)
