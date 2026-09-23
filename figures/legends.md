# Pies de figura

El texto explicativo vive aquí, no dentro de la imagen. Dentro de la figura
solo queda lo que no se puede decir en el pie: unidades, categorías, umbrales
y etiquetas de serie.

Estructura de cada pie, según convención Nature Portfolio:
`Fig. N | ` + título nominal en negrita → paneles `a`/`b`/`c` en presente y
estilo telegráfico → parámetros y procedencia → disponibilidad de datos.

---

## Fig. 1 | Cinética de senal del biosensor frente al umbral visual

**a** Reportero maduro frente al tiempo para las tres versiones del
constructo, con el plásmido inducido a 10 µM de contaminante. La linea gris
horizontal discontinua marca el umbral de color visible a simple vista
(5000 nM); la vertical punteada, los 30 min que afirma el proyecto. El punto
señala el cruce del umbral. Solo la v3 lo alcanza, y lo hace a los 50 min.

**b** Dosis-respuesta a los 30 min para tres reporteros, con las mismas lineas
de umbral: la superior es el umbral visual y la inferior el de instrumento
(50 nM). Ninguno de los tres cruza el umbral visual a ninguna concentración
ensayada, de modo que el techo a los 30 min lo impone la maduracion del
cromóforo y no la dosis de contaminante. La lectura con instrumento si es
alcanzable en toda la meseta.

Simulación determinista del sistema de ODEs de `model/kinetics.py`, integrada
con LSODA. Los parámetros provienen de literatura sobre sistemas análogos, no
de mediciones de este proyecto, y sirven para órdenes de magnitud. Paleta
Okabe-Ito.

Código fuente: `model/fig1_sensor.py`.

---

## Fig. 2 | Contención del kill switch MazEF en las tres versiones del diseño

**a**–**f** Antitoxina MazE libre (negro) y toxina MazF libre (bermellón) tras
la inducción, en escala logarítmica, para las tres versiones del constructo
(filas) con el plásmido retenido (columna izquierda) o perdido a los 60 min
(columna derecha). La línea gris horizontal discontinua marca el umbral letal
de MazF libre (50 nM); la vertical punteada, el momento de la pérdida del
plásmido. El rótulo de cada panel indica si la célula sobrevive o el tiempo
hasta la muerte. La contención exige «sobrevive» en la columna izquierda y
«muere» en la derecha. Paleta Okabe-Ito.

**a**, **b** v1 (B0032/B0034, MazF sin tag ssrA): la toxina alcanza nivel
letal con el plásmido retenido y no lo alcanza al perderlo; el
comportamiento es el opuesto al buscado. **c**, **d** v2 (mismos RBS, MazF con
tag ssrA): el tag evita el disparo espurio, pero con un apareamiento
Shine-Dalgarno de 4–5 nt la toxina nunca llega a concentración letal y la
contención se pierde. **e**, **f** v3 (RBS optimizados para *B. subtilis*,
MazF con tag ssrA): única versión que discrimina entre los dos destinos,
sobrevive con el plásmido y muere 74 min después de perderlo.

Simulación determinista del sistema de ODEs de `model/kinetics.py`, integrada
con LSODA. Los parámetros provienen de literatura sobre sistemas análogos en
*E. coli* y *B. subtilis*, no de mediciones de este proyecto, y sirven para
órdenes de magnitud.

Código fuente y datos: `model/fig2_killswitch.py`.

---

## Fig. 3 | Robustez de las conclusiones bajo incertidumbre paramétrica

**a** Índices de Sobol de primer orden (S1, tramo oscuro) y total (ST, barra
completa) para las tres conclusiones del modelo, con los parámetros de
ST ≥ 0.02 ordenados de menor a mayor. ST mucho mayor que S1 indica que el
parámetro actúa sobre todo por interacción. El umbral letal de MazF domina las
dos conclusiones de contención; el tiempo hasta señal depende de forma
repartida de la transcripción, la traducción y el umbral visual.

**b** Margen de MazF libre bajo el umbral letal en reposo frente al retardo
hasta la muerte tras la perdida del plásmido, para las muestras que si
alcanzan nivel letal dentro de 6 h. En negro las que cumplen las dos
condiciones de contención; en bermellón las que disparan con el plásmido
retenido (margen negativo). **c** Reparto de todas las muestras entre los
cuatro cuadrantes, en una barra apilada al 100% que incluye las censuradas que
el panel b no puede representar.

El modo de fallo dominante no es el disparo espurio sino la ausencia de
muerte: el 60% de las muestras es segura en reposo pero nunca alcanza
concentración letal tras el escape. Con promotores constitutivos, seguridad y
contención compiten por el mismo parametro.

n = 1920 evaluaciones del modelo sobre una muestra de Sobol de 13 parámetros,
con rangos de [valor/factor, valor*factor] según la incertidumbre de cada uno.
Paleta Okabe-Ito.

Código fuente: `model/fig3_sensitivity.py`.

---

## Extended Data

**ED1** — Barrido de la tasa de síntesis de MazF frente a la cota analítica de
viabilidad. Era la figura 04. Sostiene la misma afirmación que la Fig. 2 por
vía analítica en vez de por simulación, así que no gana un panel en la figura
principal, pero documenta de dónde sale la cota `k_F = d_F · F_letal` y cuánto
la desplaza el tag ssrA.
