# Pies de figura

El texto explicativo vive aquí, no dentro de la imagen. Dentro de la figura
solo queda lo que no se puede decir en el pie: unidades, categorías, umbrales
y etiquetas de serie.

Estructura de cada pie, según convención Nature Portfolio:
`Fig. N | ` + título nominal en negrita → paneles `a`/`b`/`c` en presente y
estilo telegráfico → parámetros y procedencia → disponibilidad de datos.

---

## Fig. 1 | Cinética de señal del biosensor frente al umbral visual

**a** Reportero maduro frente al tiempo para las tres versiones del
constructo, con el plásmido inducido a 10 µM de contaminante. La línea gris
horizontal discontinua marca el umbral de color visible a simple vista
(5000 nM); la vertical punteada, los 30 min que afirma el proyecto. El punto
señala el cruce del umbral. Solo la v3 lo alcanza, y lo hace a los 50 min.

**b** Dosis-respuesta a los 30 min para tres reporteros, con las mismas líneas
de umbral: la superior es el umbral visual y la inferior el de instrumento
(50 nM). Ninguno de los tres cruza el umbral visual a ninguna concentración
ensayada, de modo que el techo a los 30 min lo impone la maduración del
cromóforo y no la dosis de contaminante. La lectura con instrumento sí es
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
ST ≥ 0.02 ordenados de mayor a menor. ST mucho mayor que S1 indica que el
parámetro actúa sobre todo por interacción. El umbral letal de MazF domina las
dos conclusiones de contención; el tiempo hasta señal depende de forma
repartida de la transcripción, la traducción y el umbral visual.

**b** Margen de MazF libre bajo el umbral letal en reposo frente al retardo
hasta la muerte tras la pérdida del plásmido, sobre las 1920 muestras. En negro
las que cumplen las dos condiciones de contención; en bermellón las que
disparan con el plásmido retenido (margen negativo). La banda gris superior
reúne las muestras que nunca alcanzan nivel letal dentro de 6 h; su posición
horizontal conserva el margen en reposo, y el eje vertical está partido porque
esas muestras no tienen un retardo definido.

La banda llega a márgenes mayores que la nube inferior: al crecer el margen en
reposo, la toxina deja de alcanzar concentración letal. El 60% de las muestras
es segura pero no contiene, frente al 26% que cumple ambas condiciones y al 14%
que dispara sola. Seguridad y contención compiten por el mismo parámetro.

n = 1920 evaluaciones sobre una muestra de Sobol de 13 parámetros, con rangos
de [valor/factor, valor×factor] según su incertidumbre. Paleta Okabe-Ito.

Código fuente: `model/fig3_sensitivity.py`.

---

## Figuras suplementarias

Conservan cada panel por separado, en el formato original, para el material
adicional. Se generan con `model/figS_supplementary.py`.

**Fig. S1 | Curso temporal del reportero maduro.** Las tres versiones del
constructo con el plásmido inducido a 10 µM de contaminante. Es el panel a de
la Fig. 1, a columna simple.

**Fig. S2 | Dosis-respuesta a los 30 min.** Tres reporteros frente a los
umbrales visual (5000 nM) y de instrumento (50 nM). Es el panel b de la Fig. 1.

**Fig. S3 | Cota analítica de viabilidad del kill switch.** Tiempo hasta la
muerte espuria frente a la tasa de síntesis de MazF. La banda verde marca el
rango en que la célula es viable; las verticales, la cota de la v1
(0.0064 nM/s), la de la v2 con tag ssrA (0.144 nM/s) y el valor del diseño
(0.10 nM/s). El tag desplaza la cota 22× y deja el diseño del lado viable.
Sostiene por vía analítica la misma conclusión que la Fig. 2 obtiene por
simulación.

**Fig. S4 | Reparto de las muestras entre los cuatro cuadrantes de
contención.** Fracción de las 1920 muestras que cumple ambas condiciones, solo
una, o ninguna. El modo de fallo dominante es no matar, no disparar sola.
