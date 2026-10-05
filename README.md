# TEFSM_inv

Genera figuras a partir de lecturas de mV del equipo TEFSM (método de selección de frecuencias del campo eléctrico telúrico).

Basado en:
- Yang, T. et al., *Simulation of the Telluric Electrical Field Frequency Selection Method and Its Application in Mineral Water Exploration*, Water 2025, 17, 3314. https://doi.org/10.3390/w17223314 (figuras 2, 4, 6 y 13).
- Gomo, M. y Ngobe, T., *Telluric Electric Frequency Selection Method (TEFSM) in Geophysical Groundwater Exploration: Emerging Issues*, en *Aquifers – Advances in Hydrogeology*, IntechOpen. https://doi.org/10.5772/intechopen.1013979 (figura 4).

## Contenido

```
tefsm_figuras.py   programa
requirements.txt   dependencias (pip install -r requirements.txt)
datos/             lecturas del equipo (CSV)
figuras/           figuras finales, una carpeta por archivo de entrada (p. ej. L89_150m/ para 150M_L89.csv),
                   cada una con: fig13_*, fig4_horizontal_*, SP_Inversion_* y Modelo_Resistividad_*
```

## Datos de entrada

CSV del equipo, una fila por punto de medida: `L` (registro de memoria), `N` (punto, 80…97) y `freq01…freqNN` (ΔV en mV; 36 o 40 canales).
**El nombre del archivo debe empezar con el rango de profundidad (`100M_`, `150M_`, `300M_` o `500M_`)**; de ahí el programa lo toma y, si falta o no es uno de esos cuatro valores, se detiene con un aviso (o use `--prof`). Ejemplos: `150M_L89.csv` (150 m, registro 89) y `300M_L19.csv` (300 m, registro 19).
El espaciado entre electrodos es de 1 m (`--dx 1`, valor por defecto): 22 electrodos dan 18 puntos, de 0 a 17 m. `freq01` es el canal más somero y `freqNN` el más profundo; la profundidad de cada canal es lineal con su número hasta el rango configurado (`--prof`, o se lee del nombre).

## Figuras

| Comando | Resultado |
|---|---|
| `pseudo` | Figura 13: pseudo-sección normalizada log10(ΔV/ΔVmin), `fig13_*.png` |
| `horizontal` | Figura 4 (Gomo y Ngobe): EPD (log) vs distancia, una curva por profundidad, `fig4_horizontal_*.png` |
| `inversion` | `SP_Inversion_*.png` y `Modelo_Resistividad_*.png` (estación elegida) |
| `modulo` | Figuras 2, 4 y 6 de Yang et al. (|Ey| y fase; CSV con `y_m, f_hz, ey_mv[, fase_deg]`; vienen de una simulación 3D, no del equipo) |

```
python tefsm_figuras.py todas      datos/150M_L89.csv     # las cuatro figuras en figuras/L89_150m/
python tefsm_figuras.py todas      datos/300M_L19.csv     # las cuatro figuras en figuras/L19_300m/

python tefsm_figuras.py pseudo     datos/150M_L89.csv     # o una sola figura
python tefsm_figuras.py horizontal datos/150M_L89.csv
python tefsm_figuras.py inversion  datos/300M_L19.csv --punto 88
```

Cada comando guarda sus PNG en `figuras/<L##_###m>/`, según el nombre del archivo (otra carpeta: `--salida`) e imprime la ruta completa de cada archivo; los nombres por defecto son `fig13_*`, `fig4_horizontal_*`, `SP_Inversion_*` y `Modelo_Resistividad_*`. Por defecto no se abre ninguna ventana: añada `--mostrar` para ver la figura en pantalla además de guardarla.

Opciones útiles (`-h` en cada comando muestra todas):
- `--umbral` (0.1 mV): lecturas menores se consideran ruido o canal muerto.
- `pseudo`: `--zk x` marca un sondeo; `--escala-h ec12` usa hs = c·503·√(ρ/f) del artículo en lugar de la profundidad lineal.
- `pseudo` e `inversion`: sin `--punto`/`--x` usan la misma estación (la de mayor anomalía baja).
- `horizontal`: `--zona x1,x2` recuadra una zona (p. ej. de ruido); `--cada N` (def.: el necesario para ~18 curvas) / `--profundidades`.
- `inversion`: `--punto N` o `--x m`; sin ellos usa la estación de mayor anomalía baja (sin contar los bordes).

## Sobre el modelo de resistividad

Gomo y Ngobe indican que los valores del equipo no se pueden convertir en resistividad. Por eso `Modelo_Resistividad_*.png` es una **aproximación empírica y cualitativa**, no una inversión física:
- *Calculado*: ajuste regularizado (Tikhonov, 2.ª derivada) de la EPD observada contra la profundidad (`--lam`); *Suavizado*: filtro gaussiano (`--factor-suavizado`). RMS = diferencia entre el ajuste y las lecturas.
- *Etiquetas y recuadros*: el modelo muestra ~8 valores de resistividad repartidos por la curva (`--etiquetas N`: extremos, máximos y mínimos locales y profundidades regulares, sin encimarse). Un recuadro fucsia marca un tramo donde la resistividad baja respecto a la capa anterior, por lo menos 15 % (`--prominencia 0.06`, en log10); si la curva es monótona no aparece ninguno. En la pseudo-sección (`pseudo`) los recuadros fucsia son de dos tipos: trazo **continuo** grueso = los mismos tramos de menor resistividad de la estación analizada (la línea vertical negra); trazo **discontinuo** = zonas de EPD baja respecto a la mediana lateral a cada profundidad (`--umbral-anomalia`, `--min-celdas`). Solo se dibujan con datos del equipo y profundidad lineal.
- *Resistividad*: mapeo log-lineal de la EPD ajustada al intervalo `--rho-min 5` … `--rho-max 5000` Ω·m (más EPD, más resistividad). Calíbrelo con un sondeo o con resistividad eléctrica local.

### Valores de referencia en Guatemala (aportados por el equipo de campo)

Resistividades aparentes típicas, útiles para elegir `--rho-min` y `--rho-max` (no están verificadas aquí contra sus fuentes):

| Zona / material | Resistividad (Ω·m) |
|---|---|
| Pómez y cenizas secas (superficiales) | 500 – más de 2000 |
| Flujos piroclásticos y toba, saturados o semisaturados | 60 – 250 |
| Basaltos, andesitas y lavas masivas y secas | > 1000 |
| Lavas fracturadas y saturadas (recarga hídrica) | 150 – 400 |
| Rellenos aluviales y suelos areno-arcillosos | 30 – 150 |
| Sedimentos lacustres y arcillas saturadas | 5 – 30 |
| Calizas kársticas densas y secas | 800 – más de 5000 |
| Calizas con cavernas colmatadas o acuífero kárstico | 80 – 200 |
| Esquistos, gneises y serpentinitas sanos | > 800 |
| Esquistos meteorizados (horizontes arcillosos) | 40 – 120 |

Los valores por defecto (5 y 5000 Ω·m) cubren desde las arcillas lacustres hasta las calizas y rocas metamórficas sanas. Si en un sitio el rango esperado es más estrecho, ajústelos (por ejemplo `--rho-min 30 --rho-max 2500` en tierras altas volcánicas).
