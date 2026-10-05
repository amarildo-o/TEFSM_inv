# TEFSM_inv

Programa en Python que reproduce las figuras 2, 4, 6 y 13 de:

> Yang, T. et al. *Simulation of the Telluric Electrical Field Frequency Selection Method and Its Application in Mineral Water Exploration*. Water 2025, 17, 3314. https://doi.org/10.3390/w17223314

Las lecturas (mV) se ingresan desde un archivo CSV o Excel.

## Instalación

```
pip install -r requirements.txt
```

## Figuras 2, 4 y 6 (módulo de Ey y fase)

CSV con columnas `y_m, f_hz, ey_mv` (y opcionalmente `fase_deg`, que añade el panel (c) de la fig. 6):

```
python tefsm_figuras.py modulo datos_fig2.csv --nombre fig2
```

Estas figuras provienen de una simulación numérica 3D; el programa solo grafica los valores ingresados.

## Figura 13 (pseudo-sección normalizada de ΔV)

Formato del equipo (`L, N, freq01…freqNN`, en mV; ver `ejemplos/150M_L89.csv`) o formato largo `y_m, f_hz, dv_mv [, rho_ohm_m]`:

```
python tefsm_figuras.py pseudo ejemplos/150M_L89.csv
python tefsm_figuras.py pseudo ejemplos/300M_L19.csv      # sondeo a 300 m, 40 canales
```

- K = log10(ΔV/ΔVmin) (ec. 13) y hs = c·503·√(ρ/f) (ec. 12).
- En el nombre `150M_L89.csv`, `L89` es el registro de memoria del equipo y `150M` el rango de profundidad configurado (100, 150 o 300 m). El rango se toma del nombre del archivo o de `--prof`; el número de canales (36, 40…) se detecta de las columnas `freqNN`.
- `--freqs`: frecuencias (Hz) de cada columna `freqNN`. Si se omite se asumen log-espaciadas entre 5000 y 12 Hz (`--fmax`, `--fmin`), con `freq01` la más alta (más somera) y `freqNN` la más baja (más profunda).
- `--umbral` (def. 0.1 mV): las lecturas ≤ umbral se consideran ruido/canal muerto y se dejan en ΔVmin (K = 0, azul), como en la pantalla del equipo.
- `--escala-h`: `lineal` (def. con datos del equipo) calcula la profundidad como `prof · canal / n` (freq01 somero … freqNN profundo, igual que el Profile del equipo); `ec12` usa hs = c·503·√(ρ/f) del artículo.
- `--dx` (def. 1.5 m entre puntos; 22 electrodos dan 18 puntos N), `--y0`, `--y-es-n`: posición de los puntos N; `--linea`: filtra por registro L; `--zk`: marca un sondeo.

`python tefsm_figuras.py demo --salida ejemplos` genera CSV sintéticos de prueba (no son datos del artículo).

## Figuras 3 y 4 (perfil vertical y horizontal de la EPD)

Basadas en las figuras 3 y 4 de Gomo y Ngobe, *Telluric Electric Frequency Selection Method (TEFSM) in Geophysical Groundwater Exploration: Emerging Issues*, en *Aquifers – Advances in Hydrogeology* (IntechOpen, https://doi.org/10.5772/intechopen.1013979). Usan los datos del equipo (formato `L, N, freqNN`); la profundidad de cada canal es lineal hasta el rango configurado (100, 150 o 300 m).

```
# Fig. 3: EPD (mV) contra profundidad en una estacion, con litologia y venas de agua opcionales
python tefsm_figuras.py vertical ejemplos/150M_L89.csv --punto 85 --agua 25,45 --litologia ejemplos/litologia_ejemplo.csv

# Fig. 4: EPD (escala log) contra distancia horizontal, una curva por profundidad
python tefsm_figuras.py horizontal ejemplos/150M_L89.csv --zona 9,15
```

- `vertical`: `--punto N` (uno o varios N separados por coma), `--agua` (profundidades de venas de agua), `--litologia` (CSV `tope_m,base_m,nombre`; `ejemplos/litologia_ejemplo.csv` es un ejemplo inventado, no es la litología del sitio).
- `horizontal`: `--profundidades 10,20,…` o `--cada N` (1 de cada N canales, def. 2); `--zona x1,x2[,ymin,ymax]` marca con un recuadro, por ejemplo una zona afectada por ruido. Las lecturas ≤ 0 no se dibujan (escala logarítmica).
- Según ese capítulo, la EPD que mide el equipo no se convierte en resistividad; refleja el efecto de apantallamiento del subsuelo.

## Figura resumen (curvas, sección 2D, inversión y resistividad)

Cuatro paneles a partir de un CSV del equipo: curvas de frecuencia normalizadas con zonas de agua, sección 2D (PowerNorm), "inversión" de la EPD en una estación y modelo de resistividad.

```
python tefsm_figuras.py resumen ejemplos/150M_L89.csv
python tefsm_figuras.py resumen ejemplos/300M_L19.csv --punto 88 --gamma 0.3 --rho-min 5 --rho-max 500
```

**Qué es y qué no es.** Gomo y Ngobe indican que los valores del equipo no se pueden convertir en resistividad. Por eso el panel "modelo de resistividad" es una **aproximación empírica y cualitativa**, no una inversión física:

- *Calculado*: ajuste regularizado (Tikhonov, segunda derivada) de la EPD observada contra la profundidad (`--lam`); *Suavizado*: filtro gaussiano (`--factor-suavizado`). RMS = diferencia entre el ajuste y las lecturas.
- *Resistividad*: mapeo log-lineal de la EPD ajustada al intervalo `--rho-min`…`--rho-max` (más EPD, más resistividad). Calíbrelo con un sondeo o con resistividad eléctrica local (Tabla 1 de Yang et al.: arcilla limosa 20–200, zona fracturada 80–400, granito 130–14 000 Ω·m).
- *Zonas de agua* (cian y recuadros azules): EPD baja respecto a la mediana lateral a cada profundidad (`--umbral-anomalia`, def. 0.15 en log10, ≈ −30 %; `--min-celdas`).
- *Estación*: `--punto N` o `--x m`; por defecto, la de mayor anomalía baja (sin contar los bordes del perfil).
- `--gamma` (def. 0.30) es el exponente de PowerNorm de la sección 2D.

## Figuras de ejemplo

En `ejemplos/figuras/` están los PNG generados con los CSV de `ejemplos/`
(los CSV `datos_fig*.csv` son sintéticos, solo para probar el formato; no reproducen los resultados del artículo):

- `fig2.png`, `fig4.png`, `fig6.png`: módulo de Ey (curvas y pseudo-sección; fig6 incluye la fase).
- `fig13_sintetico.png`: pseudo-sección normalizada con datos sintéticos.
- `resumen_150M_L89.png` y `resumen_300M_L19.png`: figura resumen de cuatro paneles.
- `fig3_vertical_*.png` y `fig4_horizontal_*.png`: figuras 3 y 4 de Gomo y Ngobe con los datos del equipo.
- `fig13_150M_L89.png` y `fig13_300M_L19.png`: pseudo-sección con los datos reales del equipo (sondeos a 150 m y 300 m; 18 puntos cada 1.5 m, profundidad lineal con el canal).

Para regenerarlas: `python tefsm_figuras.py modulo ejemplos/datos_fig2.csv --nombre fig2 --salida ejemplos/figuras`
