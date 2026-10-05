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

Formato del equipo (`L, N, freq01…freqNN`, en mV; ver `ejemplos/150M_L93.csv`) o formato largo `y_m, f_hz, dv_mv [, rho_ohm_m]`:

```
python tefsm_figuras.py pseudo ejemplos/150M_L93.csv --freqs frecuencias.txt --dx 1 --rho 220
```

- K = log10(ΔV/ΔVmin) (ec. 13) y hs = c·503·√(ρ/f) (ec. 12).
- `L93` es el registro de memoria del equipo; `150M` es el rango de profundidad configurado (100, 150 o 300 m). Se toma del nombre del archivo o de `--prof`, y con él se calibra c para que la frecuencia más baja llegue a esa profundidad. Con `--c` se fija manualmente.
- `--freqs`: frecuencias (Hz) de cada columna `freqNN`. Si se omite se asumen log-espaciadas entre 12 y 5000 Hz (`--fmin`, `--fmax`).
- `--umbral` (def. 0.1 mV): descarta lecturas ≤ umbral como ruido/canal muerto.
- `--dx`, `--y0`, `--y-es-n`: posición de los puntos N; `--linea`: filtra por registro L; `--zk`: marca un sondeo.

`python tefsm_figuras.py demo --salida ejemplos` genera CSV sintéticos de prueba (no son datos del artículo).
