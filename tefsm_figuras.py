#!/usr/bin/env python3
"""
TEFSM_inv: figuras a partir de lecturas de mV del equipo TEFSM (selector de frecuencias).

Basado en:
  Yang et al., Water 2025, 17, 3314 (figuras 2, 4, 6 y 13), y
  Gomo y Ngobe, Aquifers - Advances in Hydrogeology, IntechOpen (figuras 3 y 4).

Comandos (datos del equipo: columnas L, N, freq01..freqNN en mV; ver datos/):
  pseudo      Figura 13: pseudo-seccion normalizada log10(dV/dVmin) vs profundidad
  vertical    Fig. 3 (Gomo y Ngobe): EPD vs profundidad en una estacion
  horizontal  Fig. 4 (Gomo y Ngobe): EPD (log) vs distancia, una curva por profundidad
  inversion   SP Inversion y Modelo de resistividad (aproximado), figuras independientes
  modulo      Figuras 2, 4 y 6 de Yang et al.: |Ey| y fase (formato largo y_m, f_hz, ey_mv[, fase_deg])

Ejemplos:
  python tefsm_figuras.py pseudo     datos/150M_L89.csv
  python tefsm_figuras.py vertical   datos/150M_L89.csv --punto 85
  python tefsm_figuras.py horizontal datos/150M_L89.csv
  python tefsm_figuras.py inversion  datos/150M_L89.csv --punto 88

El rango de profundidad (100, 150, 300 o 500 m) se toma del inicio del nombre del archivo (150M_L89.csv) o de --prof.
Las figuras se guardan en figuras/ (--salida). Use -h en cada comando para ver sus opciones.
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# ----------------------------------------------------------------------------
# Lectura de datos
# ----------------------------------------------------------------------------
ALIAS = {
    "y_m": ["y_m", "y", "x", "pos", "posicion", "distancia"],
    "f_hz": ["f_hz", "f", "freq", "frecuencia", "frequency"],
    "ey_mv": ["ey_mv", "ey", "e", "modulo", "mv", "mv_m", "ey_mv_m"],
    "fase_deg": ["fase_deg", "fase", "phase", "phase_deg"],
    "dv_mv": ["dv_mv", "dv", "v", "delta_v", "deltav", "v_mv"],
    "rho_ohm_m": ["rho_ohm_m", "rho", "resistividad"],
}


def leer(path, requeridas, opcionales=()):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep=None, engine="python", decimal=".")
    df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
    out = {}
    for canon in list(requeridas) + list(opcionales):
        for a in ALIAS[canon]:
            if a in df.columns:
                out[canon] = pd.to_numeric(df[a], errors="coerce")
                break
        else:
            if canon in requeridas:
                sys.exit(f"Falta la columna '{canon}' (alias aceptados: {ALIAS[canon]}). "
                         f"Columnas encontradas: {list(df.columns)}")
    res = pd.DataFrame(out).dropna(subset=list(requeridas))
    return res


def leer_equipo(path, args):
    """Formato del equipo: L, N, freq01..freqNN (mV). Devuelve formato largo y_m, f_hz, dv_mv."""
    ext = os.path.splitext(path)[1].lower()
    df = pd.read_excel(path) if ext in (".xlsx", ".xls") else pd.read_csv(path, sep=None, engine="python")
    df.columns = [str(c).strip().lower() for c in df.columns]
    cols = sorted(c for c in df.columns if c.startswith("freq") and c[4:].isdigit())
    if args.linea is not None:
        df = df[df["l"] == args.linea]
    nf = len(cols)
    if args.freqs:
        if os.path.isfile(args.freqs):
            fr = np.loadtxt(args.freqs, delimiter=",").ravel()
        else:
            fr = np.array([float(x) for x in args.freqs.split(",")])
        if len(fr) != nf:
            sys.exit(f"Se dieron {len(fr)} frecuencias pero el archivo tiene {nf} columnas freqNN")
    else:
        fr = np.logspace(np.log10(args.fmax), np.log10(args.fmin), nf)
        if args.escala_h == "ec12":
            print(f"AVISO: no se dieron frecuencias; se asumen {nf} valores log-espaciados "
                  f"entre {args.fmin:g} y {args.fmax:g} Hz (freq01 = la mayor, la mas somera). Use --freqs.")
    y = (df["n"].astype(float) - df["n"].astype(float).min()) * args.dx + args.y0
    if args.y_es_n:
        y = df["n"].astype(float)
    filas = []
    for k, (c, f) in enumerate(zip(cols, fr), start=1):
        filas.append(pd.DataFrame({"y_m": y.values, "f_hz": f, "canal": k,
                                   "dv_mv": pd.to_numeric(df[c], errors="coerce").values}))
    out = pd.concat(filas, ignore_index=True)
    out = out.dropna()
    valido = out["dv_mv"] > args.umbral
    # lecturas de ruido/canal muerto: se dejan en el minimo valido (K = 0), como en la pantalla del equipo
    out.loc[~valido, "dv_mv"] = out.loc[valido, "dv_mv"].min()
    out.attrs["fmin_equipo"] = float(np.min(fr))
    out.attrs["n_canales"] = nf
    return out


def a_malla(df, valor):
    """Pivota a matriz [frecuencias x posiciones]; promedia duplicados."""
    p = df.pivot_table(index="f_hz", columns="y_m", values=valor, aggfunc="mean")
    p = p.sort_index().sort_index(axis=1)
    return p.columns.values.astype(float), p.index.values.astype(float), p.values


# ----------------------------------------------------------------------------
# Paleta tipo "jet" usada en las figuras del articulo
# ----------------------------------------------------------------------------
CMAP_E = plt.get_cmap("jet")
CMAP_K = LinearSegmentedColormap.from_list(
    "k13", ["#0000ff", "#00b4ff", "#a8ffd8", "#ffff80", "#ffb000", "#ff3000", "#ff0000"])

# Frecuencias mostradas en el panel (a) de las figuras 2, 4 y 6: 10^x Hz
FREC_CURVAS = [1.1, 1.4, 1.6, 1.85, 2.0, 2.25]
ESTILOS = {1.1: (":", None), 1.4: ("-", None), 1.6: (":", "^"),
           1.85: (":", None), 2.0: (":", "x"), 2.25: ("-", "^")}


# ----------------------------------------------------------------------------
# Figuras 2, 4, 6
# ----------------------------------------------------------------------------
PROFUNDIDADES = (100, 150, 300, 500)     # rangos de profundidad del equipo (m)


def profundidad_equipo(path, prof=None):
    """Rango de profundidad (100, 150, 300 o 500 m) indicado al inicio del nombre: '150M_L89.csv' -> 150."""
    if prof is None:
        m = re.match(r"(\d+)M_", os.path.basename(path), re.I)
        prof = float(m.group(1)) if m else None
    if prof not in PROFUNDIDADES:
        sys.exit(f"No se pudo determinar la profundidad de '{os.path.basename(path)}'. El nombre debe empezar con "
                 f"100M_, 150M_, 300M_ o 500M_ (p. ej. 150M_L89.csv), o indique --prof "
                 f"({', '.join(map(str, PROFUNDIDADES))}).")
    return float(prof)


def carpeta_salida(path):
    """'150M_L89.csv' -> 'figuras/L89_150m' (una carpeta por archivo de entrada)."""
    m = re.match(r"(\d+)\s*M_L(\d+)", os.path.splitext(os.path.basename(path))[0], re.I)
    return os.path.join("figuras", f"L{m.group(2)}_{m.group(1)}m") if m else \
        os.path.join("figuras", os.path.splitext(os.path.basename(path))[0])


def etiqueta_csv(path):
    """'150M_L89.csv' -> 'L89 a 150m' (registro y rango de profundidad del nombre del archivo)."""
    nombre = os.path.basename(path)
    m = re.match(r"(\d+)\s*M_L(\d+)", os.path.splitext(nombre)[0], re.I)
    return f"L{m.group(2)} a {m.group(1)}m" if m else nombre


def figura_modulo(args):
    df = leer(args.archivo, ["y_m", "f_hz", "ey_mv"], ["fase_deg"])
    y, f, E = a_malla(df, "ey_mv")
    lgf = np.log10(f)
    paneles = 3 if "fase_deg" in df.columns and df["fase_deg"].notna().any() else 2

    fig, axs = plt.subplots(1, paneles, figsize=(5.2 * paneles, 5.2))

    # (a) curvas de modulo vs y
    ax = axs[0]
    for lg in FREC_CURVAS:
        i = int(np.argmin(np.abs(lgf - lg)))
        if abs(lgf[i] - lg) > 0.06:  # no hay dato cercano a esa frecuencia
            continue
        ls, mk = ESTILOS[lg]
        ax.plot(y, E[i], color="k", lw=0.9, ls=ls, marker=mk, ms=4, mfc="k",
                markevery=max(1, len(y) // 20))
        etiqueta = f"$10^{{{lg:g}}}$ Hz"
        ax.annotate(etiqueta, (y[-1], E[i][-1]), xytext=(-4, -12),
                    textcoords="offset points", ha="right", fontsize=8)
    ax.set_xlabel("x / m")
    ax.set_ylabel(r"$|E_y|$ / mV·m$^{-1}$")
    ax.set_xlim(y.min(), y.max())
    ax.set_title("(a)", y=-0.2)

    # (b) pseudo-seccion del modulo
    niveles = [0, .5, .75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
               15, 16, 17, 18, 19, max(20, float(np.nanmax(E)))]
    ax = axs[1]
    cf = ax.contourf(y, lgf, E, levels=niveles, cmap=CMAP_E, extend="max")
    cs = ax.contour(y, lgf, E, levels=[.75, 1, 1.5, 2, 3, 5, 10, 15],
                    colors="k", linewidths=0.6)
    ax.clabel(cs, fmt="%g", fontsize=7)
    ax.set_xlabel("x / m")
    ax.set_ylabel("lg f / Hz")
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_title("(b)", y=-0.2)
    fig.colorbar(cf, ax=ax, label=r"$|E_y|$ / mV·m$^{-1}$", ticks=niveles[:-1], pad=0.03)

    # (c) pseudo-seccion de la fase
    if paneles == 3:
        _, _, Ph = a_malla(df, "fase_deg")
        ax = axs[2]
        cf = ax.contourf(y, lgf, Ph, levels=20, cmap=CMAP_E)
        cs = ax.contour(y, lgf, Ph, levels=10, colors="k", linewidths=0.5)
        ax.clabel(cs, fmt="%g", fontsize=6)
        ax.set_xlabel("x / m")
        ax.set_ylabel("lg f / Hz")
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_title("(c)", y=-0.2)
        fig.colorbar(cf, ax=ax, label="φ / (°)", pad=0.03)

    fig.suptitle(f"Módulo de Ey - {etiqueta_csv(args.archivo)}", fontweight="bold")
    guardar(fig, args)


# ----------------------------------------------------------------------------
# Figura 13: pseudo-seccion normalizada
#   K_i = log10(dV_i / dV_min)        (ec. 13)
#   h_s = c * 503 * sqrt(rho / f)     (ec. 12)
# ----------------------------------------------------------------------------
def figura_pseudo(args):
    cab = pd.read_csv(args.archivo, nrows=0, sep=None, engine="python") \
        if not args.archivo.lower().endswith((".xlsx", ".xls")) else pd.read_excel(args.archivo, nrows=0)
    if any(str(c).strip().lower().startswith("freq") for c in cab.columns):
        df = leer_equipo(args.archivo, args)
    else:
        df = leer(args.archivo, ["y_m", "f_hz", "dv_mv"], ["rho_ohm_m"])
    if "rho_ohm_m" not in df.columns:
        df["rho_ohm_m"] = args.rho
    df["rho_ohm_m"] = df["rho_ohm_m"].fillna(args.rho)
    df = df[df["dv_mv"] > 0]

    dvmin = df["dv_mv"].min()  # minimo de todo el perfil
    df["K"] = np.log10(df["dv_mv"] / dvmin)
    # Rango de profundidad del equipo (100, 150, 300 o 500 m): del nombre ("150M_L89.csv") o --prof
    prof_eq = profundidad_equipo(args.archivo, args.prof) if "canal" in df.columns else args.prof
    c = args.c
    if c is None:
        if prof_eq is not None:
            # c tal que la frecuencia mas baja del equipo llegue a la profundidad configurada
            fmin = df.attrs.get("fmin_equipo", df["f_hz"].min())
            c = prof_eq / (503.0 * np.sqrt(args.rho / fmin))
            if args.escala_h == "ec12":
                print(f"Rango del equipo {prof_eq:g} m -> c = {c:.4f} (calibrado con f_min = {fmin:g} Hz, rho = {args.rho:g})")
        else:
            c = 1.0
    if args.escala_h == "lineal" and "canal" in df.columns and prof_eq is not None:
        # igual que el Profile del equipo: h = prof * canal / n_canales (freq01 somero ... freqNN profundo)
        df["hs"] = prof_eq * df["canal"] / (df.attrs.get("n_canales", df["canal"].max()) + args.offset_canales)
    else:
        df["hs"] = c * 503.0 * np.sqrt(df["rho_ohm_m"] / df["f_hz"])

    # interpolar K(hs) en cada posicion y sobre una malla de profundidad comun
    ys = np.sort(df["y_m"].unique())
    hmax = args.hmax or prof_eq or float(df["hs"].max())
    prof = np.linspace(0.0, hmax, 200)
    K = np.full((len(prof), len(ys)), np.nan)
    for j, yy in enumerate(ys):
        s = df[df["y_m"] == yy].sort_values("hs")
        if len(s) >= 2:
            K[:, j] = np.interp(prof, s["hs"], s["K"], right=np.nan)

    fig, ax = plt.subplots(figsize=(7.5, 6))
    tope = max(2.0, np.ceil(np.nanmax(K) * 10) / 10)
    niveles = np.round(np.arange(0, tope + 1e-9, 0.1), 1)
    cf = ax.contourf(ys, -prof, K, levels=niveles, cmap=CMAP_K, extend="neither")
    cs = ax.contour(ys, -prof, K, levels=niveles[1:-1], colors="k", linewidths=0.4)
    ax.clabel(cs, fmt="%.1f", fontsize=6)
    if args.zk is not None:
        ax.axvline(args.zk, color="red", lw=1.5)
        ax.annotate("ZK", (args.zk, 0), xytext=(0, 22), textcoords="offset points",
                    color="red", ha="center", weight="bold")
    ax.set_xlabel("x / m")
    ax.set_ylabel(r"$h_s$ / m")
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xticks(ys)
    ax.set_xticklabels([f"{v:g}" for v in ys], fontsize=8)
    ax.set_ylim(-hmax, 0)
    fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.04, shrink=0.8,
                 label=r"$\log_{10}(\Delta V/\Delta V_{min})$")
    ax.set_title(f"Pseudo-sección normalizada\n{etiqueta_csv(args.archivo)}", fontweight="bold", pad=42)
    guardar(fig, args)


# ----------------------------------------------------------------------------
# Figuras 3 y 4 de: Gomo & Ngobe, "Telluric Electric Frequency Selection Method (TEFSM) in
# Geophysical Groundwater Exploration: Emerging Issues", en Aquifers - Advances in
# Hydrogeology, IntechOpen, DOI 10.5772/intechopen.1013979
#   Fig. 3: perfil VERTICAL de la EPD (mV) contra la profundidad en una estacion,
#           con columna litologica y venas de agua opcionales.
#   Fig. 4: perfil HORIZONTAL de la EPD (escala log) con una curva por profundidad.
# Entrada: formato del equipo (L, N, freqNN), igual que la figura 13.
# ----------------------------------------------------------------------------
def matriz_equipo(args):
    """Devuelve y (m, por punto N), h (m, por canal) y V[canal, punto] en mV, sin descartar nada."""
    ext = os.path.splitext(args.archivo)[1].lower()
    df = pd.read_excel(args.archivo) if ext in (".xlsx", ".xls") else pd.read_csv(args.archivo, sep=None, engine="python")
    df.columns = [str(c).strip().lower() for c in df.columns]
    cols = sorted(c for c in df.columns if c.startswith("freq") and c[4:].isdigit())
    if not cols:
        sys.exit("Estas figuras requieren el formato del equipo (columnas L, N, freq01..freqNN)")
    if args.linea is not None:
        df = df[df["l"] == args.linea]
    df = df.sort_values("n")
    n = df["n"].astype(float).values
    y = (n - n.min()) * args.dx + args.y0
    V = df[cols].apply(pd.to_numeric, errors="coerce").values.T
    nf = V.shape[0]
    prof = profundidad_equipo(args.archivo, args.prof)
    h = prof * np.arange(1, nf + 1) / (nf + args.offset_canales)   # profundidad lineal con el canal
    args.prof = prof
    return n, y, h, V


def estacion_auto(V, umbral=0.1, umbral_anomalia=0.15):
    """Indice del punto con mayor anomalia de EPD BAJA respecto a la mediana lateral (sin los bordes)."""
    import warnings
    npts = V.shape[1]
    valido = V > umbral
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(valido, np.log10(np.where(valido, V, 1.0)), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(L, axis=1, keepdims=True)       # canales sin datos validos -> NaN
    A = L - med
    bajo = np.where(np.isnan(A), 0.0, np.minimum(A + umbral_anomalia, 0.0))
    score = bajo.sum(axis=0) / np.maximum(valido.sum(axis=0), 1)    # < 0 = zona baja
    lo, hi = (1, npts - 1) if npts > 4 else (0, npts)
    return lo + int(np.argmin(score[lo:hi]))


def figura_vertical(args):
    n, y, h, V = matriz_equipo(args)
    puntos = [float(p) for p in args.punto.split(",")] if args.punto else [n[estacion_auto(V, args.umbral, args.umbral_anomalia)]]
    fig = plt.figure(figsize=(7.5, 8))
    litologia = None
    if args.litologia:
        litologia = pd.read_csv(args.litologia)
        litologia.columns = [c.strip().lower() for c in litologia.columns]   # tope_m, base_m, nombre
        gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.02)
        ax = fig.add_subplot(gs[0])
        axl = fig.add_subplot(gs[1], sharey=ax)
    else:
        ax = fig.add_subplot(111)
    for p in puntos:
        j = int(np.argmin(np.abs(n - p)))
        ax.plot(V[:, j], h, color="#e8603c", lw=1.4, marker="D", ms=4, mfc="#c8102e",
                mec="k", mew=0.6, label=f"N = {n[j]:g}  (y = {y[j]:g} m)")
    for a in (args.agua.split(",") if args.agua else []):
        ax.axhline(float(a), color="k", ls="--", lw=1.4)
    ax.set_ylim(h.max() * 1.02, 0)
    ax.set_xlim(0, np.nanmax(V[:, [int(np.argmin(np.abs(n - p))) for p in puntos]]) * 1.1)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Diferencia de potencial eléctrico (mV)")
    ax.set_ylabel("Profundidad bajo la superficie (m)")
    ax.grid(color="#bbb", lw=0.4)
    ax.grid(which="major", color="#555", lw=0.6)
    ax.minorticks_on()
    if len(puntos) > 1:
        ax.legend(loc="lower right", fontsize=8)
    if litologia is not None:
        axl.set_xlim(0, 1)
        paleta = {}
        for i, r in enumerate(litologia.itertuples()):
            col = paleta.setdefault(r.nombre, plt.get_cmap("Set2")(len(paleta) % 8))
            axl.add_patch(plt.Rectangle((0, r.tope_m), 1, r.base_m - r.tope_m, fc=col, ec="k", lw=0.6,
                                        hatch=["", "..", "////", "xx", "\\\\"][list(paleta).index(r.nombre) % 5]))
        axl.set_title("Litología", fontsize=9)
        axl.set_xticks([])
        axl.yaxis.tick_right()
        axl.tick_params(labelright=False)
        axl.xaxis.tick_top()
        handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="k") for c in paleta.values()]
        fig.legend(handles, list(paleta), loc="lower center", ncol=len(paleta), fontsize=8, frameon=False)
    est = ""
    if len(puntos) == 1:
        jj = int(np.argmin(np.abs(n - puntos[0])))
        est = f"\nx = {y[jj]:g} m (N = {n[jj]:g})"
    fig.suptitle(f"Perfil vertical de EPD\n{etiqueta_csv(args.archivo)}{est}", fontweight="bold")
    guardar(fig, args)


def figura_horizontal(args):
    n, y, h, V = matriz_equipo(args)
    if args.profundidades:
        idx = [int(np.argmin(np.abs(h - float(d)))) for d in args.profundidades.split(",")]
    else:
        cada = args.cada or max(1, round(V.shape[0] / 18))     # unas 18 curvas, como la figura original
        idx = list(range(cada - 1, V.shape[0], cada))
    fig, ax = plt.subplots(figsize=(10, 6))
    cm = plt.get_cmap("tab20")
    for k, i in enumerate(idx):
        v = np.where(V[i] > 0, V[i], np.nan)           # la escala log no admite ceros
        ax.plot(y, v, lw=1.2, color=cm(k % 20), label=f"{h[i]:.0f} m")
    ax.set_yscale("log")
    pos = V[V > 0]
    ax.set_ylim(10 ** np.floor(np.log10(pos.min())), pos.max() * 1.5)
    ax.set_xlim(y.min(), y.max())
    ax.set_xlabel("Distancia horizontal (m)")
    ax.set_ylabel("Diferencia de potencial eléctrico (mV)")
    ax.grid(which="both", color="#ccc", lw=0.4)
    if args.zona:
        z = [float(t) for t in args.zona.split(",")]
        y0, y1 = ax.get_ylim()
        zy0, zy1 = (z[2], z[3]) if len(z) == 4 else (y0, y1)
        ax.add_patch(plt.Rectangle((z[0], zy0), z[1] - z[0], zy1 - zy0, fill=False, ls="--", ec="k", lw=1.4))
    ax.legend(ncol=3, fontsize=8, loc="lower right", framealpha=0.9)
    ax.set_title(f"Perfil horizontal de EPD\n{etiqueta_csv(args.archivo)}", fontweight="bold")
    guardar(fig, args)


# ----------------------------------------------------------------------------
# "SP Inversion" y "Modelo de resistividad (aproximado)": dos figuras independientes
# para una estacion del perfil.
#
# IMPORTANTE: la EPD (mV) del equipo NO es una resistividad (Gomo y Ngobe: "groundwater detector
# values cannot be converted to resistivity values"). Por eso aqui:
#   * "Calculado" = ajuste regularizado (Tikhonov, 2a derivada) de la EPD observada vs profundidad.
#   * Resistividad = mapeo EMPIRICO log-lineal de la EPD ajustada al intervalo [rho_min, rho_max]
#     (mas EPD => mas resistividad). Es un modelo cualitativo; calibrelo con un sondeo o con
#     resistividad electrica local (Tabla 1 del articulo de Yang et al.: arcilla limosa 20-200,
#     zona fracturada 80-400, granito 130-14000 ohm.m).
#   * Estacion por defecto: la de mayor anomalia de EPD BAJA respecto a la mediana lateral a cada
#     profundidad (sin contar los bordes del perfil).
# ----------------------------------------------------------------------------
def figura_inversion(args):
    from scipy import ndimage

    n, x, h, V = matriz_equipo(args)
    valido = V > args.umbral
    if args.x is not None:
        j = int(np.argmin(np.abs(x - args.x)))
    elif args.punto is not None:
        j = int(np.argmin(np.abs(n - args.punto)))
    else:
        j = estacion_auto(V, args.umbral, args.umbral_anomalia)
    xs = x[j]

    obs = valido[:, j]
    hv, vv = h[obs], V[obs, j]
    if obs.sum() < 4:
        sys.exit(f"La estacion N = {n[j]:g} tiene menos de 4 lecturas validas (> {args.umbral:g} mV)")
    D = np.diff(np.eye(len(vv)), 2, axis=0)
    sc = vv.max()
    m = np.linalg.solve(np.eye(len(vv)) + args.lam * D.T @ D, vv / sc) * sc    # ajuste regularizado
    suav = ndimage.gaussian_filter1d(vv, args.factor_suavizado, mode="nearest")
    rms = float(np.sqrt(np.mean((m - vv) ** 2)))

    base = os.path.splitext(os.path.basename(args.archivo))[0]
    os.makedirs(args.salida, exist_ok=True)

    # --- SP Inversion
    fig, ax = plt.subplots(figsize=(5, 7))
    ax.plot(vv, hv, ".", color="#6a5acd", ms=6, alpha=0.7, label="Observado")
    ax.plot(suav, hv, color="#17becf", lw=1.3, label=f"Suavizado (factor {args.factor_suavizado:g})")
    ax.plot(m, hv, color="red", lw=3, label=f"Calculado\nRMS={rms:.2f} ({100 * rms / sc:.1f}%)")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_ylim(args.prof, 0)
    ax.set_xlim(left=0)
    ax.set_xlabel("SP (mV)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"SP Inversion\n{etiqueta_csv(args.archivo)}\nx = {xs:g} m (N = {n[j]:g})", fontweight="bold")
    ax.grid(ls=":", alpha=0.5)
    fig.tight_layout()
    ruta = os.path.join(args.salida, f"SP_Inversion_{base}.png")
    terminar(fig, ruta, args, f"(estacion x = {xs:g} m, N = {n[j]:g})")

    # --- Modelo de resistividad (aproximado)
    t = np.clip((m - m.min()) / max(m.max() - m.min(), 1e-12), 0, 1)
    rho = args.rho_min * (args.rho_max / args.rho_min) ** t
    paso = h[1] - h[0]
    bordes = np.concatenate([hv - paso / 2, [hv[-1] + paso / 2]])
    fig, ax = plt.subplots(figsize=(5, 7))
    ax.step(np.concatenate([rho, rho[-1:]]), bordes, where="post", color="k", lw=1.8)
    ext = [0] + [i for i in range(1, len(rho) - 1) if (rho[i] - rho[i - 1]) * (rho[i + 1] - rho[i]) < 0] + [len(rho) - 1]
    for i in sorted(set(ext)):
        ax.annotate(f"{rho[i]:.3g}", (rho[i], hv[i]), fontsize=7.5, fontweight="bold",
                    bbox=dict(boxstyle="square,pad=0.15", fc="#f5e663", ec="k", lw=0.6),
                    xytext=(3, 0), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_ylim(args.prof, 0)
    ax.set_xlabel("Resistividad (ohm·m)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"Modelo de resistividad (aproximado)\n{etiqueta_csv(args.archivo)}\nx = {xs:g} m (N = {n[j]:g})", fontweight="bold")
    ax.grid(ls=":", alpha=0.5, which="both")
    fig.text(0.5, 0.005, f"Mapeo empírico de la EPD a {args.rho_min:g}–{args.rho_max:g} ohm·m; no es una inversión física",
             ha="center", fontsize=7, color="#555")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    ruta = os.path.join(args.salida, f"Modelo_Resistividad_{base}.png")
    terminar(fig, ruta, args)


def terminar(fig, ruta, args, extra=""):
    """Guarda el PNG, imprime su ruta completa y, con --mostrar, abre la figura en pantalla."""
    fig.savefig(ruta, dpi=args.dpi)
    print("Figura guardada en", os.path.abspath(ruta), extra)
    if args.mostrar:
        plt.show()
    plt.close(fig)


def guardar(fig, args):
    os.makedirs(args.salida, exist_ok=True)
    ruta = os.path.join(args.salida, f"{args.nombre}.png")
    fig.tight_layout()
    terminar(fig, ruta, args)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def comun(sp):
        sp.add_argument("--nombre", default=None, help="nombre del PNG de salida (sin extension)")
        sp.add_argument("--salida", default=None, help="carpeta de salida (def.: figuras/L89_150m, segun el archivo)")
        sp.add_argument("--dpi", type=int, default=200)
        sp.add_argument("--mostrar", action="store_true", help="abrir las figuras en una ventana (ademas de guardarlas)")

    a = sub.add_parser("modulo", help="Figuras 2, 4 y 6 (|Ey| y fase)")
    a.add_argument("archivo")
    comun(a)
    a.set_defaults(fn=figura_modulo)

    b = sub.add_parser("pseudo", help="Figura 13 (pseudo-seccion normalizada de dV)")
    b.add_argument("archivo")
    b.add_argument("--rho", type=float, default=220.0,
                   help="resistividad aparente (ohm.m) si no hay columna rho (def. 220, zona fracturada)")
    b.add_argument("--c", type=float, default=None,
                   help="coeficiente empirico c de la ec. 12 (si falta: se calibra con --prof, o 1)")
    b.add_argument("--prof", type=float, default=None, choices=PROFUNDIDADES,
                   help="rango de profundidad del equipo (m); por defecto se lee del nombre del archivo (150M_...)")
    b.add_argument("--escala-h", choices=["lineal", "ec12"], default="lineal",
                   help="(formato equipo) profundidad: 'lineal' con el canal, como el Profile del equipo (def.), "
                        "o 'ec12' = c*503*sqrt(rho/f) del articulo")
    b.add_argument("--offset-canales", dest="offset_canales", type=int, default=0,
                   help="profundidad = prof * canal / (n_canales + offset); 0 (def.) o 1 segun la convencion del software")
    b.add_argument("--hmax", type=float, default=None, help="profundidad maxima mostrada (m)")
    b.add_argument("--zk", type=float, default=None, help="posicion (m) del sondeo ZK a marcar")
    b.add_argument("--freqs", default=None,
                   help="(formato equipo) frecuencias en Hz separadas por coma, o archivo con ellas; "
                        "una por columna freqNN, en el mismo orden")
    b.add_argument("--fmin", type=float, default=12.0, help="(formato equipo) fmin si no hay --freqs")
    b.add_argument("--fmax", type=float, default=5000.0, help="(formato equipo) fmax si no hay --freqs")
    b.add_argument("--dx", type=float, default=1.0, help="(formato equipo) metros entre electrodos / puntos N consecutivos (def. 1)")
    b.add_argument("--y0", type=float, default=0.0, help="(formato equipo) posicion y (m) del menor N")
    b.add_argument("--y-es-n", action="store_true", help="(formato equipo) usar N directamente como y (m)")
    b.add_argument("--linea", type=int, default=None, help="(formato equipo) filtrar por el registro L del equipo (p. ej. 93)")
    b.add_argument("--umbral", type=float, default=0.1,
                   help="descarta lecturas <= umbral (mV) como ruido/canal muerto (def. 0.1); "
                        "evita que dV_min sea ~0 en la ec. 13")
    comun(b)
    b.set_defaults(fn=figura_pseudo)

    def comun_eq(sp):
        sp.add_argument("archivo")
        sp.add_argument("--prof", type=float, default=None, choices=PROFUNDIDADES,
                        help="rango de profundidad del equipo (m); por defecto se lee del nombre del archivo (150M_...)")
        sp.add_argument("--dx", type=float, default=1.0, help="metros entre electrodos / puntos N consecutivos (def. 1)")
        sp.add_argument("--y0", type=float, default=0.0, help="posicion y (m) del menor N")
        sp.add_argument("--linea", type=int, default=None, help="filtrar por registro L del equipo")
        sp.add_argument("--offset-canales", dest="offset_canales", type=int, default=0,
                        help="profundidad = prof * canal / (n_canales + offset); 0 (def.) o 1 segun la convencion del software")
        sp.add_argument("--umbral", type=float, default=0.1, help="lecturas <= umbral (mV) = ruido/canal muerto (def. 0.1)")
        sp.add_argument("--umbral-anomalia", dest="umbral_anomalia", type=float, default=0.15,
                        help="caida (log10) bajo la mediana lateral para elegir la estacion (def. 0.15 = -30%%)")
        comun(sp)

    v = sub.add_parser("vertical", help="Fig. 3 (Gomo y Ngobe): perfil vertical de EPD vs profundidad")
    comun_eq(v)
    v.add_argument("--punto", default=None, help="N de la estacion (o varios separados por coma); def.: la de mayor anomalia baja")
    v.add_argument("--litologia", default=None, help="CSV con columnas tope_m, base_m, nombre")
    v.add_argument("--agua", default=None, help="profundidades (m) de venas de agua, separadas por coma")
    v.set_defaults(fn=figura_vertical)

    hz = sub.add_parser("horizontal", help="Fig. 4 (Gomo y Ngobe): perfil horizontal de EPD (log), una curva por profundidad")
    comun_eq(hz)
    hz.add_argument("--profundidades", default=None, help="profundidades (m) a graficar, separadas por coma")
    hz.add_argument("--cada", type=int, default=None, help="si no hay --profundidades: 1 de cada N canales (def.: el necesario para ~18 curvas)")
    hz.add_argument("--zona", default=None, help="recuadro a marcar: x1,x2[,ymin,ymax] (p. ej. zona de ruido)")
    hz.set_defaults(fn=figura_horizontal)

    rs = sub.add_parser("inversion", help="SP Inversion y Modelo de resistividad (aproximado), dos figuras independientes")
    comun_eq(rs)
    rs.add_argument("--punto", type=float, default=None, help="N de la estacion (def.: la de mayor anomalia baja)")
    rs.add_argument("--x", type=float, default=None, help="posicion x (m) de la estacion (alternativa a --punto)")
    rs.add_argument("--lam", type=float, default=3.0, help="regularizacion del ajuste (def. 3)")
    rs.add_argument("--factor-suavizado", dest="factor_suavizado", type=float, default=3.0, help="sigma del suavizado (canales)")
    rs.add_argument("--rho-min", dest="rho_min", type=float, default=5.0, help="resistividad asignada a la EPD minima (ohm.m)")
    rs.add_argument("--rho-max", dest="rho_max", type=float, default=500.0, help="resistividad asignada a la EPD maxima (ohm.m)")
    rs.set_defaults(fn=figura_inversion)

    def completar(a):
        if not a.mostrar:
            plt.switch_backend("Agg")        # sin ventana: solo se guardan los PNG
        if a.salida is None:
            a.salida = carpeta_salida(a.archivo)
        if a.nombre is None:                 # nombre por defecto distinto para cada comando
            prefijo = {"pseudo": "fig13_", "vertical": "fig3_vertical_", "horizontal": "fig4_horizontal_",
                       "modulo": "modulo_"}.get(a.cmd, "")
            a.nombre = prefijo + os.path.splitext(os.path.basename(a.archivo))[0]

    def figura_todas(a):
        """Genera todas las figuras de un archivo del equipo en su carpeta (figuras/L89_150m)."""
        base = [a.archivo, "--salida", a.salida, "--dpi", str(a.dpi), "--dx", str(a.dx), "--y0", str(a.y0),
                "--umbral", str(a.umbral), "--offset-canales", str(a.offset_canales)]
        if a.prof:
            base += ["--prof", str(a.prof)]
        if a.linea is not None:
            base += ["--linea", str(a.linea)]
        if a.mostrar:
            base += ["--mostrar"]
        for c in ("pseudo", "vertical", "horizontal", "inversion"):
            extra = ["--punto", str(a.punto)] if (a.punto is not None and c in ("vertical", "inversion")) else []
            if c != "pseudo":                # pseudo no usa la seleccion de estacion
                extra += ["--umbral-anomalia", str(a.umbral_anomalia)]
            sa = p.parse_args([c] + base + extra)
            completar(sa)
            sa.fn(sa)

    t = sub.add_parser("todas", help="todas las figuras de un archivo del equipo, en su carpeta (figuras/L89_150m)")
    comun_eq(t)
    t.add_argument("--punto", type=float, default=None, help="N de la estacion para vertical e inversion (def.: la de mayor anomalia baja)")
    t.set_defaults(fn=figura_todas)

    args = p.parse_args()
    completar(args)
    args.fn(args)


if __name__ == "__main__":
    main()
