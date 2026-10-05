#!/usr/bin/env python3
"""
Reproduce las figuras 2, 4, 6 y 13 de:
  Yang et al., "Simulation of the Telluric Electrical Field Frequency Selection
  Method and Its Application in Mineral Water Exploration", Water 2025, 17, 3314.

Los datos (lecturas en mV) se ingresan desde un archivo CSV o Excel (.xlsx).

USO
---
  # Figuras 2, 4 y 6 (modulo de Ey: curvas + pseudo-seccion; fase si existe)
  python tefsm_figuras.py modulo  datos_fig2.csv --nombre fig2
  python tefsm_figuras.py modulo  datos_fig4.csv --nombre fig4
  python tefsm_figuras.py modulo  datos_fig6.csv --nombre fig6   # con columna fase_deg -> panel (c)

  # Figura 13 (pseudo-seccion normalizada de dV, 40 frecuencias)
  python tefsm_figuras.py pseudo  datos_L8.csv --nombre fig13b --zk 23 --c 0.1 --rho 220

  # Figuras 3 y 4 de Gomo y Ngobe (perfil vertical y horizontal de la EPD, datos del equipo)
  python tefsm_figuras.py vertical   150M_L89.csv --punto 85 --agua 25,45
  python tefsm_figuras.py horizontal 150M_L89.csv --zona 10,18

  # Archivos de ejemplo sinteticos para probar el formato
  python tefsm_figuras.py demo --salida ejemplos

FORMATO DE ENTRADA (formato "largo": una fila por lectura)
----------------------------------------------------------
  Figuras 2, 4, 6  -> columnas:  y_m, f_hz, ey_mv  [, fase_deg]
       y_m      posicion sobre el perfil (m)
       f_hz     frecuencia (Hz)
       ey_mv    |Ey| en mV/m  (modulo del campo electrico)
       fase_deg fase de Ey en grados (opcional; solo para el panel c de la fig. 6)

  Figura 13        -> columnas:  y_m, f_hz, dv_mv  [, rho_ohm_m]
       dv_mv      diferencia de potencial dV medida (mV)
       rho_ohm_m  resistividad aparente (opcional; si falta se usa --rho)

  Tambien se aceptan los alias: y/x/pos, f/freq/frecuencia, ey/e/modulo/mv,
  phase/fase, dv/v/delta_v, rho/resistividad.
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
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
    # Rango de profundidad del equipo (100/150/300 m): --prof, o se toma del nombre ("150M_L89.csv")
    prof_eq = args.prof
    if prof_eq is None:
        m = re.search(r"(\d+)\s*m", os.path.basename(args.archivo), re.I)
        prof_eq = float(m.group(1)) if m and float(m.group(1)) in (100, 150, 300) else None
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
        df["hs"] = prof_eq * df["canal"] / df.attrs.get("n_canales", df["canal"].max())
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
    prof = args.prof
    if prof is None:
        m = re.search(r"(\d+)\s*m", os.path.basename(args.archivo), re.I)
        prof = float(m.group(1)) if m else 150.0
    h = prof * np.arange(1, nf + 1) / nf          # profundidad lineal con el canal (como el equipo)
    args.prof = prof
    return n, y, h, V


def figura_vertical(args):
    n, y, h, V = matriz_equipo(args)
    puntos = [float(p) for p in args.punto.split(",")] if args.punto else [n[0]]
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
    args.nombre = args.nombre or "fig3_perfil_vertical"
    guardar(fig, args)


def figura_horizontal(args):
    n, y, h, V = matriz_equipo(args)
    if args.profundidades:
        idx = [int(np.argmin(np.abs(h - float(d)))) for d in args.profundidades.split(",")]
    else:
        idx = list(range(args.cada - 1, V.shape[0], args.cada))
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
    args.nombre = args.nombre or "fig4_perfil_horizontal"
    guardar(fig, args)


# ----------------------------------------------------------------------------
# Figura resumen (4 paneles): curvas normalizadas + zonas, seccion 2D (PowerNorm),
# "inversion" en una estacion y modelo de resistividad.
#
# IMPORTANTE: la EPD (mV) del equipo NO es una resistividad (Gomo y Ngobe: "groundwater detector
# values cannot be converted to resistivity values"). Por eso aqui:
#   * "Calculado" = ajuste regularizado (Tikhonov, 2a derivada) de la EPD observada vs profundidad.
#   * Resistividad = mapeo EMPIRICO log-lineal de la EPD ajustada al intervalo [rho_min, rho_max]
#     (mas EPD => mas resistividad). Es un modelo cualitativo; calibrelo con un sondeo o con
#     resistividad electrica local (Tabla 1 del articulo de Yang et al.: arcilla limosa 20-200,
#     zona fracturada 80-400, granito 130-14000 ohm.m).
#   * Zonas de agua = anomalias de EPD BAJA respecto a la mediana lateral a cada profundidad.
# ----------------------------------------------------------------------------
def figura_resumen(args):
    from scipy import ndimage
    from matplotlib.colors import PowerNorm

    n, x, h, V = matriz_equipo(args)
    nf, npts = V.shape
    valido = V > args.umbral
    Vn = V / np.nanmax(V)                                   # respuesta normalizada 0-1
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(valido, np.log10(np.where(valido, V, 1.0)), np.nan)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(L, axis=1, keepdims=True)       # canales sin datos validos -> NaN
    A = L - med          # anomalia (log10) vs mediana lateral
    umb_a = args.umbral_anomalia
    bajo = np.where(np.isnan(A), 0.0, np.minimum(A + umb_a, 0.0))   # solo caidas > umbral
    score = bajo.sum(axis=0) / np.maximum(valido.sum(axis=0), 1)    # <0 = zona baja
    if args.x is not None:
        j = int(np.argmin(np.abs(x - args.x)))
    elif args.punto is not None:
        j = int(np.argmin(np.abs(n - args.punto)))
    else:
        lo, hi = (1, npts - 1) if npts > 4 else (0, npts)      # se evitan los bordes del perfil
        j = lo + int(np.argmin(score[lo:hi]))               # mayor anomalia baja = candidato
    xs = x[j]

    fig = plt.figure(figsize=(14, 9.5))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.35, 0.06, 1, 1], height_ratios=[1, 1.9],
                          hspace=0.32, wspace=0.55)

    # (1) curvas normalizadas + zonas
    ax = fig.add_subplot(gs[0, :])
    fuerza = np.clip(-score / max(-score.min(), 1e-9), 0, 1) if score.min() < 0 else np.zeros(npts)
    for k in range(npts):
        if fuerza[k] > 0.05:
            ax.axvspan(x[k] - args.dx / 2, x[k] + args.dx / 2, color="cyan", alpha=0.15 + 0.5 * fuerza[k], lw=0)
    for i in range(nf):
        ax.plot(x, Vn[i], color="gray", lw=0.7)
    ax.axvline(xs, color="k", lw=2, alpha=0.8)
    ax.set_xlim(x.min() - args.dx / 2, x.max() + args.dx / 2)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Distancia x (m)")
    ax.set_ylabel("Respuesta EPD normalizada")
    ax.set_title("Curvas de frecuencia - zonas de agua subterránea (normalizado)", fontweight="bold")
    ax.grid(ls=":", alpha=0.5)

    # (2) seccion 2D con PowerNorm
    ax2 = fig.add_subplot(gs[1, 0])
    hh = np.concatenate([[0.0], h])
    Z = np.vstack([np.where(valido[:1], Vn[:1], 0.0), np.where(valido, Vn, 0.0)])
    pc = ax2.pcolormesh(x, hh, Z, shading="gouraud", cmap="jet",
                        norm=PowerNorm(gamma=args.gamma, vmin=0, vmax=1))
    ax2.axvline(xs, color="k", lw=2)
    mascara = np.nan_to_num(A, nan=0.0) < -umb_a
    lab, nl = ndimage.label(mascara)
    for k in range(1, nl + 1):
        ii, jj = np.where(lab == k)
        if len(ii) < args.min_celdas:
            continue
        y0b = hh[ii.min()] if ii.min() > 0 else h[0] - (h[1] - h[0]) / 2
        y1b = h[ii.max()] + (h[1] - h[0]) / 2
        ax2.add_patch(plt.Rectangle((x[jj.min()] - args.dx / 2, y0b), (x[jj.max()] - x[jj.min()]) + args.dx,
                                    y1b - y0b, fill=False, ec="blue", ls="--", lw=1.6))
    ax2.set_ylim(args.prof, 0)
    ax2.set_xlim(x.min(), x.max())
    ax2.set_xlabel("Distancia x (m)")
    ax2.set_ylabel("Profundidad (m)")
    ax2.set_title(f"Sección 2D EPD (γ={args.gamma:.2f})", fontweight="bold")
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label("EPD normalizada (PowerNorm)")

    # (3) "inversion" en la estacion elegida
    obs = valido[:, j]
    hv, vv = h[obs], V[obs, j]
    ax3 = fig.add_subplot(gs[1, 2])
    if obs.sum() >= 4:
        D = np.diff(np.eye(len(vv)), 2, axis=0)
        sc = vv.max()
        m = np.linalg.solve(np.eye(len(vv)) + args.lam * D.T @ D, vv / sc) * sc    # ajuste regularizado
        suav = ndimage.gaussian_filter1d(vv, args.factor_suavizado, mode="nearest")
        rms = float(np.sqrt(np.mean((m - vv) ** 2)))
        ax3.plot(vv, hv, ".", color="#6a5acd", ms=5, alpha=0.7, label="Observado")
        ax3.plot(suav, hv, color="#17becf", lw=1.3, label=f"Suavizado (factor {args.factor_suavizado:g})")
        ax3.plot(m, hv, color="red", lw=3, label=f"Calculado\nRMS={rms:.2f} ({100 * rms / sc:.1f}%)")
        ax3.legend(loc="lower right", fontsize=7)
        # (4) modelo de resistividad (mapeo empirico)
        ax4 = fig.add_subplot(gs[1, 3])
        t = np.clip((m - m.min()) / max(m.max() - m.min(), 1e-12), 0, 1)
        rho = args.rho_min * (args.rho_max / args.rho_min) ** t
        paso = h[1] - h[0]
        bordes = np.concatenate([hv - paso / 2, [hv[-1] + paso / 2]])
        ax4.step(np.concatenate([rho, rho[-1:]]), bordes, where="post", color="k", lw=1.8)
        ext = [0] + [i for i in range(1, len(rho) - 1) if (rho[i] - rho[i - 1]) * (rho[i + 1] - rho[i]) < 0] + [len(rho) - 1]
        for i in sorted(set(ext)):
            ax4.annotate(f"{rho[i]:.3g}", (rho[i], hv[i]), fontsize=6.5, fontweight="bold",
                         bbox=dict(boxstyle="square,pad=0.15", fc="#f5e663", ec="k", lw=0.6),
                         xytext=(3, 0), textcoords="offset points")
        ax4.set_xscale("log")
        ax4.set_ylim(args.prof, 0)
        ax4.set_xlabel("Resistividad (ohm·m)")
        ax4.set_ylabel("Profundidad (m)")
        ax4.set_title("Modelo de resistividad\n(aprox. empírica, ver README)", fontweight="bold", fontsize=10)
        ax4.grid(ls=":", alpha=0.5, which="both")
    ax3.set_ylim(args.prof, 0)
    ax3.set_xlim(left=0)
    ax3.set_xlabel("EPD (mV)")
    ax3.set_title(f"Inversión de EPD\nen x = {xs:g} m (N = {n[j]:g})", fontweight="bold", fontsize=10)
    ax3.grid(ls=":", alpha=0.5)
    base = os.path.splitext(os.path.basename(args.archivo))[0]
    args.nombre = "resumen_" + base if args.nombre in (None, base) else args.nombre
    os.makedirs(args.salida, exist_ok=True)
    ruta = os.path.join(args.salida, f"{args.nombre}.png")
    fig.savefig(ruta, dpi=args.dpi, bbox_inches="tight")
    print("Figura guardada en", ruta, f"(estacion x = {xs:g} m, N = {n[j]:g})")


def guardar(fig, args):
    os.makedirs(args.salida, exist_ok=True)
    ruta = os.path.join(args.salida, f"{args.nombre}.png")
    fig.tight_layout()
    fig.savefig(ruta, dpi=args.dpi)
    print("Figura guardada en", ruta)


# ----------------------------------------------------------------------------
# Datos sinteticos de ejemplo (solo para probar el programa, NO son del articulo)
# ----------------------------------------------------------------------------
def demo(args):
    os.makedirs(args.salida, exist_ok=True)
    y = np.concatenate([np.arange(-100, -20, 1.0), np.arange(-20, 20.5, 0.5), np.arange(21, 101, 1.0)])
    lgf = np.arange(1.0, 4.0001, 0.05)
    rows = []
    for l in lgf:
        base = 0.2 * 10 ** (1.2 * (l - 1.0) / 3 * 2.0) * 1.0  # crece con f
        base = 0.5 + 18.5 * ((l - 1.0) / 3.0) ** 2.2
        for yy in y:
            plato = 1 - 0.45 * np.exp(-(yy / 4.0) ** 2)            # fig. 2 (estrecho)
            esfera = 1 - 0.17 * np.exp(-(yy / 45.0) ** 2)           # fig. 4 (ancho)
            fase = 45 + 6 * (1 - np.exp(-(yy / 40.0) ** 2)) * (4 - l) / 3 - 1
            rows.append((yy, 10 ** l, base * plato, base * esfera,
                         base * plato * esfera, fase))
    d = pd.DataFrame(rows, columns=["y_m", "f_hz", "fig2", "fig4", "fig6", "fase_deg"])
    for k in ("fig2", "fig4"):
        d[["y_m", "f_hz"]].assign(ey_mv=d[k]).to_csv(
            os.path.join(args.salida, f"datos_{k}.csv"), index=False)
    d[["y_m", "f_hz", "fase_deg"]].assign(ey_mv=d["fig6"])[
        ["y_m", "f_hz", "ey_mv", "fase_deg"]].to_csv(
        os.path.join(args.salida, "datos_fig6.csv"), index=False)

    # Fig. 13: 40 frecuencias entre 12 y 5000 Hz, y = 15..29 m, anomalia en 23 m
    ys = np.arange(15, 30, 1.0)
    fs = np.logspace(np.log10(12), np.log10(5000), 40)
    rng = np.random.default_rng(0)
    r = []
    for fq in fs:
        for yy in ys:
            v = 3 + 30 * (np.log10(fq) - 1) / 2.7
            v *= 1 - 0.8 * np.exp(-((yy - 23.5) / 2.5) ** 2) * (0.4 + 0.6 * np.exp(-((np.log10(fq) - 1.8) / 0.4) ** 2))
            r.append((yy, fq, max(v * (1 + 0.03 * rng.standard_normal()), 0.3)))
    pd.DataFrame(r, columns=["y_m", "f_hz", "dv_mv"]).to_csv(
        os.path.join(args.salida, "datos_fig13_L8.csv"), index=False)
    print("Ejemplos escritos en", args.salida)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def comun(sp):
        sp.add_argument("--nombre", default=None, help="nombre del PNG de salida (sin extension)")
        sp.add_argument("--salida", default="figuras", help="carpeta de salida")
        sp.add_argument("--dpi", type=int, default=200)

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
    b.add_argument("--prof", type=float, default=None, choices=[100, 150, 300],
                   help="rango de profundidad configurado en el equipo (m); por defecto se lee del nombre del archivo")
    b.add_argument("--escala-h", choices=["lineal", "ec12"], default="lineal",
                   help="(formato equipo) profundidad: 'lineal' con el canal, como el Profile del equipo (def.), "
                        "o 'ec12' = c*503*sqrt(rho/f) del articulo")
    b.add_argument("--hmax", type=float, default=None, help="profundidad maxima mostrada (m)")
    b.add_argument("--zk", type=float, default=None, help="posicion (m) del sondeo ZK a marcar")
    b.add_argument("--freqs", default=None,
                   help="(formato equipo) frecuencias en Hz separadas por coma, o archivo con ellas; "
                        "una por columna freqNN, en el mismo orden")
    b.add_argument("--fmin", type=float, default=12.0, help="(formato equipo) fmin si no hay --freqs")
    b.add_argument("--fmax", type=float, default=5000.0, help="(formato equipo) fmax si no hay --freqs")
    b.add_argument("--dx", type=float, default=1.5, help="(formato equipo) metros entre puntos N consecutivos (def. 1.5)")
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
        sp.add_argument("--prof", type=float, default=None, choices=[100, 150, 300],
                        help="rango de profundidad del equipo (m); por defecto se lee del nombre del archivo")
        sp.add_argument("--dx", type=float, default=1.5, help="metros entre puntos N (def. 1.5)")
        sp.add_argument("--y0", type=float, default=0.0, help="posicion y (m) del menor N")
        sp.add_argument("--linea", type=int, default=None, help="filtrar por registro L del equipo")
        comun(sp)

    v = sub.add_parser("vertical", help="Fig. 3 (Gomo y Ngobe): perfil vertical de EPD vs profundidad")
    comun_eq(v)
    v.add_argument("--punto", default=None, help="N de la estacion (o varios separados por coma); def. el primero")
    v.add_argument("--litologia", default=None, help="CSV con columnas tope_m, base_m, nombre")
    v.add_argument("--agua", default=None, help="profundidades (m) de venas de agua, separadas por coma")
    v.set_defaults(fn=figura_vertical)

    hz = sub.add_parser("horizontal", help="Fig. 4 (Gomo y Ngobe): perfil horizontal de EPD (log), una curva por profundidad")
    comun_eq(hz)
    hz.add_argument("--profundidades", default=None, help="profundidades (m) a graficar, separadas por coma")
    hz.add_argument("--cada", type=int, default=2, help="si no hay --profundidades: 1 de cada N canales (def. 2)")
    hz.add_argument("--zona", default=None, help="recuadro a marcar: x1,x2[,ymin,ymax] (p. ej. zona de ruido)")
    hz.set_defaults(fn=figura_horizontal)

    rs = sub.add_parser("resumen", help="Figura resumen: curvas, seccion 2D, inversion en una estacion y modelo de resistividad")
    comun_eq(rs)
    rs.add_argument("--punto", type=float, default=None, help="N de la estacion a invertir (def.: la de mayor anomalia baja)")
    rs.add_argument("--x", type=float, default=None, help="posicion x (m) de la estacion (alternativa a --punto)")
    rs.add_argument("--gamma", type=float, default=0.30, help="exponente de PowerNorm de la seccion 2D (def. 0.30)")
    rs.add_argument("--umbral", type=float, default=0.1, help="lecturas <= umbral (mV) = ruido/canal muerto (def. 0.1)")
    rs.add_argument("--umbral-anomalia", dest="umbral_anomalia", type=float, default=0.15,
                    help="caida (log10) bajo la mediana lateral para marcar una zona baja (def. 0.15 = -30%%)")
    rs.add_argument("--min-celdas", dest="min_celdas", type=int, default=4, help="celdas minimas de un recuadro (def. 4)")
    rs.add_argument("--lam", type=float, default=3.0, help="regularizacion del ajuste (def. 3)")
    rs.add_argument("--factor-suavizado", dest="factor_suavizado", type=float, default=3.0, help="sigma del suavizado (canales)")
    rs.add_argument("--rho-min", dest="rho_min", type=float, default=5.0, help="resistividad asignada a la EPD minima (ohm.m)")
    rs.add_argument("--rho-max", dest="rho_max", type=float, default=500.0, help="resistividad asignada a la EPD maxima (ohm.m)")
    rs.set_defaults(fn=figura_resumen)

    c = sub.add_parser("demo", help="genera CSV sinteticos de ejemplo")
    c.add_argument("--salida", default="ejemplos")
    c.set_defaults(fn=demo)

    args = p.parse_args()
    if args.cmd != "demo" and args.nombre is None:
        args.nombre = os.path.splitext(os.path.basename(args.archivo))[0]
    args.fn(args)


if __name__ == "__main__":
    main()
