#!/usr/bin/env python3
"""
Plota as variaveis de output_mp_offline.nc
(mp_offline_driver: MPAS Thompson + MYNN, forcado por DEPHY-SCM GOAMAZON/SINGLE).

Uso:
    python plot_mp_offline.py output_mp_offline.nc
    python plot_mp_offline.py output_mp_offline.nc -o figuras --ztop 18 --show

Gera PNGs em <outdir>:
    01_termodinamica.png   theta, qv, pressao
    02_hidrometeoros.png   qc, qr, qi, qs, qg
    03_pbl_mynn.png        tke, mixing_length, exch_h, exch_m
    04_subgrade_edmf.png   cldfra_bl, qc_bl, qi_bl, mf_a, mf_w
    05_superficie.png      precipitacao, neve, pblh
    06_relampagos.png      flash_rate (total, CG, IC)

Sem --show, usa o backend Agg (funciona em terminal/cluster sem display).
"""
import argparse
import os
import sys

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")  # sem display

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Configuracao dos paineis tempo-altura:
#   nome -> (titulo, fator, unidade, cmap, modo)
#   modo: "seq" (sequencial), "anom" (anomalia vs t=0, divergente),
#         "pos" (sequencial mascarando ~zero), "div" (divergente centrado em 0)
# ---------------------------------------------------------------------------
G = 1e3  # kg/kg -> g/kg
U_GKG = "g/kg"

GRUPOS_2D = {
    "01_termodinamica": [
        ("theta", "theta - theta(t=0)", 1.0, "K", "RdBu_r", "anom"),
        ("qv", "Razao de mistura de vapor (qv)", G, U_GKG, "YlGnBu", "seq"),
        ("pressure", "Pressao", 1e-2, "hPa", "viridis", "seq"),
    ],
    "02_hidrometeoros": [
        ("qc", "Agua de nuvem (qc)", G, U_GKG, "Blues", "pos"),
        ("qr", "Chuva (qr)", G, U_GKG, "GnBu", "pos"),
        ("qi", "Gelo de nuvem (qi)", G, U_GKG, "Purples", "pos"),
        ("qs", "Neve (qs)", G, U_GKG, "PuBu", "pos"),
        ("qg", "Graupel (qg)", G, U_GKG, "Oranges", "pos"),
    ],
    "03_pbl_mynn": [
        ("tke", "TKE (MYNN)", 1.0, "m2/s2", "magma_r", "pos"),
        ("mixing_length", "Comprimento de mistura (el_pbl)", 1.0, "m", "cividis", "seq"),
        ("exch_h", "Difusividade de calor (exch_h)", 1.0, "m2/s", "inferno_r", "pos"),
        ("exch_m", "Difusividade de momento (exch_m)", 1.0, "m2/s", "inferno_r", "pos"),
    ],
    "04_subgrade_edmf": [
        ("cldfra_bl", "Fracao de nuvem subgrade (cldfra_bl)", 1.0, "0-1", "Greys", "pos"),
        ("qc_bl", "qc subgrade (qc_bl)", G, U_GKG, "Blues", "pos"),
        ("qi_bl", "qi subgrade (qi_bl)", G, U_GKG, "Purples", "pos"),
        ("mf_a", "Fracao de area das correntes ascendentes (mf_a)", 1.0, "0-1", "YlOrRd", "pos"),
        ("mf_w", "w das correntes ascendentes (mf_w)", 1.0, "m/s", "RdBu_r", "div"),
    ],
}


def abrir(path):
    ds = xr.open_dataset(path, decode_times=False)  # "seconds since t0" nao e CF
    t_h = ds["time"].values / 3600.0
    z_km = ds["lev"].values / 1000.0
    return ds, t_h, z_km


def plot_tempo_altura(ax, t_h, z_km, campo, cmap, modo, ztop):
    """campo: array (time, lev) ja convertido de unidade."""
    campo = np.array(campo, dtype=float)
    if modo == "anom":
        campo = campo - campo[0:1, :]
    if modo == "pos":
        campo = np.ma.masked_less_equal(campo, 1e-9 * max(np.nanmax(campo), 1e-30))

    vmin = vmax = None
    if modo in ("anom", "div"):
        m = np.nanmax(np.abs(campo)) or 1.0
        vmin, vmax = -m, m

    pm = ax.pcolormesh(t_h, z_km, campo.T, cmap=cmap, shading="nearest",
                       vmin=vmin, vmax=vmax, rasterized=True)
    ax.set_ylim(0, ztop if ztop else z_km.max())
    ax.set_ylabel("Altura (km)")
    return pm


def fig_2d(ds, t_h, z_km, nome, paineis, ztop, outdir, show):
    paineis = [p for p in paineis if p[0] in ds]
    if not paineis:
        return
    n = len(paineis)
    fig, axs = plt.subplots(n, 1, figsize=(11, 2.8 * n + 0.6), sharex=True,
                            constrained_layout=True, squeeze=False)
    for ax, (var, titulo, fator, unid, cmap, modo) in zip(axs[:, 0], paineis):
        pm = plot_tempo_altura(ax, t_h, z_km, ds[var].values * fator, cmap, modo, ztop)
        ax.set_title(titulo, loc="left", fontsize=10)
        fig.colorbar(pm, ax=ax, pad=0.01, label=unid)
    axs[-1, 0].set_xlabel("Tempo (h desde t0)")
    salvar(fig, outdir, nome, show)


def largura_barra(t_h):
    return np.median(np.diff(t_h)) if len(t_h) > 1 else 0.01


def fig_superficie(ds, t_h, outdir, show):
    fig, axs = plt.subplots(3, 1, figsize=(11, 8), sharex=True, constrained_layout=True)
    w = largura_barra(t_h)

    # Precipitacao
    ax = axs[0]
    if "rainncv" in ds:
        ax.bar(t_h, ds["rainncv"].values, width=w, color="tab:blue", alpha=0.6,
               label="rainncv (incremental)")
    ax.set_ylabel("Incremental (mm)")
    if "rainnc" in ds:
        ax2 = ax.twinx()
        ax2.plot(t_h, ds["rainnc"].values, color="k", lw=1.5, label="rainnc (acumulada)")
        ax2.set_ylabel("Acumulada (mm)")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)
    ax.set_title("Precipitacao de superficie", loc="left", fontsize=10)

    # Neve (equivalente em agua liquida)
    ax = axs[1]
    if "snowncv" in ds:
        ax.bar(t_h, ds["snowncv"].values, width=w, color="tab:cyan", alpha=0.6,
               label="snowncv (incremental)")
    ax.set_ylabel("Incremental (mm LWE)")
    if "snownc" in ds:
        ax2 = ax.twinx()
        ax2.plot(t_h, ds["snownc"].values, color="k", lw=1.5, label="snownc (acumulada)")
        ax2.set_ylabel("Acumulada (mm LWE)")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)
    ax.set_title("Neve + gelo de nuvem (equivalente em agua liquida, NAO e profundidade)",
                 loc="left", fontsize=10)

    # PBLH
    ax = axs[2]
    if "pblh" in ds:
        ax.plot(t_h, ds["pblh"].values, color="tab:green", lw=1.5)
    ax.set_ylabel("pblh (m)")
    ax.set_xlabel("Tempo (h desde t0)")
    ax.set_title("Altura da camada limite (MYNN)", loc="left", fontsize=10)

    for a in axs:
        a.grid(alpha=0.3)
    salvar(fig, outdir, "05_superficie", show)


def fig_relampagos(ds, t_h, outdir, show):
    vars_ = [v for v in ("flash_rate", "flash_rate_cg", "flash_rate_ic") if v in ds]
    if not vars_:
        return
    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    estilos = {"flash_rate": ("Total", "k", 1.8),
               "flash_rate_cg": ("Nuvem-solo (CG)", "tab:red", 1.2),
               "flash_rate_ic": ("Intra-nuvem (IC)", "tab:blue", 1.2)}
    for v in vars_:
        lab, cor, lw = estilos[v]
        ax.plot(t_h, ds[v].values, label=lab, color=cor, lw=lw)
    ax.set_xlabel("Tempo (h desde t0)")
    ax.set_ylabel("flashes/min")
    ax.set_title("Taxa de relampagos diagnosticada", loc="left", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend()
    salvar(fig, outdir, "06_relampagos", show)


def salvar(fig, outdir, nome, show):
    path = os.path.join(outdir, nome + ".png")
    fig.savefig(path, dpi=150)
    print("salvo: " + path)
    if not show:
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arquivo", help="NetCDF de saida (output_mp_offline.nc)")
    ap.add_argument("-o", "--outdir", default="figuras_mp_offline")
    ap.add_argument("--ztop", type=float, default=None,
                    help="altura maxima do eixo vertical em km (padrao: topo do modelo)")
    ap.add_argument("--show", action="store_true", help="mostrar as figuras na tela")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ds, t_h, z_km = abrir(args.arquivo)
    print("%d passos, %d niveis, %.2f a %.2f h, topo ~ %.1f km" % (
        ds.sizes["time"], ds.sizes["lev"], t_h.min(), t_h.max(), z_km.max()))

    for nome, paineis in GRUPOS_2D.items():
        fig_2d(ds, t_h, z_km, nome, paineis, args.ztop, args.outdir, args.show)
    fig_superficie(ds, t_h, args.outdir, args.show)
    fig_relampagos(ds, t_h, args.outdir, args.show)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
