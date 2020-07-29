import os
import matplotlib.pyplot as plt, pandas as pd, numpy as np
from astropy import units as u, constants as c

from gilly.helpers import get_merged_M13_CKS

from aesthetic.plot import format_ax, set_style, savefig

def plot_cks_rp_vs_mcquillan_prot():

    mdf, fp18_df = get_merged_M13_CKS()

    set_style()

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf.Prot, s=2, c='k', zorder=3)
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel('M13 P$_\mathrm{rot}$ [day]')
        ax.set_xscale(scale); ax.set_yscale(scale)
        outpath = f'../results/cks_rotation_period/cks-VII_rp_vs_mcquillan_prot_{scale}.png'
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf.Prot, s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('M13 P$_\mathrm{rot}$ [day]')
    ax.set_xscale('log'); ax.set_yscale('log')
    outpath = '../results/cks_rotation_period/cks-VII_period_vs_mcquillan_prot.png'
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    sel = ~pd.isnull(mdf.Prot)
    N_withProt = len(mdf[sel])
    N_noProt = len(fp18_df) - N_withProt
    ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=2, c='k', zorder=3,
               label='Yes M13 P$_\mathrm{rot}$ '+f'({N_withProt})')
    ax.scatter(fp18_df.VIIp_Per, fp18_df.VIIp_Rp, s=1, c='darkgray', zorder=2,
               label='No M13 P$_\mathrm{rot}$ '+f'({N_noProt})')
    ax.legend(fontsize='small')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log'); ax.set_yscale('log')
    format_ax(ax)
    outpath = '../results/cks_rotation_period/cks-VII_period_vs_cks-VII_prad_M13_match.png'
    savefig(f, outpath, writepdf=0)
