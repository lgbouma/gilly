import os
from os.path import join
import matplotlib.pyplot as plt, pandas as pd, numpy as np, matplotlib as mpl
from astropy import units as u, constants as c

from gilly.helpers import (
    get_merged_rot_CKS, get_merged_gyroage_CKS
)
from gilly.paths import RESULTSDIR

from aesthetic.plot import set_style, savefig

def plot_cks_rp_vs_gyroage(Prot_source='M15', gyro_source='A19'):
    """
    Prot_source: M13 or M15
    gyro_source: MH08 or A19
    """
    mdf = get_merged_gyroage_CKS(Prot_source=Prot_source,
                                 gyro_source=gyro_source)

    outdir = join(RESULTSDIR, f'{Prot_source}rot_{gyro_source}gyroage')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    set_style()

    agelabel = (f'Gyro-Age [{Prot_source} '+'P$_\mathrm{rot}$ +'
                +f'{gyro_source}; yr]')

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf.gyroage_yr, s=2, c='k', zorder=3)
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel(agelabel)
        ax.set_xscale(scale); ax.set_yscale('log')
        ylim = ax.get_ylim()
        ymax = 1.5e10
        ax.set_ylim([ylim[0], ymax])
        print(f'WRN! {len(mdf[mdf.gyroage_yr > ymax])} older than max plotted gyroage')
        outpath = join(outdir, f'cks-VII_rp_vs_gyroage_{scale}.png')
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf.gyroage_yr, s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel(agelabel)
    ax.set_xscale('log'); ax.set_yscale('log')
    ylim = ax.get_ylim()
    ymax = 1.5e10
    ax.set_ylim([ylim[0], ymax])
    outpath = join(outdir, 'cks-VII_period_vs_gyroage.png')
    savefig(f, outpath, writepdf=0)


    plt.close('all')
    f, ax = plt.subplots(figsize=(4.3,3))
    sel = ~pd.isnull(mdf.Prot)

    norm = mpl.colors.Normalize(vmin=8.5, vmax=10)
    cax = ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=7,
                     c=np.log10(mdf[sel].gyroage_yr),
                     cmap='plasma_r', zorder=3, norm=norm,
                     edgecolors='k', linewidths=0.2)

    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log'); ax.set_yscale('log')

    cbar = f.colorbar(cax, extend='both')
    cbar.set_label('log$_{10}$'+agelabel, fontsize='small')
    cbar.minorticks_off()

    outpath = join(outdir, 'cks-VII_period_vs_cks-VII_prad_gyroagecolors.png')
    savefig(f, outpath, writepdf=0)


def plot_cks_rp_vs_prot(Prot_source='M15'):
    """
    Prot_source: M13 or M15
    """

    mdf, fp18_df = get_merged_rot_CKS(Prot_source=Prot_source)

    set_style()

    outdir = join(RESULTSDIR, f'cks_{Prot_source}_rotation_period')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf.Prot, s=2, c='k', zorder=3)
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel(f'{Prot_source} '+'P$_\mathrm{rot}$ [day]')
        ax.set_xscale(scale); ax.set_yscale(scale)
        outpath = join(outdir, f'cks-VII_rp_vs_{Prot_source}_prot_{scale}.png')
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf.Prot, s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel(f'{Prot_source} '+'P$_\mathrm{rot}$ [day]')
    ax.set_xscale('log'); ax.set_yscale('log')
    outpath = join(outdir, f'cks-VII_period_vs_{Prot_source}_prot.png')
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    sel = ~pd.isnull(mdf.Prot)
    N_withProt = len(mdf[sel])
    N_noProt = len(fp18_df) - N_withProt
    ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=2, c='k', zorder=3,
               label=f'Yes {Prot_source} '+'P$_\mathrm{rot}$ '+f'({N_withProt})')
    ax.scatter(fp18_df.VIIp_Per, fp18_df.VIIp_Rp, s=1, c='darkgray', zorder=2,
               label=f'No {Prot_source} '+'P$_\mathrm{rot}$ '+f'({N_noProt})')
    ax.legend(fontsize='small')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log'); ax.set_yscale('log')
    outpath = join(outdir, f'cks-VII_period_vs_cks-VII_prad_{Prot_source}_match.png')
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4.3,3))
    sel = ~pd.isnull(mdf.Prot)
    N_withProt = len(mdf[sel])
    N_noProt = len(fp18_df) - N_withProt

    l = f'Yes {Prot_source}'+' P$_\mathrm{rot}$ '+f'({N_withProt})'
    norm = mpl.colors.Normalize(vmin=5, vmax=25)
    cax = ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=7, c=mdf[sel].Prot,
                     cmap='plasma_r', zorder=3, label=l, norm=norm,
                     edgecolors='k', linewidths=0.2)

    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log'); ax.set_yscale('log')

    cbar = f.colorbar(cax, extend='both')
    cbar.set_label(f'{Prot_source}'+' P$_\mathrm{rot}$ [day]')
    cbar.minorticks_off()

    outpath = join(outdir, f'cks-VII_period_vs_cks-VII_prad_{Prot_source}-Protcolors.png')
    savefig(f, outpath, writepdf=0)
