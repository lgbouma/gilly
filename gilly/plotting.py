import os
from os.path import join
import matplotlib.pyplot as plt, pandas as pd, numpy as np, matplotlib as mpl
from astropy import units as u, constants as c

from gilly.helpers import get_merged_M13_CKS, get_M13_CKS_gyro
from gilly.paths import RESULTSDIR

from aesthetic.plot import set_style, savefig

def plot_cks_rp_vs_m13_gyroage():

    mdf, fp18_df = get_M13_CKS_gyro()

    outdir = join(RESULTSDIR, 'gyroage')
    set_style()

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf.gyroage_yr, s=2, c='k', zorder=3)
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel('Gyro-Age [M13 P$_\mathrm{rot}$ + MH08; yr]')
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
    ax.set_ylabel('Gyro-Age [M13 P$_\mathrm{rot}$ + MH08; yr]')
    ax.set_xscale('log'); ax.set_yscale('log')
    ylim = ax.get_ylim()
    ymax = 1.5e10
    ax.set_ylim([ylim[0], ymax])
    outpath = join(outdir, 'cks-VII_period_vs_gyroage.png')
    savefig(f, outpath, writepdf=0)


    plt.close('all')
    f, ax = plt.subplots(figsize=(4.3,3))
    sel = ~pd.isnull(mdf.Prot)
    N_withProt = len(mdf[sel])
    N_noProt = len(fp18_df) - N_withProt

    l = 'Yes M13 P$_\mathrm{rot}$ '+f'({N_withProt})'
    norm = mpl.colors.Normalize(vmin=8.5, vmax=10)
    cax = ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=7,
                     c=np.log10(mdf[sel].gyroage_yr),
                     cmap='plasma_r', zorder=3, label=l, norm=norm,
                     edgecolors='k', linewidths=0.2)

    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log'); ax.set_yscale('log')

    cbar = f.colorbar(cax, extend='both')
    cbar.set_label('log$_{10}$Gyro-Age [M13 P$_\mathrm{rot}$ + MH08; yr]',
                   fontsize='small')
    cbar.minorticks_off()

    outpath = join(outdir, 'cks-VII_period_vs_cks-VII_prad_M13_match_gyroagecolors.png')
    savefig(f, outpath, writepdf=0)




def plot_cks_rp_vs_mcquillan_prot():

    mdf, fp18_df = get_merged_M13_CKS()

    set_style()

    outdir = join(RESULTSDIR, 'cks_rotation_period')

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf.Prot, s=2, c='k', zorder=3)
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel('M13 P$_\mathrm{rot}$ [day]')
        ax.set_xscale(scale); ax.set_yscale(scale)
        outpath = join(outdir, f'cks-VII_rp_vs_mcquillan_prot_{scale}.png')
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf.Prot, s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('M13 P$_\mathrm{rot}$ [day]')
    ax.set_xscale('log'); ax.set_yscale('log')
    outpath = join(outdir, 'cks-VII_period_vs_mcquillan_prot.png')
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
    outpath = join(outdir, 'cks-VII_period_vs_cks-VII_prad_M13_match.png')
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4.3,3))
    sel = ~pd.isnull(mdf.Prot)
    N_withProt = len(mdf[sel])
    N_noProt = len(fp18_df) - N_withProt

    l = 'Yes M13 P$_\mathrm{rot}$ '+f'({N_withProt})'
    norm = mpl.colors.Normalize(vmin=10, vmax=20)
    cax = ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=7, c=mdf[sel].Prot,
                     cmap='plasma_r', zorder=3, label=l, norm=norm,
                     edgecolors='k', linewidths=0.2)

    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log'); ax.set_yscale('log')

    cbar = f.colorbar(cax, extend='both')
    cbar.set_label('M13 P$_\mathrm{rot}$ [day]')
    cbar.minorticks_off()

    outpath = join(outdir, 'cks-VII_period_vs_cks-VII_prad_M13_match_M13-Protcolors.png')
    savefig(f, outpath, writepdf=0)
