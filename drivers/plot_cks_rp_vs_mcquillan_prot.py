import os
import matplotlib.pyplot as plt, pandas as pd, numpy as np
from astropy import units as u, constants as c

from gilly.cks_age_exploration import (
    _get_cks_data, _apply_cks_VII_filters, _get_McQuillan13_data
)

from aesthetic.plot import format_ax, set_style, savefig

def plot_cks_rp_vs_mcquillan_prot():

    # note: ideally, would do this in a more simple way... _get_cks_data is
    # clearly deprecated.

    df = _get_cks_data()
    fp18_df = df[_apply_cks_VII_filters(df)]
    m13_df = _get_McQuillan13_data()

    fp18_df['II_id_kic'] = fp18_df['II_id_kic'].astype(str)
    m13_df['KIC'] = m13_df['KIC'].astype(str)

    mdf = fp18_df.merge(m13_df, how='inner', left_on='II_id_kic', right_on='KIC')

    #
    # make plots
    #
    set_style()

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf.Prot, s=2, c='k')
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel('M13 P$_\mathrm{rot}$ [day]')
        ax.set_xscale(scale)
        ax.set_yscale(scale)
        outpath = f'../results/cks_rotation_period/cks-VII_rp_vs_mcquillan_prot_{scale}.png'
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf.Prot, s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('M13 P$_\mathrm{rot}$ [day]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    outpath = '../results/cks_rotation_period/koi_period_vs_mcquillan_prot.png'
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf.VIIp_Rp, s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    format_ax(ax)
    outpath = '../results/cks_rotation_period/koi_period_vs_giso_prad_M13_match.png'
    savefig(f, outpath, writepdf=0)


if __name__ == "__main__":
    plot_cks_rp_vs_mcquillan_prot()
