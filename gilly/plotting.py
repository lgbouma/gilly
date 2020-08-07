"""
contents:

    plot_spec_vs_phot_Prot
    plot_spec_Prot_hist

    plot_cdf_rp_agecut_KS_test
    plot_hist_rp_agecut
    plot_stsnr_vs_gyroage
    plot_cks_rp_vs_gyroage
    plot_cks_rp_vs_prot
    plot_spec_vs_phot_Prot
"""
import os
from os.path import join
import matplotlib.pyplot as plt, pandas as pd, numpy as np, matplotlib as mpl
from astropy import units as u, constants as c

from gilly.helpers import (
    get_merged_rot_CKS, get_merged_gyroage_CKS, _add_Christiansen12_CDPP,
    _get_spec_vsini_phot_Prot_overlap
)
from gilly.paths import DATADIR, RESULTSDIR

from aesthetic.plot import set_style, savefig
from numpy import array as arr
from matplotlib.lines import Line2D

def plot_spec_Prot_hist(Prot_source='M15', spec_source='cks'):
    """
    Prot_source (str): M13 or M15
    spec_source (str): always 'cks', because that's the only vsini source.
        [nb. this is the vsini dataset from J. Winn, received 2 Aug 2020]
    """

    if spec_source != 'cks':
        raise NotImplementedError

    mdf, cdf = _get_spec_vsini_phot_Prot_overlap(Prot_source=Prot_source)
    sel = cdf['id_kic'].isin(arr(mdf['id_kic']))
    scdf = cdf[~sel]

    print(f'{len(cdf)} with vsini')
    print(f'{len(mdf[(~pd.isnull(mdf.spec_Prot)) & (~pd.isnull(mdf.Prot))])} with vsini & M15 prot')
    print(f'{len(scdf)} with vsini & and no M15 prot')

    outdir = join(RESULTSDIR, 'spec_rot')

    set_style()

    f, ax = plt.subplots(figsize=(4,3))

    bins = np.arange(0, 52.5, 2.5)

    N_m = len(mdf)
    N_n = len(scdf)

    ax.hist(scdf.spec_Prot, bins=bins, cumulative=False, fill=False, density=False,
            weights=np.ones(N_n)/N_n,
            histtype='step', label=f'No M15 match ({len(scdf)})')

    ax.hist(mdf.spec_Prot, bins=bins, cumulative=False, fill=False, density=False,
            weights=np.ones(N_m)/N_m,
            histtype='step', label=f'Has M15 match ({len(mdf)})')

    # Create new legend handles but use the colors from the existing ones
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
    ax.legend(handles=new_handles, labels=labels, fontsize='small')

    txt = (
        r'$\mu_{T_{\mathrm{eff}}} \pm \sigma_{T_{\mathrm{eff}}}$ '+' for P$_{\mathrm{rot}}^{\mathrm{spec}} < 25$ d:\n'
        'No M15: '+f'{int(scdf[scdf.spec_Prot<25].cks_steff.mean())} $\pm$ {int(scdf[scdf.spec_Prot<25].cks_steff.std())}\n'
        'Has M15: '+f'{int(mdf[mdf.spec_Prot<25].VIIs_Teff.mean())} $\pm$ {int(mdf[mdf.spec_Prot<25].VIIs_Teff.std())}\n'
    )
    ax.text(0.95, 0.50, txt,
            transform=ax.transAxes, ha='right', va='center',
            fontsize='small')

    ax.set_ylabel('Fraction per bin')
    ax.set_xlabel('CKS-VII P$_{\mathrm{rot}}^{\mathrm{spec}}$ [day]')

    outpath = join(
        outdir,
        f'spec_Prot_hist_prot{Prot_source}.png'
    )
    savefig(f, outpath, writepdf=0)


def plot_spec_vs_phot_Prot(Prot_source='M15', spec_source='cks'):
    """
    Prot_source (str): M13 or M15
    spec_source (str): always 'cks', because that's the only vsini source.
        [nb. this is the vsini dataset from J. Winn, received 2 Aug 2020]
    """

    if spec_source != 'cks':
        raise NotImplementedError

    mdf, _  = _get_spec_vsini_phot_Prot_overlap(Prot_source=Prot_source)

    outdir = join(RESULTSDIR, 'spec_rot')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    set_style()

    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.Prot, mdf.spec_Prot, s=2, c='k', zorder=3)
    _, xmax = ax.get_xlim()
    span = np.linspace(2, xmax, 100)
    ax.plot(span, span, zorder=-1, lw=1, ls='--', c='lightgray',
            label='$i=90^\circ$')
    ax.plot(span, span/np.sin(np.deg2rad(45)), zorder=-1, lw=1, ls=':', c='lightgray',
            label='$i=45^\circ$')
    ax.legend(fontsize='small')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(f'{Prot_source} phot '+'P$_\mathrm{rot}$ [day]')
    ax.set_ylabel(f'CKS-VII spec '+'P$_\mathrm{rot}$ [day]')
    outpath = join(outdir, f'cks-VII_specrot_vs_{Prot_source}_prot.png')
    savefig(f, outpath, writepdf=0)

    sel = (mdf.Prot < 10)
    sdf = mdf[sel]
    sdf['specProt/photProt'] = sdf.spec_Prot/sdf.Prot
    sdf = sdf.sort_values(by='specProt/photProt', ascending=False)
    selcols = ['VIIp_KOI', 'Prot', 'spec_Prot']
    print(sdf[selcols][:20].to_string(index=False))


def plot_cdf_rp_agecut_KS_test(agecut=1e9, Prot_source='M15',
                               gyro_source='A19'):
    """
    Prot_source (str): M13 or M15
    gyro_source (str): MH08 or A19
    agecut (float): Age cut in years

    Plot CDF, do KS test, calculate p-values.
    """

    mdf = get_merged_gyroage_CKS(Prot_source=Prot_source,
                                 gyro_source=gyro_source)

    odf = mdf[(mdf.gyroage_yr > agecut) & (mdf.VIIp_Rp < 6)] # old
    ydf = mdf[(mdf.gyroage_yr <= agecut) & (mdf.VIIp_Rp < 6)] # young
    N_o = len(odf)
    N_y = len(ydf)

    #
    # make the plot!
    #
    outdir = join(
        RESULTSDIR, f'cdf_rp_agecut_KS_test'
    )
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    set_style()

    f, ax = plt.subplots(figsize=(4,3))

    bins = np.linspace(mdf.VIIp_Rp.min(), mdf.VIIp_Rp.max(), 1000)

    cnt, bin_edges = np.histogram(arr(odf.VIIp_Rp), bins=bins, normed=True)
    cdf = np.cumsum(cnt)
    ax.plot(bin_edges[:-1], cdf/cdf[-1],
            label=f'Age > {agecut/1e9:.1f} Gyr ({N_o})')

    cnt, bin_edges = np.histogram(arr(ydf.VIIp_Rp), bins=bins, normed=True)
    cdf = np.cumsum(cnt)
    ax.plot(bin_edges[:-1], cdf/cdf[-1],
            label=f'Age < {agecut/1e9:.1f} Gyr ({N_y})')

    ax.set_xscale('log')

    from scipy.stats import ks_2samp
    D, p_value = ks_2samp(arr(odf.VIIp_Rp), arr(ydf.VIIp_Rp))

    txt = f'D={D:.2e}, p={p_value:.2e}'

    ax.text(0.05, 0.95, txt,
            transform=ax.transAxes, ha='left', va='top',
            fontsize='small')

    ax.legend(fontsize='small', loc='lower right')

    ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')

    outpath = join(
        outdir,
        f'cdf_Rp_KS_cut{agecut/1e9:.1f}_prot{Prot_source}_gyro{gyro_source}.png'
    )
    savefig(f, outpath, writepdf=0)




def plot_hist_rp_agecut(agecut=1e9, Prot_source='M15', gyro_source='A19',
                        fix_ylim=None):
    """
    Prot_source (str): M13 or M15
    gyro_source (str): MH08 or A19
    agecut (float): Age cut in years

    Plot histogram. $$$?
    """

    mdf = get_merged_gyroage_CKS(Prot_source=Prot_source,
                                 gyro_source=gyro_source)

    odf = mdf[(mdf.gyroage_yr > agecut) & (mdf.VIIp_Rp < 6)] # old
    ydf = mdf[(mdf.gyroage_yr <= agecut) & (mdf.VIIp_Rp < 6)] # young
    N_o = len(odf)
    N_y = len(ydf)

    # Following Fulton et al 2017, we define a “super-Earth” as a
    # planet with a radius of 1–1.75 Re, and a “sub-Neptune” as
    # having a radius of 1.75–4.0 Re.

    # old: number of superNeptunes /  number of super Earths
    o_NsN = len(odf[(odf.VIIp_Rp > 1.75) & (odf.VIIp_Rp <= 4.00)])
    o_NsE = len(odf[(odf.VIIp_Rp > 1.00) & (odf.VIIp_Rp <= 1.75)])
    unc_o_NsN = np.sqrt(len(odf[(odf.VIIp_Rp > 1.75) & (odf.VIIp_Rp <= 4.00)]))
    unc_o_NsE = np.sqrt(len(odf[(odf.VIIp_Rp > 1.00) & (odf.VIIp_Rp <= 1.75)]))

    o_NsN_div_NsE = o_NsN / o_NsE
    unc_o_div = o_NsN_div_NsE * np.sqrt( (unc_o_NsN/o_NsN)**2 + (unc_o_NsE/o_NsE)**2 )

    # ditto for young
    y_NsN = len(ydf[(ydf.VIIp_Rp > 1.75) & (ydf.VIIp_Rp <= 4.00)])
    y_NsE = len(ydf[(ydf.VIIp_Rp > 1.00) & (ydf.VIIp_Rp <= 1.75)])
    unc_y_NsN = np.sqrt(len(ydf[(ydf.VIIp_Rp > 1.75) & (ydf.VIIp_Rp <= 4.00)]))
    unc_y_NsE = np.sqrt(len(ydf[(ydf.VIIp_Rp > 1.00) & (ydf.VIIp_Rp <= 1.75)]))

    y_NsN_div_NsE = y_NsN / y_NsE
    unc_y_div = y_NsN_div_NsE * np.sqrt( (unc_y_NsN/y_NsN)**2 + (unc_y_NsE/y_NsE)**2 )

    txt = (
        'N$_{\mathrm{subNeptune}}$/N$_{\mathrm{superEarth}}$:'
        +
        f'\nOld: {o_NsN_div_NsE:.2f} $\pm$ {unc_o_div:.2f}'
        +
        f'\nYoung: {y_NsN_div_NsE:.2f} $\pm$ {unc_y_div:.2f}'
    )

    #
    # make the plot!
    #
    outdir = join(
        RESULTSDIR, f'hist_rp_agecut'
    )
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    set_style()

    f, ax = plt.subplots(figsize=(4,3))

    bins = np.arange(0, 6.25, 0.25)

    ax.hist(odf.VIIp_Rp, bins=bins, cumulative=False, fill=False, density=False,
            weights=np.ones(N_o)/N_o,
            histtype='step', label=f'Age > {agecut/1e9:.1f} Gyr ({N_o})')

    ax.hist(ydf.VIIp_Rp, bins=bins, cumulative=False, fill=False, density=False,
            weights=np.ones(N_y)/N_y,
            histtype='step', label=f'Age < {agecut/1e9:.1f} Gyr ({N_y})')

    if fix_ylim:
        ax.set_ylim([0,0.28])

    ymin, ymax = ax.get_ylim()
    ax.vlines(
        1.74, ymin, ymax, colors='darkgray', alpha=1,
        linestyles='--', zorder=-2, linewidths=0.5
    )
    ax.set_ylim((ymin, ymax))

    # Create new legend handles but use the colors from the existing ones
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

    ax.legend(handles=new_handles, labels=labels, fontsize='small')

    ax.text(0.95, 0.50, txt,
            transform=ax.transAxes, ha='right', va='center',
            fontsize='small')

    ax.set_ylabel('Fraction per bin')
    ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')

    outpath = join(
        outdir,
        f'hist_Rp_cut{agecut/1e9:.1f}_prot{Prot_source}_gyro{gyro_source}.png'
    )
    savefig(f, outpath, writepdf=0)




def plot_stsnr_vs_gyroage(Prot_source='M15', gyro_source='A19', cdpp_id='rms3'):
    """
    Prot_source: M13 or M15
    gyro_source: MH08 or A19
    cdpp_id: rms3, rms6, or rms12

    Makes a Berger+2020 style single-transit SNR vs age plot.  The method is
    identical to Berger+2020. CDPP is taken from Christiansen 2012.
    Single-transit SNR is calculated by doing (Rp/R*)^2 * (1/CDPP). This was
    cited as being from Petigura+2018. [.....]
    """
    mdf = get_merged_gyroage_CKS(Prot_source=Prot_source,
                                 gyro_source=gyro_source)
    mdf = _add_Christiansen12_CDPP(mdf)

    Rp = 1*u.Rearth
    Rstar = arr(mdf.VIIs_R)*u.Rsun
    Rp_Rs = (Rp/Rstar).cgs.value
    mdf['stsnr'] = (
        Rp_Rs**2
        *
        (1/arr(1e-6*mdf[cdpp_id]))
    )

    scols = ['gyroage_yr', 'stsnr']
    sdf = mdf[scols]

    outdir = join(
        RESULTSDIR, f'singletransitSNR_{Prot_source}rot_{gyro_source}age'
    )
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    set_style()

    f, ax = plt.subplots(figsize=(4,3))

    ax.scatter(sdf.gyroage_yr, sdf.stsnr, s=2, c='k', zorder=3)

    agelabel = (f'Gyro-Age [{Prot_source} '+'P$_\mathrm{rot}$ +'
                +f'{gyro_source}; yr]')
    ax.set_xlabel(agelabel, fontsize='small')

    ax.set_ylabel(f'Single-transit SNR [C12 CDPP{cdpp_id.replace("rms","")}]',
                  fontsize='small')

    ax.set_xscale('log'); ax.set_yscale('log')

    outpath = join(outdir,
                   f'singletransitsnr_{cdpp_id}_rot{Prot_source}_gyro{gyro_source}.png')
    savefig(f, outpath, writepdf=0)



def plot_cks_rp_vs_gyroage(Prot_source='M15', gyro_source='A19'):
    """
    Prot_source: M13 or M15
    gyro_source: MH08 or A19

    Make scatter plots of {plantet radius, orbital period} against
    gyrochronology ages in a few different projections.
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
    Prot_key = 'Prot' if Prot_source in ['M13','M15'] else 'spec_Prot'
    sel = ~pd.isnull(mdf[Prot_key])

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

    Make scatter plots of {plantet radius, orbital period} against
    rotation periods in a few different projections.
    """

    mdf, fp18_df = get_merged_rot_CKS(Prot_source=Prot_source)
    Prot_key = 'Prot' if Prot_source in ['M13','M15'] else 'spec_Prot'

    set_style()

    outdir = join(RESULTSDIR, f'cks_{Prot_source}_rotation_period')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for scale in ['linear','log']:
        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(mdf.VIIp_Rp, mdf[Prot_key], s=2, c='k', zorder=3)
        ax.set_xlabel('CKS-VII R$_\mathrm{p}$ [R$_\oplus$]')
        ax.set_ylabel(f'{Prot_source} '+'P$_\mathrm{rot}$ [day]')
        ax.set_xscale(scale); ax.set_yscale(scale)
        outpath = join(outdir, f'cks-VII_rp_vs_{Prot_source}_prot_{scale}.png')
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    ax.scatter(mdf.VIIp_Per, mdf[Prot_key], s=2, c='k')
    ax.set_xlabel('CKS-VII P$_\mathrm{orb}$ [day]')
    ax.set_ylabel(f'{Prot_source} '+'P$_\mathrm{rot}$ [day]')
    ax.set_xscale('log'); ax.set_yscale('log')
    outpath = join(outdir, f'cks-VII_period_vs_{Prot_source}_prot.png')
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))
    sel = ~pd.isnull(mdf[Prot_key])
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
    sel = ~pd.isnull(mdf[Prot_key])
    N_withProt = len(mdf[sel])
    N_noProt = len(fp18_df) - N_withProt

    l = f'Yes {Prot_source}'+' P$_\mathrm{rot}$ '+f'({N_withProt})'
    norm = mpl.colors.Normalize(vmin=5, vmax=25)
    cax = ax.scatter(mdf[sel].VIIp_Per, mdf[sel].VIIp_Rp, s=7,
                     c=mdf[sel][Prot_key],
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
