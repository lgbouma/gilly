import matplotlib.pyplot as plt, pandas as pd, numpy as np
import os

from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
from astropy.io.votable import parse

from numpy import array as nparr

from gilly.cks_age_exploration import (
    _get_cks_data, _apply_cks_IV_metallicity_study_filters
)

from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

from aesthetic.plot import format_ax, set_style, savefig

def main():

    set_style()

    df = _get_cks_data()
    sel = _apply_cks_IV_metallicity_study_filters(df)
    p17_df = df[sel]

    vot = parse('../data/McQuillan_2013_ApJ_775L.vot')
    m13 = vot.get_first_table()
    m13_df = m13.to_table().to_pandas()

    p17_df['id_kic'] = p17_df['id_kic'].astype(str)
    m13_df['KIC'] = m13_df['KIC'].astype(str)

    mdf = p17_df.merge(m13_df, how='inner', left_on='id_kic', right_on='KIC')

    for scale in ['linear','log']:
        f, ax = plt.subplots()
        ax.scatter(mdf.giso_prad, mdf.Prot, s=2, c='k')
        #FIXME which planet radius is shown? should be CKS VII
        ax.set_xlabel('CKS Planet Size [R$_\oplus$]')
        ax.set_ylabel('M13 Rotation Period [day]')
        ax.set_xscale(scale)
        ax.set_yscale(scale)
        format_ax(ax)
        outpath = f'../results/cks_rotation_period/cks_rp_vs_mcquillan_prot_{scale}.png'
        savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots()
    ax.scatter(mdf.koi_period, mdf.Prot, s=2, c='k')
    ax.set_xlabel('KOI Period [day]')
    ax.set_ylabel('M13 Rotation Period [day]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    format_ax(ax)
    outpath = '../results/cks_rotation_period/koi_period_vs_mcquillan_prot.png'
    savefig(f, outpath, writepdf=0)

    plt.close('all')
    f, ax = plt.subplots()
    ax.scatter(mdf.koi_period, mdf.giso_prad, s=2, c='k')
    ax.set_xlabel('KOI Period [day]')
    ax.set_ylabel('CKS Planet Size [R$_\oplus$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    format_ax(ax)
    outpath = '../results/cks_rotation_period/koi_period_vs_giso_prad_M13_match.png'
    savefig(f, outpath, writepdf=0)




if __name__ == "__main__":
    main()
