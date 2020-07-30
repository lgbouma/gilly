"""
Contents:

    get_M13_CKS_gyro
    get_merged_M13_CKS

    _get_McQuillan13_data
    _get_cks_data
    _apply_cks_VII_filters
"""
import matplotlib.pyplot as plt, pandas as pd, numpy as np
import os

from os.path import join
from numpy import array as arr

from astropy.table import Table
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u, astropy.constants as c
from astropy.io.votable import parse

from gilly.paths import DATADIR, RESULTSDIR

def get_merged_M13_CKS():
    """
    Return a DataFrame containing the Fulton & Petigura (2018) CKS stellar and
    planetary parameters, after the selection function has been applied. Merge
    it against McQuillan+13 rotation periods.
    """

    df = _get_cks_data()
    fp18_df = df[_apply_cks_VII_filters(df)]
    m13_df = _get_McQuillan13_data()

    fp18_df['II_id_kic'] = fp18_df['II_id_kic'].astype(str)
    m13_df['KIC'] = m13_df['KIC'].astype(str)

    mdf = fp18_df.merge(m13_df, how='inner', left_on='II_id_kic', right_on='KIC')

    return mdf, fp18_df

def get_M13_CKS_gyro():
    """
    Same as get_merged_M13_CKS, but then calculate the gyro ages using the
    Mamajek & Hillenbrand (2008) relation.

    This relation uses B-V instead of Teff.

    Two possible approaches:

        1. Interpolate Teff->B-V from the Mamajek & Pecaut table.

        2. Query exoplanet archive... hope that KIC has B and V? Not viable.
        Because nasa exoplanet archive KOI table has griz, JHKs, Kepler mags.
        (Ugh. Maybe KOI IDs --> Gaia IDs --> TIC IDs --> any available TIC B &
        V).

    This implements #1.
    """

    df = _get_cks_data()
    fp18_df = df[_apply_cks_VII_filters(df)]
    m13_df = _get_McQuillan13_data()

    fp18_df['II_id_kic'] = fp18_df['II_id_kic'].astype(str)
    m13_df['KIC'] = m13_df['KIC'].astype(str)

    mdf = fp18_df.merge(m13_df, how='inner', left_on='II_id_kic', right_on='KIC')

    from cdips.utils.mamajek import get_interp_BmV_from_Teff
    BmV = get_interp_BmV_from_Teff(arr(mdf.VIIs_Teff))

    mdf['BmV'] = BmV

    #
    # Mamajek & Hillenbrand 08 used
    #
    # Prot = a * [ (B-V)_0 - c ]^{b} * (t/Myr)^{n}
    #
    # and found some numbers. So
    #
    # (t/Myr) = [  (Prot/a) / [ (B-V)_0 - c ]^{b} ]^{1/n}
    #

    a = 0.407
    a_err = 0.021
    b = 0.325
    b_err = 0.024
    c = 0.495
    c_err = 0.010
    n = 0.566
    n_err = 0.008

    t_Myr = ( (mdf.Prot / a) / ( (mdf.BmV - c)**(b) ) )**(1/n)

    t_yr = t_Myr * 1e6

    mdf['gyroage_yr'] = t_yr

    return mdf, fp18_df



def _get_McQuillan13_data():

    vot = parse('../data/McQuillan_2013_ApJ_775L.vot')
    m13 = vot.get_first_table()
    m13_df = m13.to_table().to_pandas()

    return m13_df


def _get_cks_data():
    '''
    Returns:

        Dataframe with CKS VII planet and stellar data, supplemented by CKS II
        transit fit parameters (e.g., b), and Furlan+2017 dilution column.

    Description:

        There are 1305 CKS spectra of KOIs with 2025 planet candidates. The
        overlapping CKS stellar samples are: (1) magnitude-limited, Kp < 14.2,
        (2) multi-planet hosts (3) USP hosts (4) HZ hosts (5) other.

        The CDS page of Fulton & Petigura (2018) at
        http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/156/264 gives:
            Table 2: stellar properties (1189 rows)
            Table 3: planet properties (1901 rows)
            Table 4: planet detection statistics (907 rows)

        We couldn't just use _only_ these planet properties, because the full
        set of transit parameters, including for instance the impact parameter,
        was not published in VII, but was in II.
    '''

    hdul = fits.open(join(DATADIR, 'Fulton_2018_CKS_VII_cds_table.fits'))

    df_VII_planets = Table(hdul[2].data).to_pandas()
    renamedict = {c:'VIIp_'+c for c in list(df_VII_planets.columns)}
    df_VII_planets.rename(columns=renamedict, inplace=True)

    df_VII_stars = Table(hdul[1].data).to_pandas()
    renamedict = {c:'VIIs_'+c for c in list(df_VII_stars.columns)}
    df_VII_stars.rename(columns=renamedict, inplace=True)

    df_VII_planets['VIIp_id_starname'] = arr(
        df_VII_planets['VIIp_KOI'].str.slice(start=0,stop=6)
    ).astype(str)
    df_VII_stars['VIIs_KOI'] = arr(df_VII_stars['VIIs_KOI']).astype(str)

    df_m = pd.merge(df_VII_planets, df_VII_stars, how='left',
                    left_on='VIIp_id_starname', right_on='VIIs_KOI')

    df_II = pd.read_csv('../data/cks_physical_merged.csv')
    renamedict = {c:'II_'+c for c in list(df_II.columns)}
    df_II.rename(columns=renamedict, inplace=True)

    # pull non-changed stellar parameters from CKS II (needed for filters)
    subcols = ['kic_kepmag', 'id_starname', 'cks_fp', 'koi_impact',
               'koi_count', 'koi_dor', 'koi_model_snr', 'id_koicand',
               'id_kic', 'koi_steff', 'koi_slogg', 'koi_dor_err1',
               'koi_dor_err2']
    for ix, s in enumerate(subcols):
        subcols[ix] = 'II_'+s

    df = pd.merge(df_m, df_II[subcols], how='left', left_on='VIIp_KOI',
                  right_on='II_id_koicand')

    assert np.array_equal(arr(df['VIIp_id_starname']),
                          arr(df['II_id_starname']))
    df['id_starname'] = df['II_id_starname']
    df = df.drop('VIIp_id_starname', 1)
    df = df.drop('II_id_starname', 1)

    return df


def _apply_cks_VII_filters(df):
    """
    Apply filters from Section 4.2 of CKS VII (Fulton & Petigura 2018)
    """

    # 1. Magnitude-limited CKS subsample, Kp<14.2
    sel = np.isfinite(arr(df['II_kic_kepmag']))
    sel &= arr(df['II_kic_kepmag']) < 14.2

    # 2. Stellar radius
    sel &= np.log10(df['VIIs_R']) < ( (df['VIIs_Teff'] - 5500)/(4000) + 0.2 )

    # 3. Teff=4700-6500K
    sel &= arr(df['VIIs_Teff']) < 6500
    sel &= arr(df['VIIs_Teff']) > 4700

    # 4. Isochrone parallax within 4 sigma of Gaia parallax.
    plxdiff = arr(df.VIIs_plx) - arr(df.VIIs_plxspec)
    plxerr = np.sqrt(df['VIIs_e_plx']**2 + df['VIIs_E_plxspec']**2)
    sel &= (arr(plxdiff/plxerr) < 4)

    # 5. Stellar dilution (Gaia)
    assert np.all(np.isfinite(arr(df.VIIs_r8)))
    sel &= arr(df['VIIs_r8']) < 1.1

    # 6. Stellar dilution (imaging), Furlan+17 radius correction factor <5%, or
    # not given in Furlan+2017 table.
    sel &= ( (arr(df['VIIs_RCF']) < 1.05) | (~np.isfinite(df['VIIs_RCF']) ) )

    # 7. not a false positive
    is_fp = arr(df['II_cks_fp'])
    sel &= ~is_fp

    # 8. not grazing (b<0.9)
    sel &= np.isfinite(arr(df['II_koi_impact']))
    sel &= arr(df['II_koi_impact']) < 0.9

    return sel
