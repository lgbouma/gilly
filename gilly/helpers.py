"""
Contents:

    get_merged_rot_CKS
    get_merged_gyroage_CKS

    _get_spec_vsini_phot_Prot_overlap

    _add_Christiansen12_CDPP
    _get_McQuillan13_data
    _get_Mazeh15_data
    _get_Santos21_data
    _get_cks_data
    _apply_cks_VII_filters
    _get_Curtis20_data

    PleiadesInterpModel
    PraesepeInterpModel
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
from scipy.interpolate import interp1d

from gilly.paths import DATADIR, RESULTSDIR
from gilly.gyrochronology import (
    MamajekHillenbrand08_gyro, Angus19_gyro, SpadaLanzafame20_gyro
)
from cdips.utils.mamajek import get_interp_BmV_from_Teff

def _get_spec_vsini_phot_Prot_overlap(Prot_source='M15'):

    rdf, fp18_df = get_merged_rot_CKS(Prot_source=Prot_source)

    cdf = _get_CKS_vsini_data()

    cdf.id_kic = cdf.id_kic.astype(np.int64)
    rdf.II_id_kic = rdf.II_id_kic.astype(np.int64)

    mdf = rdf.merge(cdf, how='left', left_on='II_id_kic', right_on='id_kic')
    assert len(mdf) == len(rdf)

    # vsini ~= 2πR*/Prot (when sini=1)
    vsini = arr(mdf.cks_svsini)*u.km/u.s
    rstar = arr(mdf.VIIs_R)*u.Rsun
    spec_Prot = (2*np.pi*rstar / vsini).to(u.day).value
    mdf['spec_Prot'] = spec_Prot

    return mdf, cdf


def get_merged_rot_CKS(Prot_source='M15'):

    df = _get_cks_data()
    fp18_df = df[_apply_cks_VII_filters(df)]
    if Prot_source == 'M13':
        rot_df = _get_McQuillan13_data()
        kic_key = 'KIC'
    elif Prot_source == 'M15':
        rot_df = _get_Mazeh15_data()
        kic_key = 'KIC'
    elif Prot_source == 'VSINI':
        rot_df = _get_CKS_vsini_data(snrcut=True)
        kic_key = 'id_kic'
    elif Prot_source == 'S21':
        rot_df = _get_Santos21_data()
        kic_key = 'KIC'
    else:
        raise NotImplementedError

    fp18_df['II_id_kic'] = fp18_df['II_id_kic'].astype(str)
    rot_df[kic_key] = rot_df[kic_key].astype(str)

    mdf = fp18_df.merge(rot_df, how='inner', left_on='II_id_kic',
                        right_on=kic_key)

    return mdf, fp18_df


def get_merged_gyroage_CKS(Prot_source='M15', gyro_source='MH08'):
    """
    Same as get_merged_rot_CKS, but also calculate the gyro ages using the
    Mamajek & Hillenbrand (2008) relation (or the Angus+19 one, or the Spada &
    Lanzafame 2020 one).
    """

    mdf, _ = get_merged_rot_CKS(Prot_source=Prot_source)
    Prot_key = 'Prot' if Prot_source in ['M13','M15','S21'] else 'spec_Prot'
    Prot = arr(mdf[Prot_key])

    if gyro_source == 'MH08':
        BmV = get_interp_BmV_from_Teff(arr(mdf.VIIs_Teff))
        mdf['BmV'] = BmV
        t_yr = MamajekHillenbrand08_gyro(arr(mdf.BmV), Prot)
        mdf['gyroage_yr'] = t_yr

    elif gyro_source == 'SL20':
        BmV = get_interp_BmV_from_Teff(arr(mdf.VIIs_Teff))
        mdf['BmV'] = BmV
        t_yr = SpadaLanzafame20_gyro(arr(mdf.BmV), Prot)
        mdf['gyroage_yr'] = t_yr
        mdf = mdf[mdf.gyroage_yr != 10e9]

    elif gyro_source == 'A19':

        gaiasupppath = join(
            DATADIR, f'prot{Prot_source}_gyro{gyro_source}_Gaiasupp.csv'
        )
        if not os.path.exists(gaiasupppath):
            # merge the CKS + rotation KOIs against Gaia. Use Megan Bedell's
            # https://gaia-kepler.fun/ 1 arcsecond crossmatch.
            t = Table.read(
                join(DATADIR, 'kepler_dr2_1arcsec.fits'), format='fits'
            )
            g_df = t.to_pandas()

            g_df['KIC'] = g_df.kepid.astype(str)

            selcols = ['KIC', 'source_id', 'astrometric_excess_noise',
                       'phot_g_mean_mag', 'phot_bp_mean_mag',
                       'phot_rp_mean_mag']

            g_df = g_df[selcols]

            kic_key = 'KIC' if Prot_source in ['M13','M15','S21'] else 'id_kic'
            mdf[kic_key] = mdf[kic_key].astype(str)

            # NOTE: Bedell's crossmatch was KIC against Gaia DR2. So, some KIC
            # entries have multiple source_ids. (i.e., chance alignments or
            # binaries within 1 arcsec). Usually the bright star dominates --
            # take that one.
            g_df = g_df.sort_values(by=['KIC','phot_g_mean_mag'])
            g_df = g_df.drop_duplicates(subset='KIC', keep='first')

            mdf_gaiasupp = pd.merge(mdf, g_df, left_on=kic_key, right_on='KIC', how='left')
            print(gaiasupppath)
            print(f'Pre merge: {len(mdf)}')
            print(f'Gaia supp: {len(mdf_gaiasupp)}')

            mdf_gaiasupp.to_csv(gaiasupppath, index=False)

        if Prot_source in ['M13', 'M15', 'S21', 'VSINI']:
            N_missed = 0
        else:
            raise NotImplementedError

        mdf = pd.read_csv(gaiasupppath)
        print(42*'-'+f'\nWRN! Missed {N_missed} sources during Gaia merge\n'+42*'-')

        BpmRp = arr(mdf.phot_bp_mean_mag) - arr(mdf.phot_rp_mean_mag)
        mdf['BpmRp'] = BpmRp

        t_yr = Angus19_gyro(BpmRp, Prot)
        mdf['gyroage_yr'] = t_yr

    return mdf


def _add_Christiansen12_CDPP(mdf):

    t = Table.read(
        join(DATADIR, 'Christiansen_2012_Kepler_CDPP_table1.fits'), format='fits'
    )
    c_df = t.to_pandas()

    # duplicates have different program IDs (PID), either "STKS" or "EX_STKS".
    # not important! all stellar params are the same.
    c_df = c_df.drop_duplicates(subset='KIC', keep='first')

    outdf = mdf.merge(c_df, how='left', on='KIC')

    assert len(outdf) == len(mdf)

    return outdf



def _get_CKS_vsini_data(snrcut=False):
    """
    if snrcut: require vsini > 3 km/s
    """

    cdf = pd.read_csv(join(DATADIR, '20200806_cksgaia-planets-filtered.csv'))

    if snrcut:
        cdf = cdf[cdf.cks_svsini > 3]

    vsini = arr(cdf.cks_svsini)*u.km/u.s
    rstar = arr(cdf.giso_srad)*u.Rsun
    spec_Prot = (2*np.pi*rstar / vsini).to(u.day).value
    cdf['spec_Prot'] = spec_Prot

    return cdf



def _get_McQuillan13_data():

    vot = parse(join(DATADIR,'McQuillan_2013_ApJ_775L.vot'))
    m13 = vot.get_first_table()
    m13_df = m13.to_table().to_pandas()

    return m13_df



def _get_Mazeh15_data():

    hdul = fits.open(join(DATADIR,'Mazeh_2015_ApJ_801_3_table1.fits'))

    t = Table(hdul[1].data)

    m15_df = t.to_pandas()

    # C: Centroid motion shows transit and rotational modulation on different
    # stars 
    # F: false positive (EB) identified in Mazeh+15
    # R: 1 if rejected by visual exmaination stage
    # M1: inconsistent periods across quarters
    # M2: peaks not high enough depending on temperature and period
    # 
    # ignore:
    # G: probable Giant flag, because the CKS Gaia parameters already select
    # for this.
    # T: Teff outside 3500-6500K, because the CKS Gaia parameters already select
    # for this.
    # 
    # 

    print(f'Mazeh15 starting with {len(m15_df)} KOIs')

    sel = (
        (~m15_df.C.astype(bool)) &
        (~m15_df.F.astype(bool)) &
        (~m15_df.R.astype(bool)) &
        (~m15_df.M1.astype(bool)) &
        (~m15_df.M2.astype(bool))
    )

    print(f'Mazeh15 with {len(m15_df[sel])} KOIs post-selection')

    return m15_df[sel]


def _get_Santos21_data():

    t = Table.read(
        join(DATADIR, 'Santos_2021_apjsac033ft1_mrt.txt'),
        format='cds'
    )

    df = t.to_pandas()

    return df


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

    df_II = pd.read_csv(os.path.join(DATADIR, 'cks_physical_merged.csv'))
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


def _get_Curtis20_data(cluster):
    """
    cluster == "pleiades" or "praesepe"
    """

    t = Table.read(
        join(DATADIR, 'Curtis_2020_apjabbf58t5_mrt.txt'),
        format='cds'
    )
    if cluster == 'pleiades':
        df = t[t['Cluster'] == 'Pleiades'].to_pandas()
    elif cluster == 'praesepe':
        df = t[t['Cluster'] == 'Praesepe'].to_pandas()
    else:
        raise NotImplementedError

    return df


def PleiadesInterpModel(Teff):

    df = pd.read_csv(
        join(DATADIR, 'Curtis_2020_webplotdigitize_Pleiades_upper.csv')
    )
    df = df.sort_values(by='Teff')

    fn_teff_to_prot = interp1d(
        df.Teff, df.Prot, kind='quadratic',
        bounds_error=True, fill_value=np.nan
    )

    Protmod = fn_teff_to_prot(Teff)

    return Protmod


def PraesepeInterpModel(Teff):

    df = pd.read_csv(
        join(DATADIR, 'Curtis_2020_webplotdigitize_Praesepe_mid.csv')
    )
    df = df.sort_values(by='Teff')

    fn_teff_to_prot = interp1d(
        df.Teff, df.Prot, kind='quadratic',
        bounds_error=True, fill_value=np.nan
    )

    Protmod = fn_teff_to_prot(Teff)

    return Protmod
