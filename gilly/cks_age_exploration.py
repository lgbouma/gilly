'''
We look at the CKS (paper 7) isochrone ages.

We also get the photometric rotation period measurements of CKS stars.
... and then calculate gyrochronological ages from them (& B-V colors).
(???)
'''
import matplotlib.pyplot as plt, pandas as pd, numpy as np
import os

from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u, astropy.constants as c

def arr(x):
    return np.array(x)


def _get_cks_data(merge_vs_gaia=None):
    '''
    Returns:

        dataframe with CKS II planet data, supplemented by CKS VII stellar
        data, supplemented by Furlan+2017 dilution column.

    Description:

        There are 1305 CKS spectra of KOIs with 2025 planet candidates. The
        overlapping CKS stellar samples are: (1) magnitude-limited, Kp < 14.2,
        (2) multi-planet hosts (3) USP hosts (4) HZ hosts (5) other.  Column
        data are in ../data/cks-column-definitions.txt.

        df_II data are from the CKS website, 2018/07/30, and were released in
        paper II.

        df_VII data are currently (2018/08/04) secret, from CKS paper VII.

        since cks VII table only has star data, I merge it against cks II
        table, using id_starname column, to get planet data too.

        I then also merge against Furlan+2017 to get the dilution column,
        "Avg".
    '''

    df_VII_planets = pd.read_csv(
        '../data/Fulton_2018_CKS_VII_planet_properties.csv')

    df_VII_stars = pd.read_csv(
        '../data/Fulton_2018_CKS_VII_star_properties.csv')

    df_VII_planets['id_starname'] = np.array(
        df_VII_planets['id_koicand'].str.slice(start=0,stop=6)
    )

    df_II = pd.read_csv('../data/cks_physical_merged.csv')

    df_m = pd.merge(df_VII_planets, df_VII_stars, how='left', on='id_starname',
                  suffixes=('_VII_p', '_VII_s'))

    # pull non-changed stellar parameters from CKS II (needed for filters)
    subcols = ['kic_kepmag', 'id_starname', 'cks_fp', 'koi_impact',
               'koi_count', 'koi_dor', 'koi_model_snr', 'id_koicand',
               'id_kic', 'koi_steff', 'koi_slogg', 'koi_dor_err1',
               'koi_dor_err2']

    df = pd.merge(df_m, df_II[subcols], how='left', on='id_koicand',
                  suffixes=('_VII', '_II'))

    assert np.array_equal(arr(df['id_starname_VII']),
                          arr(df['id_starname_II']))
    df['id_starname'] = df['id_starname_VII']
    df = df.drop('id_starname_VII', 1)
    df = df.drop('id_starname_II', 1)

    # compute a/Rstar using "a" from the period and CKS-VII stellar mass, and
    # then divide by the "R" from CKS-VII.  The value of “a/R” that comes from
    # the light curve is pretty ratty.  First, do it directly. Then use the
    # uncertainties package to linearly propagate the errors (take the mean
    # error to make it tractable, otherwise we need to pull out the posteriors,
    # if they're public.)
    period = np.array(df['koi_period'])*u.day
    Mstar = np.array(df['giso_smass'])*u.Msun
    a = ( period**2 * (c.G * Mstar)/(4*np.pi**2) )**(1/3)
    Rstar = np.array(df['giso_srad'])*u.Rsun
    df['cks_VII_dor'] = (a/Rstar).cgs.value

    from uncertainties import unumpy
    period = (np.array(df['koi_period'])*u.day).cgs.value
    period_err = ( np.mean( np.array( [
        np.array(df['koi_period_err1']).astype(float),
        np.abs(np.array(df['koi_period_err2']).astype(float))
        ]), axis=0 )*u.day).cgs.value

    Mstar = (np.array(df['giso_smass'])*u.Msun).cgs.value
    Mstar_err = ( np.mean( np.array( [
        np.array(df['giso_smass_err1']).astype(float),
        np.abs(np.array(df['giso_smass_err2']).astype(float))
        ]), axis=0 )*u.Msun).cgs.value

    Rstar = (np.array(df['giso_srad'])*u.Rsun).cgs.value
    Rstar_err = ( np.mean( np.array( [
        np.array(df['giso_srad_err1']).astype(float),
        np.abs(np.array(df['giso_srad_err2']).astype(float))
        ]), axis=0 )*u.Rsun).cgs.value

    u_period = unumpy.uarray(period, period_err)
    u_Mstar = unumpy.uarray(Mstar, Mstar_err)
    u_Rstar = unumpy.uarray(Rstar, Rstar_err)

    u_a = ( u_period**2 * ((c.G).cgs.value * u_Mstar)/(4*np.pi**2) )**(1/3)

    u_abyRstar = u_a / u_Rstar

    np.testing.assert_array_almost_equal(
        np.array([v.n for v in u_abyRstar]), np.array(df['cks_VII_dor']),
        decimal=8)

    cks_VII_dor_err = np.array([v.s for v in u_abyRstar])
    df['cks_VII_dor_err1'] = cks_VII_dor_err
    df['cks_VII_dor_err2'] = cks_VII_dor_err

    # Get gaia astrometric excess noise significance. see:
    # https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
    # By default, the 1 arcsecond crossmatch leads to some KIC IDs having
    # multiple entries in the Gaia catalog's crossmatch.  We first make
    # `cgk_idm`, the CKS-gaia-kepler ID merge. This has duplicates when there
    # are multiple <1 arcsec crossmatch overlaps with same KICID.  For purposes
    # of assigning a gaia astrometric excess noise significance *ONLY*, we
    # assign the closest Gaia matching star in separation to the given KIC
    # coordinates.

    if merge_vs_gaia:
        gaia_kepler_fun_dir = '/home/luke/local/gaia-kepler_fun_crossmatch/'
        fun_file = 'kepler_dr2_1arcsec.fits'
        gk = Table.read(gaia_kepler_fun_dir+fun_file, format='fits')

        gcols = ['astrometric_excess_noise_sig','kepid','kepler_gaia_ang_dist',
                 'teff', 'logg', 'kepmag', 'phot_g_mean_mag', 'teff_val']
        gkp = gk[gcols].to_pandas()

        cgk_idm = pd.merge(df, gkp, how='left', left_on=['id_kic'],
                           right_on=['kepid'], indicator=True)

        cks_kicids = arr(df['id_kic'])
        cgk_idm_kicids = arr(cgk_idm['id_kic'])

        u_cks, inv_cks, counts_cks = np.unique(cks_kicids,
                                               return_inverse=True,
                                               return_counts=True)

        u_cgk_idm, inv_cgk_idm, counts_cgk_idm = np.unique(cgk_idm_kicids,
                                                           return_inverse=True,
                                                           return_counts=True)

        different_ids = u_cks[counts_cgk_idm - counts_cks != 0]

        # take the closest star in separation (with the same assigned kicid) as the
        # host star, for Gaia astrometric excess noise identification purposes.
        cgk_idm_s = cgk_idm.sort_values(['id_koicand', 'kepler_gaia_ang_dist'])
        cgk_idm_cut = cgk_idm_s.drop_duplicates(subset=['id_koicand'], keep='first')

        del df
        df = cgk_idm_cut

    ##########################################
    ## # I checked that I did the Gaia merge right...
    ## df_scols = ['cks_fp', 'fur17_rcorr_avg', 'koi_steff', 'koi_slogg',
    ##             'kic_kepmag']

    ## baz_scols = ['cks_fp', 'fur17_rcorr_avg', 'koi_steff', 'koi_slogg',
    ##             'kic_kepmag', 'kepler_gaia_ang_dist']

    ## for different_id in different_ids:

    ##     print(df[df['id_kic']==different_id][df_scols])

    ##     print(gkp[gkp['kepid']==different_id])

    ##     print(baz[baz['kepid']==different_id][baz_scols])

    ##     print(50*'#')
    ##     print('\n')
    ##########################################

    # logg is not included, so we must calculate it
    M = np.array(df['giso_smass'])*u.Msun
    R = np.array(df['giso_srad'])*u.Rsun
    g = (c.G * M / R**2).cgs
    df['giso_slogg'] = np.array(np.log10(g.cgs.value))

    f_ages = np.isfinite(df['giso_slogage'])
    f_p_ages = np.isfinite(df['giso_slogage_err1'])
    f_m_ages = np.isfinite(df['giso_slogage_err2'])

    sel = f_ages & f_p_ages & f_m_ages


    df['selected'] = sel

    ##########################################
    ## I previously matched against Furlan17; but 'fur17_rcorr_avg' has the
    ## relevant column already.

    ## # get Furlan+ 2017 table. Format in id_starname in same 0-padded foramt
    ## # as CKS.
    ## furlan_fname = '../data/Furlan_2017_table9.csv'
    ## if not os.path.exists(furlan_fname):
    ##     download_furlan_radius_correction_table()
    ## f17 = pd.read_csv(furlan_fname)
    ## f17['id_starname'] = arr( ['K'+'{0:05d}'.format(koi) for koi in
    ##                            arr(f17['KOI'])] )
    ## # sanity check that above works, i.e. that "id_starname" is just the KOI
    ## # number with strings padded.
    ## assert np.array_equal(arr(df['id_starname']),
    ##                       arr([n[:-3] for n in arr(df['id_koicand'])]) )

    ## df = pd.merge(df, f17, how='left', on='id_starname', suffixes=('cks','f17'))
    ##########################################

    return df


def _apply_cks_IV_metallicity_study_filters(df):
    '''
    given df from _get_cks_data, return the boolean indices for the
    subset of interest.

    (See Petigura+ 2018, CKS IV, table 1)
    '''

    sel = np.isfinite(arr(df['giso_prad']))
    sel &= arr(df['giso_prad']) < 32
    sel &= np.isfinite(arr(df['giso_slogage']))

    # apply all filters from table 1 of Petigura+ 2018
    # Kp<14.2
    sel &= np.isfinite(arr(df['kic_kepmag']))
    sel &= arr(df['kic_kepmag']) < 14.2
    # Teff=4700-6500K
    sel &= np.isfinite(arr(df['cks_steff']))
    sel &= arr(df['cks_steff']) < 6500
    sel &= arr(df['cks_steff']) > 4700
    # logg=3.9-5dex
    sel &= np.isfinite(arr(df['giso_slogg']))
    sel &= arr(df['giso_slogg']) < 5
    sel &= arr(df['giso_slogg']) > 3.9
    # P<350days
    sel &= np.isfinite(arr(df['koi_period']))
    sel &= arr(df['koi_period']) < 350
    # not a false positive
    is_fp = arr(df['cks_fp'])
    sel &= ~is_fp
    # not grazing (b<0.9)
    sel &= np.isfinite(arr(df['koi_impact']))
    sel &= arr(df['koi_impact']) < 0.9
    # radius correction factor <5%, logical or not given in Furlan+2017 table.
    sel &= ( (arr(df['fur17_rcorr_avg']) < 1.05) | (~np.isfinite(df['fur17_rcorr_avg']) ) )

    return sel


def _apply_cks_IV_filters_plus_gaia_astrom_excess(df):
    '''
    also cut on gaia astrometric excess
    '''

    sel = np.isfinite(arr(df['giso_prad']))
    sel &= arr(df['giso_prad']) < 32
    sel &= np.isfinite(arr(df['giso_slogage']))

    # apply all filters from table 1 of Petigura+ 2018
    # Kp<14.2
    sel &= np.isfinite(arr(df['kic_kepmag']))
    sel &= arr(df['kic_kepmag']) < 14.2
    # Teff=4700-6500K
    sel &= np.isfinite(arr(df['cks_steff']))
    sel &= arr(df['cks_steff']) < 6500
    sel &= arr(df['cks_steff']) > 4700
    # logg=3.9-5dex
    sel &= np.isfinite(arr(df['giso_slogg']))
    sel &= arr(df['giso_slogg']) < 5
    sel &= arr(df['giso_slogg']) > 3.9
    # P<350days
    sel &= np.isfinite(arr(df['koi_period']))
    sel &= arr(df['koi_period']) < 350
    # not a false positive
    is_fp = arr(df['cks_fp'])
    sel &= ~is_fp
    # not grazing (b<0.9)
    sel &= np.isfinite(arr(df['koi_impact']))
    sel &= arr(df['koi_impact']) < 0.9
    # radius correction factor <5%, logical or not given in Furlan+2017 table.
    sel &= ( (arr(df['fur17_rcorr_avg']) < 1.05) | (~np.isfinite(df['fur17_rcorr_avg']) ) )

    # of the 1944 objects in the initial objects, 1606 have non-zero
    # astrometric excess noise values.
    # the associated distribution stats are:
    #
    # In [7]: df[df['astrometric_excess_noise_sig']!=0]['astrometric_excess_noise_sig'].describe()
    # Out[7]: 
    # count    3.380000e+02
    # mean     3.710397e+02
    # std      1.164272e+03
    # min      6.486338e-16
    # 25%      1.330256e-15
    # 50%      6.488688e-01
    # 75%      3.895363e+01
    # max      9.709389e+03
    # Name: astrometric_excess_noise_sig, dtype: float64
    # 
    # so we'll follow Rizzuto+ 2018 (ZEIT VIII) and use D > 10sigma as the
    # cutoff.
    sel &= (
            ( arr(df['astrometric_excess_noise_sig']) >= 0 ) &
            ( arr(df['astrometric_excess_noise_sig']) < 10 )
    )

    return sel
