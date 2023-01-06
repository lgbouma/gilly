"""
Given color and Prot, what is age?
    MamajekHillenbrand08_gyro(BmV, Prot)
    SpadaLanzafame20_gyro(BmV=None, Prot=None, makeplot=True)
    Angus19_gyro(BpmRp, Prot)

Given color and age, what is Prot?
    MamajekHillenbrand08_gyro_Prot
    Angus19_gyro_Prot
    SpadaLanzafame20_gyro_Prot
"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib as mpl
from gilly.paths import DATADIR, RESULTSDIR

def MamajekHillenbrand08_gyro_Prot(BmV, t_Myr):

    if BmV > 0.95 or BmV < 0.5:
        BmV = np.nan

    a = 0.407
    a_err = 0.021
    b = 0.325
    b_err = 0.024
    c = 0.495
    c_err = 0.010
    n = 0.566
    n_err = 0.008

    Prot = a * (BmV - c)**b * (t_Myr)**n

    return Prot


def MamajekHillenbrand08_gyro(BmV, Prot):
    """
    Mamajek & Hillenbrand 08 used

        Prot = a * [ (B-V)_0 - c ]^{b} * (t/Myr)^{n}

    and fitted some numbers. So

        (t/Myr) = [  (Prot/a) / [ (B-V)_0 - c ]^{b} ]^{1/n}.

    This function returns age given Prot and B-V.
    """

    a = 0.407
    a_err = 0.021
    b = 0.325
    b_err = 0.024
    c = 0.495
    c_err = 0.010
    n = 0.566
    n_err = 0.008

    t_Myr = ( (Prot / a) / ( (BmV - c)**(b) ) )**(1/n)

    t_yr = t_Myr * 1e6

    return t_yr


def find_nearest(array, value, verbose=False, return_index=False):

    idx = (np.abs(np.array(array)-value)).argmin()

    if verbose:
        print(f"{value:.5f} --> {array[idx]:.5f}")

    if not return_index:
        return array[idx]

    if return_index:
        return array[idx], idx


def SpadaLanzafame20_gyro_Prot(BmV, age_myr):

    from scipy.interpolate import interp2d

    sl20path = os.path.join(DATADIR, 'Spada_Lanzafame_2020_BmV_table.csv')
    df = pd.read_csv(sl20path, sep='&')

    df = df.sort_values(by='BmV')
    BmV_arr = np.array(df.BmV).astype(float)
    age_gyr = np.array(df.columns[1:-1]).astype(float)
    Prot_array = np.array(
        df[[c for c in df.columns if (c!='BmV' and c!='4.57')]]
    ).astype(float)

    fn = interp2d(age_gyr, BmV_arr, Prot_array, kind='cubic', bounds_error=False,
                  fill_value=0)

    Prot_new = fn(age_myr/(1e3), BmV)

    if np.sum(Prot_new) == 0:
        Prot_new == np.nan

    return float(Prot_new)


def SpadaLanzafame20_gyro(BmV=None, Prot=None, makeplot=True):
    """
    Spada & Lanzafame 2020 give gyrochrones in B-V (and mass) vs age space.
    Will need to interpolate between them on a 2D grid.
    """

    sl20path = os.path.join(DATADIR, 'Spada_Lanzafame_2020_BmV_table.csv')
    df = pd.read_csv(sl20path, sep='&')

    from scipy.interpolate import interp2d, interp1d
    from scipy.optimize import newton

    df = df.sort_values(by='BmV')
    BmV_arr = np.array(df.BmV).astype(float)
    age_gyr = np.array(df.columns[1:-1]).astype(float)
    Prot_array = np.array(
        df[[c for c in df.columns if (c!='BmV' and c!='4.57')]]
    ).astype(float)

    fn = interp2d(age_gyr, BmV_arr, Prot_array, kind='cubic', bounds_error=False,
                  fill_value=0)

    N_pts = 100
    BmV_new = np.linspace(min(BmV_arr), max(BmV_arr), N_pts)
    age_new = np.linspace(min(age_gyr), max(age_gyr), 2*N_pts)
    Prot_new = fn(age_new, BmV_new)

    if makeplot:
        from aesthetic.plot import set_style, savefig
        set_style()
        fig, ax = plt.subplots(figsize=(4,3))
        X,Y = np.meshgrid(BmV_new, age_new)
        #cmap = plt.cm.viridis
        cmap = plt.cm.Paired
        norm = mpl.colors.Normalize(vmin=np.nanmin(Prot_new),
                                    vmax=np.nanmax(Prot_new))
        im = ax.pcolormesh(X, Y, Prot_new.T, cmap=cmap, norm=norm,
                           shading='flat') # vs 'gouraud'
        ax.set_xlabel('B-V')
        ax.set_ylabel('Age [Gyr]')
        cbar = fig.colorbar(im, orientation='vertical', extend='min',
                            label='Prot [day]')
        cbar.ax.tick_params(labelsize=6, direction='out')
        outpath = os.path.join(RESULTSDIR,'SpadaLanzafame2020',
                               'interpolation_check.png')
        savefig(fig, outpath, writepdf=False)

    age_gyrs = []

    for _prot, _bmv in zip(Prot, BmV):

        nearest_BmV, idx = find_nearest(BmV_new, _bmv, verbose=True,
                                        return_index=True)

        # now it's a 1-dimensional optimization problem
        fn_1d = interp1d(age_new, Prot_new[idx, :], kind='quadratic',
                         bounds_error=False, fill_value=0)

        interp_fn2 = lambda x: fn_1d(x) - _prot

        age_init_gyr = 2
        try:
            root_age_gyr = newton(interp_fn2, age_init_gyr)
        except:
            if _prot > max(Prot_new[idx, :]):
                root_age_gyr = 10
            else:
                root_age_gyr = 0

        age_gyrs.append(root_age_gyr)

    t_yr = np.array(age_gyrs)*1e9

    return t_yr




def Angus19_gyro(BpmRp, Prot):
    """
    Run a ~*SIMPLE*~ gyrochronology model from Angus+2019.
    Requires `stardate`: https://github.com/RuthAngus/stardate

    If you want more ~*PRECISE*~ ages, it is better to use spectroscopic
    parameters for the joint isochrone + gyrochrone fit, and to run MCMC. See
    the README from the repo above.
    """

    from stardate.lhf import age_model

    log10_period = np.log10(Prot)

    log10_age_yrs = []
    for p, _bpmrp in zip(log10_period, BpmRp):
        log10_age_yrs.append(age_model(p, _bpmrp))

    t_yr = 10**np.array(log10_age_yrs)

    return t_yr

def Angus19_gyro_Prot(BpmRp, t_Myr):
    """
    Run a ~*SIMPLE*~ gyrochronology model from Angus+2019.
    Requires `stardate`: https://github.com/RuthAngus/stardate

    If you want more ~*PRECISE*~ ages, it is better to use spectroscopic
    parameters for the joint isochrone + gyrochrone fit, and to run MCMC. See
    the README from the repo above.
    """

    from stardate.lhf import gyro_model_praesepe

    log10_age_yrs = np.log10(t_Myr*1e6)
    log10_period = gyro_model_praesepe(log10_age_yrs, BpmRp)

    return 10**log10_period
