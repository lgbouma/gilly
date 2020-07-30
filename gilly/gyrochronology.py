import numpy as np

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
