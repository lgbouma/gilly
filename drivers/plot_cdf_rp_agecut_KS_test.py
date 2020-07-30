import numpy as np
from gilly.plotting import plot_cdf_rp_agecut_KS_test

for a in 1e9*np.arange(0.7,1.6,0.1):
    plot_cdf_rp_agecut_KS_test(agecut=a, Prot_source='M15', gyro_source='A19')
