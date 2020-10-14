import numpy as np
from gilly.plotting import plot_hist_rp_agecut

for a in 1e9*np.arange(0.7,1.6,0.1):
    plot_hist_rp_agecut(agecut=a, Prot_source='VSINI', gyro_source='SL20',
                        fix_ylim=1)
    plot_hist_rp_agecut(agecut=a, Prot_source='VSINI', gyro_source='A19',
                        fix_ylim=1)

for a in 1e9*np.arange(0.7,1.6,0.1):
    plot_hist_rp_agecut(agecut=a, Prot_source='M15', gyro_source='SL20',
                        fix_ylim=1)
    plot_hist_rp_agecut(agecut=a, Prot_source='M15', gyro_source='A19',
                        fix_ylim=1)
