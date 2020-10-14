o=../results/best_20201014/

if [ ! -d $o ]
then
  echo 'making' $o
  mkdir $o
fi

cp ../results/cks_M15_rotation_period/cks-VII_period_vs_cks-VII_prad_M15_match.png $o/00.png
cp ../results/cks_M15_rotation_period/cks-VII_period_vs_cks-VII_prad_M15-Protcolors.png $o/01.png
cp ../results/cks_M15_rotation_period/cks-VII_rp_vs_M15_prot_log.png $o/02.png

cp ../results/M15rot_A19gyroage/cks-VII_period_vs_cks-VII_prad_gyroagecolors.png $o/03a.png
cp ../results/M15rot_A19gyroage/cks-VII_rp_vs_gyroage_log.png $o/04a.png

#cp ../results/cdf_rp_agecut_KS_test/cdf_Rp_KS_cut0.7_protM15_gyroA19.png $o/06.png
cp ../results/cdf_rp_agecut_KS_test/cdf_Rp_KS_cut1.0_protM15_gyroA19.png $o/05a.png

# cp ../results/hist_rp_agecut/hist_Rp_cut0.7_protM15_gyroA19.png $o/08.png
cp ../results/hist_rp_agecut/hist_Rp_cut1.0_protM15_gyroA19.png $o/06a.png

cp ../results/singletransitSNR_M15rot_A19age/singletransitsnr_rms3_rotM15_gyroA19.png $o/07.png

cp ../results/M15rot_SL20gyroage/cks-VII_period_vs_cks-VII_prad_gyroagecolors.png $o/03b.png
cp ../results/M15rot_SL20gyroage/cks-VII_rp_vs_gyroage_log.png $o/04b.png
cp ../results/cdf_rp_agecut_KS_test/cdf_Rp_KS_cut0.8_protM15_gyroSL20.png $o/05b.png
cp ../results/hist_rp_agecut/hist_Rp_cut0.8_protM15_gyroSL20.png $o/06b.png

