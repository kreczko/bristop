
# make subdirectories
mkdir nofit
mkdir gaus_free12j_fix34j
mkdir gaus_mean_0.3to0.6_12j_fix34j
mkdir gaus_mean_0.4to0.6_12j_fix34j


# plot without fit
echo "plot without fit"
root -q -b run_fit_none.cc
mv -f *QCD_reliso*.pdf ./nofit/



# Fit: gaus, no constrain in 1,2j, constrain in 3,4j
echo "Fit gaus"
root -q -b run_fit.cc >& fit_res_gaus.txt
so animate_gaus.sh
mv -f *pdf *txt *gif ./gaus_free12j_fix34j/


# Fit gaus, constrain mean in 1,2j mean:0.3-1.6, fix 3,4j
echo "Fit gaus, constrain mean in 1,2j (0.3-1.6)"
root -q -b run_fit_mean12j_geq03.cc >& fit_res_mean12j_geq03.txt
so animate_gaus.sh
mv -f *pdf *txt *gif ./gaus_mean_0.3to0.6_12j_fix34j/


# Fit gaus, constrain mean in 1,2j mean:0.4-1.6, fix 3,4j
echo "Fit gaus, constrain mean in 1,2j (0.4-1.6)"
root -q -b run_fit_mean12j_geq04.cc >& fit_res_mean12j_geq04.txt
so animate_gaus.sh
mv -f *pdf *txt *gif ./gaus_mean_0.4to0.6_12j_fix34j/

echo done
