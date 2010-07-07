#root -q QCDPure.cc

# perform fit for gaus and pol3 with various fit ranges
root -q run_fit.cc >& fit_res.txt

# variation of est given dff fit ranges
root -q 'plot_qcd_estimate.cc("landau",1)'
#root -q 'plot_qcd_estimate.cc("gaus",0)'
#root -q 'plot_qcd_estimate.cc("pol3",1)'
#root -q 'plot_qcd_estimate.cc("pol3",0)'

#animate
so animate_landau.sh
