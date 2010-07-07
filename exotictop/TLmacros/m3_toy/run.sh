cp ../../FullRun_btag_new_mc_mixture.root .
root -q -b run_m3.C >& m3.txt
root -q -b 'plot_nfit_pull.C'
