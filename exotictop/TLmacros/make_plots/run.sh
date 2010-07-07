
root -q 'plot_njet.C'

root -q 'plot_reliso_01to4j_rev1.C("data")'
root -q 'plot_reliso_01to4j_rev1.C("QCD")'

root -q 'plot_iso_aes_compo_nes_3in1.C'

root -q 'plot_iso_QCD_sig_reg_superimpose_1234j.C'

root -q 'plot_QCD_signal_control.C'


# plot iso, met, mtw, dphi
root -q 'plot_vars_ind.C'


root -q 'plot_explore_events_sel.C'
root -q 'plot_explore2D_events_sel.C'

root -q 'plot_planB.C'

 
### plot reliso NES/AES in the same y scale
cd QCD_iso_setyscale
root -q 'plot_reliso_01to4j_rev1.C'
cd ..


### QCD reliso poisson error
cd QCD_reliso_poisson_error
root -q 'plot_reliso_poisson_error.C'
cd ..


### ele/muon variables
root -q -b 'plot_lepton_var.C'


### correlation iso:met
#########################
cd QCD_correlation_iso_met
#mkdir both_barrel_plus_endcap
#mkdir barrel_only
#mkdir endcap_only
echo "project iso from isoVmet"
root -q -b 'plot_isoVmet_iso.C("both","t1")'
root -q -b 'plot_isoVmet_iso.C("barrel","t1")'
root -q -b 'plot_isoVmet_iso.C("endcap","t1")'
echo "project met from isoVmet"
root -q -b 'plot_isoVmet_met.C("both","t1")'
root -q -b 'plot_isoVmet_met.C("barrel","t1")'
root -q -b 'plot_isoVmet_met.C("endcap","t1")'
cd ..

echo m3
cd m3
root -q -b 'plot_m3_stack_norm_2in1.cc'
root -q -b 'plot_m3_norm_wzstop.cc'
root -q -b 'plot_m3_show_QCD_dominant_in_controlreg.cc'
cd ..

echo "done :)"
