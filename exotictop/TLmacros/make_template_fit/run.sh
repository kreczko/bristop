
## use EME+BCE template
#########################
root -q -b 'do_fit.C' >& fit_res.txt

#source animate_plots.sh
rm -rf  template_from_enriched_qcd
mkdir template_from_enriched_qcd

for dir in \
    B3_e20_1T_1j \
    B3_e20_1T_2j \
    B3_e20_1T_3j \
    B3_e20_1T_4mj \
    B3_e20_4T ;
do
    echo moving $dir
    mkdir -p $dir
    mv -f *$dir*.* $dir
    mv $dir  ./template_from_enriched_qcd/
done

mv fit_res.txt ./template_from_enriched_qcd/


## use dijet template
#######################
root -q -b 'do_fit.C("dijet")' >& fit_res_dijet.txt


for dir in B3_e20_1T_1j
do
    echo moving $dir
    mkdir -p $dir
    mv -f *$dir*.* $dir
done

rm -rf template_from_dijet_pt15
mkdir template_from_dijet_pt15
mv B3_e20_1T_1j  ./template_from_dijet_pt15/
mv fit_res_dijet.txt   ./template_from_dijet_pt15/

echo done
