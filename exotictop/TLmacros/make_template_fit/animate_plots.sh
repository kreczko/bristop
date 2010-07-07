# produce animated pictures of QCD fits
for i in \
    B3_e20__1T4j \
#    B3_e20__4T \
#    B2_e20__1T \
#    B2_e20__4T;
do
    convert -delay 100 \
	MLfit__${i}__r0.2to1.0.gif \
	MLfit__${i}__r0.2to1.2.gif \
	MLfit__${i}__r0.2to1.4.gif \
	MLfit__${i}__r0.3to1.1.gif \
	MLfit__${i}__r0.3to1.3.gif \
	MLfit__${i}__r0.3to1.5.gif \
	MLfit__${i}__r0.4to1.2.gif \
	MLfit__${i}__r0.4to1.4.gif \
	MLfit__${i}__r0.4to1.6.gif \
	-loop 0 play_${i}.gif
done
