#!/bin/csh -f
#USAGE
# lnpirate.csh piratehome
# e.g. lnpirate.csh /home/hwang/whsoft/pirate/v2.0
# Hua Wang, Nov. 23, 2010

if ($#argv<1) then
  set pirate = `echo $PIRATE`
else
  set pirate = $1
endif

cat << ENDLIST >! piratelist
extractfaultonprof.m
extractprof.m
getpar.m
gmt2mat_faults.m
hdr2rsc.m
ifghdrlooks.m
imagecrop.m
ll2utm.m
ll2xy.m
looks.m
make_prof.m
make_vcms.m
make_vcms_ratemap.m
matcrop.m
plotifg.m
plotifgs.m
profdefine.m
profll2xy.m
readmat.m
readparfile.m
rsc2hdr.m
writemat.m
cvdcalc.m
pendiffexp.m
make_vcm.m
ENDLIST

foreach func (`cat piratelist`)
  rm -f ./pilib/$func
  ln -s $pirate/$func ./pilib
end
rm -f piratelist
