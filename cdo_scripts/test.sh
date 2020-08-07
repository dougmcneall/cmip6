#!/bin/bash


#ls /data/users/hadaw/cmip6/areaavg/tas/*.nc

for i in /data/users/hadaw/cmip6/areaavg/tas/*historical_r1i1p1*;
do
echo $i
cdo -ntime $i
#ncdump -t $i
#cdo remapbil,target.grd  "$i" > "tas_${i}_bil.nc"
done


#cdo -seltimestep,154/165 /data/users/hadaw/cmip6/areaavg/tas/tas_Amon_UKESM1-0-LL_historical_r1i1p1f2_gn_area_avg.nc ./test.nc
