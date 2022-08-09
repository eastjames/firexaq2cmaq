#!/bin/bash

#conda activate jpy01
python firexaq2cmaq.py 20190801 20190807 108NHEMI2 10

for f in ./gridded/*.nc
do
    nccopy -7 -d 1 -c TSTEP/1,LAY/1,ROW/187,COL/187 $f ${f%.*}.compressed.nc
done

for f in ./gridded/*.compressed.nc
do
    echo $f ${f%.*.*}.nc
    mv -n $f ${f%.*.*}.nc
done
