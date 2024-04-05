#!/bin/bash

# get the data

ddir=data/ams/
mkdir -p $ddir
curl -L https://www-air.larc.nasa.gov/img/tmp/WWW-AIR_1712276024307.zip -o ${ddir}/rawdata.zip
unzip ${ddir}/rawdata.zip -d ${ddir}

ddir=data/smoke_flag/
mkdir -p $ddir
curl -L https://www-air.larc.nasa.gov/img/tmp/WWW-AIR_1712278888480.zip -o ${ddir}/rawdata.zip
unzip ${ddir}/rawdata.zip -d ${ddir}

#conda activate jpy01
python run/firexaq2cmaq.py 20190801 20190807 108NHEMI2 10

for f in ./gridded/*.nc
do
    nccopy -7 -d 1 -c TSTEP/1,LAY/1,ROW/187,COL/187 $f ${f%.*}.compressed.nc
done

for f in ./gridded/*.compressed.nc
do
    echo $f ${f%.*.*}.nc
    mv -n $f ${f%.*.*}.nc
done
