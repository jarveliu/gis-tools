#! /bin/bash

#extract-os-terrain.bash - extract a rect dtm from Ordinance Service data set.
#
# Usage: extract-os-terrain.bash llx lly urx ury
# coordinates as in british_national_grid_regions.bash
# outputs to files dtm.tiff and dtm_100.tiff in different resolutions
#
# NB! Beware of bad coding with hardcoded filenames and search path for data.
# Depends on british_national_grid_regions.bash
# 
# Johan Arvelius, SMHI, 2017-01-17

terrainfiledir=~/gisslask/UK/os_terrain_50/data/

function filenames()
{
    for reg in $(british_national_grid_regions.bash $1 $2 $3 $4)
    do
	echo $reg.asc
    done
}

function providefile()
{
    lf=${1,,}
    filename="${lf%.asc}*.zip"
    if [ ! -f "$filename" ]
    then
	unzip $terrainfiledir${lf:0:2}/$filename
    fi
}

for f in $(filenames $1 $2 $3 $4)
do
    if [ ! -f "$f" ]
    then
	providefile $f
    fi
done


gdal_merge.py $(filenames $1 $2 $3 $4) -o dtm.tiff -of GTiff -ul_lr $1 $4 $3 $2
gdalwarp dtm.tiff dtm_100.tiff -tr 100 100
