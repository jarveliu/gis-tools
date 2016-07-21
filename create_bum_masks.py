#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys, ogr
import numpy as np
from rasterize import BUMMask, Raster, rasterize_parser

def main():
    parser = argparse.ArgumentParser(parents=[rasterize_parser,BUMMask.parser])
    parser.add_argument('--rastergeometry', nargs=6,
                        default=[222000,7674000,1000,-1000,740,1570],
                        help='origin_x origin_y step_x step_y No_of_pixels_x No_of_pixels_y')
    su = parser.parse_args()
    driver = ogr.GetDriverByName('ESRI Shapefile')
    inDataSource = driver.Open(su.infile,0)
    inLayer = inDataSource.GetLayer()
    myRaster = Raster(su=su,srs=inLayer.GetSpatialRef())

    for i in range(inLayer.GetFeatureCount()):
        inFeature = inLayer.GetFeature(i)
        geom = inFeature.GetGeometryRef()
        myRaster.updatePixelsShape(inFeature)

    myMask = BUMMask(su, myRaster)
    for setid,pixset in myRaster.pixelsets.iteritems():
        myMask.addPixelSet(pixset)
    myMask.createMasks()
    myMask.printmasks()

if __name__ == "__main__":
    main()
