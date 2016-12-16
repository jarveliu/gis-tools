#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import ogr
from rasterize import Raster, rasterize_parser

def createRasterFromPolygons(su,inLayer):
    myRaster = Raster(su=su, srs=inLayer.GetSpatialRef())
    for i in range(inLayer.GetFeatureCount()):
        inFeature = inLayer.GetFeature(i)
        myRaster.updatePixelsShape(inFeature)
    myRaster.updateMatrix(su.field)
    return myRaster

def createRasterFromPoints(su,inLayer):
    myRaster = Raster(su=su, srs=inLayer.GetSpatialRef())
    myRaster. updateMatrixFromPointsLayer(inLayer,su.field)
    return myRaster

def main():
    parser = argparse.ArgumentParser(parents=[rasterize_parser, Raster.parser])
    parser.add_argument('field')
    su = parser.parse_args()

    driver = ogr.GetDriverByName('ESRI Shapefile')
    inDataSource = driver.Open(su.infile,0)
    inLayer = inDataSource.GetLayer()
    
    inFeature = inLayer.GetFeature(0)
    geom = inFeature.GetGeometryRef()
    if geom.GetGeometryName() == 'POINT':
        myRaster = createRasterFromPoints(su,inLayer)
    elif geom.GetGeometryName() == 'POLYGON':
        myRaster = createRasterFromPolygons(su,inLayer)

    myRaster.printFile()

if __name__ == "__main__":
    main()
