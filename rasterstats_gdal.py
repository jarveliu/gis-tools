#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy
import gdal

from gdalconst import GA_ReadOnly

class Raster:
    def __init__(self, su=None, infile=None, nodatavalue=None):
        if infile==None:
            self.infile = su.infile
        else:
            self.infile = infile

        if nodatavalue == None:
            self.nodatavalue = su.nodatavalue
        else:
            self.nodatavalue = nodatavalue

        self.dataset = gdal.Open(self.infile, GA_ReadOnly)
        if self.dataset == None:
            raise IOError("No such file "+self.infile)

        self.band = self.dataset.GetRasterBand(1)
        self.arr = self.band.ReadAsArray()
        if self.nodatavalue == None:
            self.raster = self.arr
        else:
            self.raster = numpy.ma.masked_values(self.arr,self.nodatavalue)
        
    def sum(self):
        return self.raster.sum()

    def mean(self):
        return self.raster.mean()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--nodatavalue', type=int)
    su = parser.parse_args()
    rast = Raster(su=su)
    print rast.sum(), rast.mean()

if __name__ == "__main__":
    main()
