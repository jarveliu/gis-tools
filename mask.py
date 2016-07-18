#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ogr, os, sys, osr
import numpy as np
from math import floor

class Pixelset:
    def __init__(self,raster):
        self.pixels = dict()
        self.raster = raster
        #                   rowmin, rowmax, colmin, colmax
        self.boundingbox = [raster.ny, -1, raster.nx, -1]

    def getPixel(self,rowcol):
        row = rowcol[0]
        col = rowcol[1]
        if row not in self.pixels.keys():
            self.pixels[row] = dict()
            
        if col in self.pixels[row].keys():
            pix = self.pixels[row][col]
        else: 
            pixcoords = self.raster.calcPixelCoordinates(row,col)
            pix = Pixel(row,col,pixcoords)
            self.addPixel(pix,row,col)
        return pix
    
    def addPixel(self,pixel,row,col):
        if row not in self.pixels.keys():
            self.pixels[row] = dict()
        if col in self.pixels[row].keys():
            raise exception

        self.pixels[row][col] = pixel

        if row < self.boundingbox[0]: 
            self.boundingbox[0] = row
        if row > self.boundingbox[1]: 
            self.boundingbox[1] = row
        if col < self.boundingbox[2]:
            self.boundingbox[2] = col
        if col > self.boundingbox[3]:
            self.boundingbox[3] = col

    def addEnclosedPixels(self, geom):
        candidates = dict()
        #import pdb; pdb.set_trace()
        for i in range(self.boundingbox[0]+1, self.boundingbox[1]):
            if i in self.pixels.keys():
                keys = self.pixels[i].keys()
                cands = list()
                prel = list()
                for col in range(keys[0],keys[-1]):
                    if col+1 not in keys:
                        prel.append(col)
                    elif prel:
                        for val in prel:
                            cands.append(val)
                        prel = list()
                for cand in cands:
                    pixcoords = self.raster.calcPixelCoordinates(i,cand)
                    pix = Pixel(i,cand,pixcoords)
                    centerpt = pix.getCenterPoint(self.raster.srs)
                    print i, cand, pixcoords
                    import pdb; pdb.set_trace()
                    if centerpt.Within(geom):
                        self.addPixel(pix)
                        print hit
    
    # Unused, could increase speed if incorporated
    #def findLooseEnd(self, candidates):
    #    for i,row in candidate.iteritems():
    #        for col in row:
    #            self.lookForNeighbours(i,col)

            

    def updateMaskMatrix(self,ma):
        for (row,rowdict) in self.pixels.iteritems():
            for (col,pixel) in rowdict.iteritems():
                ma[row][col] = 1
        return ma
            
class Raster:
    def __init__(self, 
                 filename = "test5.asc",
                 ox = 222000,
                 oy = 7674000,
                 sx = 1000,
                 sy = -1000,
                 nx = 740,
                 ny = 1570,
                 srs = None):
        self.filename = filename
        self.ox = ox
        self.oy = oy
        self.sx = sx
        self.sy = sy
        self.nx = nx
        self.ny = ny
        if srs == None:
            srs = osr.SpatialReference()
            srs.SetWellKnownGeogCS("EPSG:3006")
        self.srs = srs

        self.ma = np.zeros( (ny,nx) )
        #self.pixels = Pixelset(self)
    
    def printFile(self):
        ofh = open(self.filename,'w')
        ofh.write("ncols        %d\n" % self.nx)
        ofh.write("nrows        %d\n" % self.ny)
        ofh.write("xllcorner    %d\n" % self.ox)
        if self.sy > 0:
            ofh.write("yllcorner    %d\n" % self.oy)
        else:
            yll = self.oy + self.ny*self.sy
            ofh.write("yllcorner    %d\n" % yll)
        ofh.write("cellsize     %d\n" % self.sx)
        for row in self.ma:
            for val in row:
                ofh.write(' %d' % val)
            ofh.write("\n")
        ofh.close()

    def updatePixelsShape(self, geom):
        pixels = Pixelset(self)
        prevPixelNo = self.pixelNumber(geom.GetPoint(0))
        currentPixel = pixels.getPixel(prevPixelNo)
        corners = list()
        for i in range(0, geom.GetPointCount()):
            pixelNo = self.pixelNumber(geom.GetPoint(i))

            if pixelNo != prevPixelNo:
                self.updatePixelsEdge(pixels, geom.GetPoint(i-1), geom.GetPoint(i))
                currentPixel = pixels.getPixel(pixelNo)
                prevPixelNo = pixelNo
            currentPixel.addCorner(geom.GetPoint(i))
        pixels.addEnclosedPixels(geom)
        pixels.updateMaskMatrix(self.ma)

    def updatePixelsEdge(self,pixels,geom1,geom2):
        pixelNo1 = self.pixelNumber(geom1)
        pixelNo2 = self.pixelNumber(geom2)
        pixel1 = pixels.getPixel(pixelNo1)
        #pixel2 = pixels.getPixel(pixelNo2)

        #Hur många skärningar i resp led letar vi efter?
        nox = pixelNo2[1] - pixelNo1[1]
        noy = pixelNo2[0] - pixelNo1[0]

        if nox != 0:
            sigx = nox/abs(nox)
            # Riktningskoeff för linjen
            kx = (geom1[1]-geom2[1])/(geom1[0]-geom2[0])
        if noy != 0:
            sigy = noy/abs(noy)
            # inv rikningskoeff för att undvika div med noll
            ky = (geom1[0]-geom2[0])/(geom1[1]-geom2[1])

        nox = abs(nox)
        noy = abs(noy)
        
        #Räkna från mittpunkten av pixel1
        (x0,y0) = pixel1.getCenterCoords()

        for i in range(nox):
            x = x0 + (i+0.5)*self.sx*sigx
            y = kx*(x-geom1[0]) + geom1[1]
            # Get row number from pixelNumber(). Dont trust col nr as it is on
            # an edge.
            (row,scrat) = self.pixelNumber((x,y))
            col1 = pixelNo1[1] + i*sigx
            col2 = pixelNo1[1] + (i+1)*sigx
            p1 = pixels.getPixel((row,col1))
            p1.addEdge((x,y),3+sigx)
            p2 = pixels.getPixel((row,col2))
            p2.addEdge((x,y),3-sigx)

        for i in range(noy):
            y = y0 + (i+0.5)*self.sy*sigy
            x = ky*(y-geom1[1]) + geom1[0]
            # Get col number from pixelNumber(). Dont trust row nr as it is on
            # an edge.
            (scrat,col) = self.pixelNumber((x,y))
            row1 = pixelNo1[0] + i*sigy
            row2 = pixelNo1[0] + (i+1)*sigy
            p1 = pixels.getPixel((row1,col))
            p1.addEdge((x,y),2+sigy)
            p2 = pixels.getPixel((row2,col))
            p2.addEdge((x,y),2-sigy)

    def pixelNumber(self, point):
        row = int(floor((point[1]-self.oy)/self.sy))
        col = int(floor((point[0]-self.ox)/self.sx))
        return row,col

    def calcPixelCoordinates(self, row, col):
        x0 = self.ox + col*self.sx
        y0 = self.oy + row*self.sy
        return [x0, y0, self.sx, self.sy]

    def setval(self, point, val=1):
        self.ma[point[0],point[1]] = val

class Pixel:
    def __init__(self,row,col,coords):
        self.row = row
        self.col = col
        self.cornerpoints = list()
        self.edgepoints = list()
        self.points = list()
        self.egdecross = 0
        self.x0 = coords[0]
        self.y0 = coords[1]
        self.sx = coords[2]
        self.sy = coords[3]

    def __eq__(self,other):
        if isinstance(other,Pixel):
            return self.row == other.row and self.col == other.col
        elif isinstance(other,tuple) or isinstance(other,list):
            return self.row == other[0] and self.col == other[1]
        return NotImplemented

    def __ne__(self,other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        return hash(str(self.row)+"."+str(self.col))

    def getCenterCoords(self):
        xc = self.x0 + self.sx/2
        yc = self.y0 + self.sy/2
        return xc,yc

    def getCenterPoint(self,srs):
        (x,y) = self.getCenterCoords()
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AssignSpatialReference(srs)
        pt.SetPoint(0, x, y)
        return pt

    def addCorner(self,coords):
        self.cornerpoints.append(len(self.points))
        self.points.append(coords)

    def addEdge(self,coords,edge):
        self.edgepoints.append(len(self.points))
        self.points.append(coords)

    def setEnclosedPartial(self,partial):
        self.enclosedArea = self.getArea() * partial

    def getSquarePolygon(self):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint_2D(x0,y0)
        ring.AddPoint_2D(x0+sx,y0)
        ring.AddPoint_2D(x0+sx,y0+sy)
        ring.AddPoint_2D(x0,y0+sy)
        ring.AddPoint_2D(x0,y0)
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        return poly

    def getEnclosedPolygon(self,geom):
        squarePoly = self.getSquarePolygon()
        return squarePoly.Intersection(geom)

    def getEnclosedArea(self,geom):
        return self.getEnclosedPolygon(geom).GetArea()
        
    def getArea(self):
        return abs(sx*sy)

def main(argv):
    if len(argv) < 1:
        print "usage: mask.py infile outfile"
        sys.exit(1)
	
    inShapeFile = sys.argv[1]
    outRasterFile = 'test'
    outSuffix = 'asc'
    #outRasterFile = sys.argv[2]
    #if outRasterFile[-4:] == ".shp":
    #    outRasterFile =  outRasterFile[:-4]
    maskFile(inShapeFile,outRasterFile)

def maskFile(inShapeFile,outRasterFile):	
    driver = ogr.GetDriverByName('ESRI Shapefile')
    inDataSource = driver.Open(inShapeFile,0)
    inLayer = inDataSource.GetLayer()
	
    # Remove output shapefile if it already exists
    #if os.path.exists(outRasterFile+'.'+outSuffix):
        #driver.DeleteDataSource(outRasterFile+"."+outSuffix)
	
    #outDataSource = driver.CreateDataSource(outRasterFile+"."+outSuffix)
    #srs = osr.SpatialReference()
    #srs.ImportFromEPSG(3006)
    maskArray = Raster(srs=inLayer.GetSpatialRef())

    for i in range(inLayer.GetFeatureCount()):
        inFeature = inLayer.GetFeature(i)
        geometries = inFeature.GetGeometryRef()
        for j in range(geometries.GetGeometryCount()):
            shapelim = geometries.GetGeometryRef(j)
            maskArray.updatePixelsShape(shapelim)

    #import pdb; pdb.set_trace()
    print "Printing result to file", maskArray.filename
    maskArray.printFile()

if __name__ == "__main__":
    main(sys.argv[1:])
