#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ogr, os, sys, osr
import numpy as np
from math import floor
from copy import deepcopy

class Pixelset:
    def __init__(self,raster,feature=None):
        self.pixels = dict()
        self.raster = raster
        #                   rowmin, rowmax, colmin, colmax
        self.boundingbox = [raster.ny, -1, raster.nx, -1]

        if feature is not None:
            self.updatePixelsShape(feature)

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
            raise Exception

        self.pixels[row][col] = pixel

        if row < self.boundingbox[0]: 
            self.boundingbox[0] = row
        if row > self.boundingbox[1]: 
            self.boundingbox[1] = row
        if col < self.boundingbox[2]:
            self.boundingbox[2] = col
        if col > self.boundingbox[3]:
            self.boundingbox[3] = col

    def getPixelStruct(self):
        retstruct = dict()
        for (rownr,row) in self.pixels.iteritems():
            retstruct[rownr] = set(row.keys())
        return retstruct

    def updatePixelsShape(self, feature):
        geom = feature.GetGeometryRef()
        for i in range(geom.GetGeometryCount()):
            ring = geom.GetGeometryRef(i)
            prevPixelNo = self.raster.pixelNumber(ring.GetPoint(0))
            currentPixel = self.getPixel(prevPixelNo)
            for j in range(0, ring.GetPointCount()):
                pixelNo = self.raster.pixelNumber(ring.GetPoint(j))

                if pixelNo != prevPixelNo:
                    self.updatePixelsEdge(ring.GetPoint(j-1), ring.GetPoint(j))
                    currentPixel = self.getPixel(pixelNo)
                    prevPixelNo = pixelNo
                currentPixel.addCorner(ring.GetPoint(j))
        self.addEnclosedPixels(geom)
        #self.updateMaskMatrix(self.ma)

    def updatePixelsEdge(self,geom1,geom2):
        pixelNo1 = self.raster.pixelNumber(geom1)
        pixelNo2 = self.raster.pixelNumber(geom2)
        pixel1 = self.getPixel(pixelNo1)
        #pixel2 = self.getPixel(pixelNo2)

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
            x = x0 + (i+0.5)*self.raster.sx*sigx
            y = kx*(x-geom1[0]) + geom1[1]
            # Get row number from pixelNumber(). Dont trust col nr as it is on
            # an edge.
            (row,scrat) = self.raster.pixelNumber((x,y))
            col1 = pixelNo1[1] + i*sigx
            col2 = pixelNo1[1] + (i+1)*sigx
            p1 = self.getPixel((row,col1))
            p1.addEdge((x,y),3+sigx)
            p2 = self.getPixel((row,col2))
            p2.addEdge((x,y),3-sigx)

        for i in range(noy):
            y = y0 + (i+0.5)*self.raster.sy*sigy
            x = ky*(y-geom1[1]) + geom1[0]
            # Get col number from pixelNumber(). Dont trust row nr as it is on
            # an edge.
            (scrat,col) = self.raster.pixelNumber((x,y))
            row1 = pixelNo1[0] + i*sigy
            row2 = pixelNo1[0] + (i+1)*sigy
            p1 = self.getPixel((row1,col))
            p1.addEdge((x,y),2+sigy)
            p2 = self.getPixel((row2,col))
            p2.addEdge((x,y),2-sigy)

    def addEnclosedPixels(self, geom):
        for i in range(self.boundingbox[0]+1, self.boundingbox[1]):
            if i in self.pixels.keys():
                sortedkeys = sorted(self.pixels[i].keys())
                envelope = range(sortedkeys[0], sortedkeys[-1])
                cands = filter(lambda x:x not in sortedkeys,envelope)
                for cand in cands:
                    pixcoords = self.raster.calcPixelCoordinates(i,cand)
                    pix = Pixel(i,cand,pixcoords)
                    centerpt = pix.getCenterPoint(self.raster.srs)
                    if centerpt.Within(geom):
                        self.addPixel(pix,i,cand)
    
    # Unused, could increase speed if incorporated
    #def findLooseEnd(self, candidates):
    #    for i,row in candidate.iteritems():
    #        for col in row:
    #            self.lookForNeighbours(i,col)

            

    def updateMaskMatrix(self,ma):
        for (row,rowdict) in self.pixels.iteritems():
            for (col,pixel) in rowdict.iteritems():
                ma[row][col] = id(self)
        return ma

class Mask:
    def __init__(self, raster):
        self.raster = raster
        self.rows = raster.ny
        self.cols = raster.nx
        self.areas = dict()
        self.masks = dict()
        self.addedIds = list()

    def addPixelSet(self, pixelset):
        self.addedIds.append(id(pixelset))
        #self.sets[id(pixelset)] = pixelset
        #self.sets.append(pixelset)
        newarea = Area(self, pixelset)
        adsorbed = list()
        for areaId,area in self.areas.iteritems():
            if newarea.padOverlap(area):
                #import pdb; pdb.set_trace()
                newarea.merge(area)
                adsorbed.append(areaId)
        self.areas[newarea.getId()] = newarea
        for areaId in adsorbed:
            del self.areas[areaId]

    def padAreas(self):
        for areaId,area in self.areas.iteritems():
            area.padSelf()

    def createMasks(self):
        self.padAreas()
        for fitarea in self.areas:
            fit = False
            for nr in sorted(self.masks, reverse=True):
                fit = True
                for fittedarea in self.masks[nr]:
                    if self.areas[fitarea].doublePadOverlap(self.areas[fittedarea]):
                        fit = False
                        break
                if fit == True:
                    self.masks[nr].append(fitarea)
                    break
            if fit == False:
                self.masks[len(self.masks)] = list()
                self.masks[len(self.masks)-1].append(fitarea)
                
    def printmasks(self):
        for maskNo in self.masks:
            self.printmask(maskNo)
                    
    def printmask(self,maskNo):
        ma = np.zeros( (self.raster.ny,self.raster.nx) )
        for areaId in self.masks[maskNo]:
            self.areas[areaId].printToMatrix(ma)
            self.areas[areaId].printPadToMatrix(ma)
        filename = "mask"+str(maskNo)+".asc"
        self.raster.printFile(ma,filename)
            
class Area:
    def __init__(self,mask,origin):
        self.mask = mask
        self.boundingbox = origin.boundingbox
        self.pixelstruct = origin.getPixelStruct()
        self.idnr = id(origin)

    def getId(self):
        return self.idnr

    def padset(self,theset):
        retset = deepcopy(theset)
        for nr in theset:
            retset.add(nr-1)
            retset.add(nr+1)
        return retset

    def padSelf(self):
        self.pad = self.padStruct(self.pixelstruct)

    def padStruct(self,struct):
        #import pdb; pdb.set_trace()
        pad = dict()
        sortedrows = sorted(struct.keys())
        prevpadset = self.padset(struct[sortedrows[0]])
        pad[sortedrows[0]-1] = deepcopy(prevpadset)
        for rownr in sortedrows:
            pad[rownr] = deepcopy(prevpadset)
            prevpadset = self.padset(struct[rownr])
            pad[rownr-1] |= prevpadset
            pad[rownr] |= prevpadset
        pad[sortedrows[-1]+1] = deepcopy(prevpadset)

        for rownr in sortedrows:
            pad[rownr] -= struct[rownr]
        return pad

    def getOutPaddedStruct(self):
        retstruct = self.padStruct(self.pixelstruct)
        for rownr,row in self.pixelstruct.iteritems():
            retstruct[rownr] |= row
        return retstruct

    def getDoublePad(self):
        struct = self.getOutPaddedStruct()
        return self.padStruct(struct)

    #def lists_overlap(a, b):
    #    for i in a:
    #        if i in b:
    #            return True
    #    return False

    def overlap(self,other):
        if self.boundingbox[0] > other.boundingbox[1] or \
           self.boundingbox[1] < other.boundingbox[0] or \
           self.boundingbox[2] > other.boundingbox[3] or \
           self.boundingbox[3] < other.boundingbox[2]:
            return False
        for (rownr,row) in self.pixelstruct.iteritems():
            if rownr in other.pixelstruct.keys() and \
               bool(row & other.pixelstruct[rownr]):
                return True
        return False
    
    def padOverlap(self,other):
        if self.boundingbox[0]-1 > other.boundingbox[1] or \
           self.boundingbox[1]+1 < other.boundingbox[0] or \
           self.boundingbox[2]-1 > other.boundingbox[3] or \
           self.boundingbox[3]+1 < other.boundingbox[2]:
            return False
        padded = self.getOutPaddedStruct()
        for (rownr,row) in padded.iteritems():
            if rownr in other.pixelstruct.keys() and \
               bool(row & other.pixelstruct[rownr]):
                return True
        return False

    def doublePadOverlap(self,other):
        if self.boundingbox[0]-2 > other.boundingbox[1]+1 or \
           self.boundingbox[1]+2 < other.boundingbox[0]-1 or \
           self.boundingbox[2]-2 > other.boundingbox[3]+1 or \
           self.boundingbox[3]+2 < other.boundingbox[2]-1:
            return False
        outerpad = self.getDoublePad()
        for (rownr,row) in outerpad.iteritems():
            if rownr in other.pad.keys() and \
               bool(row & other.pad[rownr]):
                return True
        return False
    
    def merge(self,other):
        self.boundingbox[0] = min(self.boundingbox[0], other.boundingbox[0])
        self.boundingbox[1] = max(self.boundingbox[1], other.boundingbox[1])
        self.boundingbox[2] = min(self.boundingbox[2], other.boundingbox[2])
        self.boundingbox[3] = max(self.boundingbox[3], other.boundingbox[3])
        for (rownr,row) in other.pixelstruct.iteritems():
            if rownr in self.pixelstruct.keys():
                self.pixelstruct[rownr] |= row
            else:
                self.pixelstruct[rownr] = row
    
    def printToMatrix(self,ma):
        for rownr,row in self.pixelstruct.iteritems():
            for col in row:
                ma[rownr][col] = 1

    def printPadToMatrix(self,ma):
        for rownr,row in self.pad.iteritems():
            for col in row:
                ma[rownr][col] = 2


class Raster:
    def __init__(self, 
                 filename = "test7.asc",
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
        self.mask = Mask(self)
        self.pixelsets = dict()
    
    def printFile(self, matrix=None, filename=None):
        if matrix == None:
            matrix = self.ma
        if filename == None:
            filename = self.filename
        ofh = open(filename,'w')
        ofh.write("ncols        %d\n" % self.nx)
        ofh.write("nrows        %d\n" % self.ny)
        ofh.write("xllcorner    %d\n" % self.ox)
        if self.sy > 0:
            ofh.write("yllcorner    %d\n" % self.oy)
        else:
            yll = self.oy + self.ny*self.sy
            ofh.write("yllcorner    %d\n" % yll)
        ofh.write("cellsize     %d\n" % self.sx)
        for row in matrix:
            for val in row:
                ofh.write(' %d' % val)
            ofh.write("\n")
        ofh.close()
    
    def updatePixelsShape(self, feature):
        pixels = Pixelset(self, feature)
        self.pixelsets[id(pixels)] = pixels
        #pixels.updateMaskMatrix(self.ma)

    def pixelNumber(self, point):
        row = int(floor((point[1]-self.oy)/self.sy))
        col = int(floor((point[0]-self.ox)/self.sx))
        return row,col

    def calcPixelCoordinates(self, row, col):
        x0 = self.ox + col*self.sx
        y0 = self.oy + row*self.sy
        return [x0, y0, self.sx, self.sy]

    def setVal(self, point, val=1):
        self.ma[point[0]][point[1]] = val

    def setMaskVal(self, point, val=1):
        self.Mask

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
    myRaster = Raster(srs=inLayer.GetSpatialRef())

    for i in range(inLayer.GetFeatureCount()):
        inFeature = inLayer.GetFeature(i)
        geom = inFeature.GetGeometryRef()
        #maskArray.updatePixelsShape(geom)
        myRaster.updatePixelsShape(inFeature)
        #for j in range(geom.GetGeometryCount()):
        #    shapelim = geom.GetGeometryRef(j)
        #    maskArray.updatePixelsShape(shapelim)

    for setid,pixset in myRaster.pixelsets.iteritems():
        myRaster.mask.addPixelSet(pixset)
    myRaster.mask.createMasks()
    myRaster.mask.printmasks()
    ma = np.zeros( (myRaster.ny, myRaster.nx) )
    #for setid,pixset in myRaster.pixelsets.iteritems():
    #    pixset.updateMaskMatrix(ma)
    #    myRaster.printFile(ma,"oxelpixel.asc")

    #import pdb; pdb.set_trace()
    #print "Printing result to file", maskArray.filename
    #maskArray.printFile()

if __name__ == "__main__":
    main(sys.argv[1:])
