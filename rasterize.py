#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import ogr, os, sys, osr
import numpy as np
from math import floor
from copy import deepcopy

rasterize_parser = argparse.ArgumentParser(add_help=False)
rasterize_parser.add_argument('infile')

class Pixelset:
    """A set of pixels in a raster.

    Attributes:
    pixels -- a dictionary of dictionaries, with integer keys containing row
              and column numbers, containing Pixel objects.
    raster -- a Raster object the Pixelset belongs to
    feature -- a shapefile feature object that used to instansiate the 
               Pixelset
    boundingbox -- bb of the Pixelset as 4 int array, [rowmin, rowmax, 
                   colmin, colmax], kept to improve speed in comparisons.
    """
    
    def __init__(self,raster,feature=None):
        """With a pointer to a raster it belongs to.

        positional argument:
        raster -- a Raster object the Pixelset belongs to

        keyword argument:
        feature -- a shapefile feature object that should be used to 
                   instansiate the Pixelset
        """
        self.pixels = dict()
        self.raster = raster
        #                   rowmin, rowmax, colmin, colmax
        self.boundingbox = [raster.ny, -1, raster.nx, -1]
        self.feature = feature
        if feature is not None:
            self.updatePixelsShape(feature)

    def getPixels(self):
        """Generator for all Pixels in the set"""
        for rownr, row in self.pixels.iteritems():
            for colnr, pix in row.iteritems():
                yield pix

    def getPixel(self,rowcol):
        """Return a Pixel object for given coordinates.

        If a Pixel of the given row and column values exists it is returned,
        if not a new empty Pixel object is added to set and returned.

        positional argument:
        rowcol -- two integer list with row and column numbers.
        """
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
        """Add a Pixel object to the set.

        positional arguments:
        pixel -- Pixel object to add.
        row -- (int) row number
        col -- (int) col number

        exception:
        Exception -- if a pixel with those row and column numbers already 
                     exists in set.
        """
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
        """Return a dict of lists with row and column numbers of pixels.

        A lightweight representation of the object."""
        retstruct = dict()
        for (rownr,row) in self.pixels.iteritems():
            retstruct[rownr] = set(row.keys())
        return retstruct

    def updatePixelsShape(self, feature):
        """Update the set with pixels corresponding to feature.

        positional argument:
        feature -- ogr feature object"""
        geom = feature.GetGeometryRef()
        self.updatePixelsRecursion(geom)
        self.addEnclosedPixels(geom)

    def updatePixelsRecursion(self,geom):
        """Recursion function for updatePixelsShape.

        positional argument:
        geom -- ogr geometry object, may be a ring or a geometry object 
                containing one or more rings at any level of nesting.
        """
        for i in range(geom.GetGeometryCount()):
            self.updatePixelsRecursion(geom.GetGeometryRef(i))
        if geom.GetGeometryCount() == 0:
            self.updatePixelsRing(geom)

    def updatePixelsRing(self, ring):
        """Update the set with pixels corresponding to ogr ring.

        positional argument:
        ring -- ogr ring object.
        """
        prevPixelNo = self.raster.pixelNumber(ring.GetPoint(0))
        currentPixel = self.getPixel(prevPixelNo)
        for j in range(0, ring.GetPointCount()):
            pixelNo = self.raster.pixelNumber(ring.GetPoint(j))

            if pixelNo != prevPixelNo:
                self.updatePixelsEdge(ring.GetPoint(j-1), ring.GetPoint(j))
                currentPixel = self.getPixel(pixelNo)
                prevPixelNo = pixelNo
            currentPixel.addCorner(ring.GetPoint(j))

    def updatePixelsEdge(self,geom1,geom2):
        """Update the set with pixels on a linte between two points.

        positional arguments:
        geom1 and geom2 -- coordinates for the two endpoints of the line in 
                           the coordinate system used in Raster

        note:
        Calculates the exact positions of all pixel edge crossnings. This is
        at the moment waisted effort as the information is not used in the 
        area coverage calculation.
        """
        pixelNo1 = self.raster.pixelNumber(geom1)
        pixelNo2 = self.raster.pixelNumber(geom2)
        pixel1 = self.getPixel(pixelNo1)

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
        """Update the set with the pixels fully enclosed in geometry.

        To be run after all pixels containing to corners and edges of 
        the geometry is already identified and added.

        positional argument:
        geom -- ogr geometry object
        """
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

    #def updateMaskMatrix(self,ma):
    #    for (row,rowdict) in self.pixels.iteritems():
    #        for (col,pixel) in rowdict.iteritems():
    #            ma[row][col] = id(self)
    #    return ma

class BUMMask:
    """A mask type used in BUM

    A special mask type where a set of shape areas are read from a file,
    all pixel in contact with a shape from the input are marked, and all 
    sets of pixels masked that is in direct connection to each other are
    merged together in areas (Area objects). Around those clusters are a 
    bufferzone of one pixel. In the end a subset of such clusters that does
    not overlap or get in contact to each other are written to file, a "BUM
    mask". 

    Attributes:
    su -- setup object with settings from command line attributes.
    raster -- a Raster object.
    rows -- number of rows in raster
    cols -- number of columns in raster.
    areas -- dict of Area objects.
    masks -- dict of lists of Area objects that form a BUM mask.
    addedIds -- list of Id numbers for the Pixelset objects added.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--outfilebasename', default="mask")
    def __init__(self, su, raster):
        """Initialize with setup object and Raster object"""
        self.su = su
        self.raster = raster
        self.rows = raster.ny
        self.cols = raster.nx
        self.areas = dict()
        self.masks = dict()
        self.addedIds = list()

    def addPixelSet(self, pixelset):
        """Add a Pixelset object"""
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
        """Store a buffer zone around each area."""
        for areaId,area in self.areas.iteritems():
            area.padSelf()

    def createMasks(self):
        """Create the BUM masks (Attribute masks)"""
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
        """Print all BUM masks to asc files."""
        for maskNo in self.masks:
            self.printmask(maskNo)
                    
    def printmask(self,maskNo):
        """Print a BUM mask to an asc file.

        positional argument:
        maskNo -- zerobased number of mask.
        """
        ma = np.zeros( (self.raster.ny,self.raster.nx) )
        for areaId in self.masks[maskNo]:
            self.areas[areaId].printToMatrix(ma)
            self.areas[areaId].printPadToMatrix(ma)
        filename = self.su.outfilebasename+str(maskNo)+".asc"
        self.raster.printFile(ma,filename)
            
class Area:
    """A set of pixels containing all pixels in contact with eachother.

    The class for structures building up a BUM mask. Sets of pixels read
    from pixelsets but pixelsets overlaping eachother or containing pixels
    neighbouring eachother are merged together. 

    Attributes:
    mask -- a BUMMask object
    boundingbox -- similar to Pixelset boundingbox
    pixelstruct -- dict of lists of row and column numbers of contained pix
    idnr -- ID number taken from the pixelset the area is built from
    """

    def __init__(self,mask,origin):
        """Initialized with a BUMMask and a Pixelstruct object."""
        self.mask = mask
        self.boundingbox = origin.boundingbox
        self.pixelstruct = origin.getPixelStruct()
        self.idnr = id(origin)

    def getId(self):
        """Return unique ID number for the object."""
        return self.idnr

    def padset(self,theset):
        """Take a set of integers return a set with also adjusent numbers"""
        retset = deepcopy(theset)
        for nr in theset:
            retset.add(nr-1)
            retset.add(nr+1)
        return retset

    def padSelf(self):
        """Store a self.pad pad with pixels neigbouring pixelstruct."""
        self.pad = self.padStruct(self.pixelstruct)

    def padStruct(self,struct):
        """Calculate a pad with pixels neigbouring a structure.

        positional argument:
        struct -- a dict of sets of integers where dict keys are row numbers
                  and sets contains column numbers for pixels.
        """
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
        """Return a structure that is a join of pixelstruct and pad"""
        retstruct = self.padStruct(self.pixelstruct)
        for rownr,row in self.pixelstruct.iteritems():
            retstruct[rownr] |= row
        return retstruct

    def getDoublePad(self):
        """Return a structure that is padded twice."""
        struct = self.getOutPaddedStruct()
        return self.padStruct(struct)

    def overlap(self,other):
        """Return true if other Area overlap with self else false.

        positional argument:
        other -- another Area object
        """
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
        """Return true if other Area overlap with selfs pad else false.

        positional argument:
        other -- another Area object
        """
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
        """Return true if other Areas pad overlap with selfs doublepad.

        positional argument:
        other -- another Area object
        """
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
        """Merge other Area into self."""
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
        """Print ones to pixels covered by the area to matrix"""
        for rownr,row in self.pixelstruct.iteritems():
            for col in row:
                ma[rownr][col] = 1

    def printPadToMatrix(self,ma):
        """Print twos to pixels covered by the 1px pad around area to matrix"""
        for rownr,row in self.pad.iteritems():
            for col in row:
                ma[rownr][col] = 2


class Raster:
    """A printable raster geometry and collection of Pixelset objects.

    Class variable:
    argparse.ArgumentParser, may be used as parent to script parsers.

    Attributes:
    srs -- osr spatial reference object
    ox -- origo x lower left coordinate
    oy -- origo y lower left coordinate
    sx -- x direction step between pixels
    sy -- y direction step between pixels
    nx -- number of columns in raster
    ny -- number of rows in raster
    pixelsets -- dict containing Pixelset objects, keys are Pixelset IDs.
    ma -- numpy array containing raster numerical values.
    """
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--rastergeometry', nargs=6,
                        default=[222000,7674000,1000,-1000,740,1570],
                        help='origin_x origin_y step_x step_y No_of_pixels_x No_of_pixels_y')
    parser.add_argument('outfile')
    def __init__(self, 
                 su = None,
                 srs = None):
        """Constructor
        
        positional arguments:
        su -- argparse.Namespace object with settings. Should contain 
              rastergeometry and outfile attributes.
        srs -- osr spatial reference object.
        """
        
        self.ox = su.rastergeometry[0]
        self.oy = su.rastergeometry[1] 
        self.sx = su.rastergeometry[2] 
        self.sy = su.rastergeometry[3] 
        self.nx = su.rastergeometry[4] 
        self.ny = su.rastergeometry[5] 

        if hasattr(su,'filename') \
           or isinstance(su,argparse.Namespace) and su.__contains__('outfile'):
            self.filename = su.outfile
        if srs == None:
            srs = osr.SpatialReference()
            srs.SetWellKnownGeogCS("EPSG:3006")
        self.srs = srs

        self.pixelsets = dict()

    def updateMatrix(self,field):
        """create the ma attribute with values from pixelsets

        positional argument:
        field -- name of the field in shapefile go get values from.
        """
        self.ma = np.zeros( (self.ny,self.nx) )
        for setId,ps in self.pixelsets.iteritems():
            geom = ps.feature.GetGeometryRef()
            value = ps.feature.GetField(field)
            area = geom.GetArea()
            density = value/area
            for pixel in ps.getPixels():
                self.ma[pixel.row][pixel.col] \
                    += pixel.getEnclosedArea(geom)*density
    
    def printFile(self, matrix=None, filename=None):
        """print matrix values to an ascii raster file

        keyword arguments:
        matrix -- matrix to print, default is self.ma
        filename -- name of output file, default is self.filename
        """
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
        """create a new Pixelset and add to pixelsets attribute

        positional argument:
        feature -- ogr feature object to add 
        """
        pixels = Pixelset(self, feature)
        self.pixelsets[id(pixels)] = pixels
        #pixels.updateMaskMatrix(self.ma)

    def pixelNumber(self, point):
        """return the row and column number for pixel covering point
        
        positional argument:
        point -- two number list with x and y coordinates for point
        """
        row = int(floor((point[1]-self.oy)/self.sy))
        col = int(floor((point[0]-self.ox)/self.sx))
        return row,col

    def calcPixelCoordinates(self, row, col):
        """return pixel corner coordinates and size for pixel in row, col

        positional arguments:
        row -- row number for pixel.
        col -- col number for pixel.

        return values:
        4 numbers list: [x0, y0, delta_x, delta_y]
        """
        x0 = self.ox + col*self.sx
        y0 = self.oy + row*self.sy
        return [x0, y0, self.sx, self.sy]

    def setVal(self, point, val=1):
        """Set pixel value for a pixel

        positional argument:
        point -- 2 numbers list: row, col

        keyword argument:
        val -- value, defaults to 1
        """
        self.ma[point[0]][point[1]] = val

class Pixel:
    """Single pixel object.

    Attributes:
    row -- row number in raster
    col -- column number in raster
    cornerpoints -- list of ogr geometry object corner points in pixel
    edgepoints -- list of points on pixel edge with crossings by lines
                  connecting ogr geometry object corners.
    points -- list of ogr geometry object points in pixel 
    edgecross --
    x0 -- x-coordinate of the corner closest to the raster origo
    y0 -- y-coordinate of the corner closest to the raster origo
    sx -- "delta x" x0+sx gives the x coordinate of the corner furtest from
          the raster origo.
    sy -- similar in y direction.
    """
    
    def __init__(self,row,col,coords):
        """constructor

        Attributes:
        row -- row number in raster
        col -- column number in raster
        coords -- 4 number list with coordinates [x0, y0, sx, sy]
        """
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
        """return center of pixel as 2 numbers, x,y"""
        xc = self.x0 + self.sx/2
        yc = self.y0 + self.sy/2
        return xc,yc

    def getCenterPoint(self,srs):
        """return center of pixel as ogr Geometry point object

        positional argument:
        srs -- osr spatial reference object for return point
        """
        (x,y) = self.getCenterCoords()
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AssignSpatialReference(srs)
        pt.SetPoint(0, x, y)
        return pt

    def addCorner(self,coords):
        """add a corner point to the points attribute

        positional argument:
        coords -- coordinates
        """
        self.cornerpoints.append(len(self.points))
        self.points.append(coords)

    def addEdge(self,coords,edge):
        self.edgepoints.append(len(self.points))
        self.points.append(coords)

    def setEnclosedPartial(self,partial):
        self.enclosedArea = self.getArea() * partial

    def getSquarePolygon(self):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint_2D(self.x0,self.y0)
        ring.AddPoint_2D(self.x0+self.sx,self.y0)
        ring.AddPoint_2D(self.x0+self.sx,self.y0+self.sy)
        ring.AddPoint_2D(self.x0,self.y0+self.sy)
        ring.AddPoint_2D(self.x0,self.y0)
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        return poly

    def getEnclosedPolygon(self,geom):
        squarePoly = self.getSquarePolygon()
        return squarePoly.Intersection(geom)

    def getEnclosedArea(self,geom):
        encpoly = self.getEnclosedPolygon(geom)
        ##Here the warning comes about GetArea called on non-surface
        retval = encpoly.GetArea()
        return retval
        
    def getArea(self):
        return abs(sx*sy)

def main():
    parser = argparse.ArgumentParser(parents=[rasterize_parser, Raster.parser])
    parser.add_argument('field')
    su = parser.parse_args()

    driver = ogr.GetDriverByName('ESRI Shapefile')
    inDataSource = driver.Open(su.infile,0)
    inLayer = inDataSource.GetLayer()
    myRaster = Raster(su=su, srs=inLayer.GetSpatialRef())
    for i in range(inLayer.GetFeatureCount()):
        #print i
        #if i == 144:
        #    import pdb; pdb.set_trace()
        inFeature = inLayer.GetFeature(i)
        myRaster.updatePixelsShape(inFeature)
    myRaster.updateMatrix(su.field)
    myRaster.printFile()

if __name__ == "__main__":
    main()
