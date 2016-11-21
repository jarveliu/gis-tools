#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Split an excelfile in csv files of equal size
#
# The input excelfile may contain several sheets but will loop through those in sequence
# and the rows will be written to files with names based on outfilebase and a sequential
# number. First row in each sheet is supposed to be a header. If the first row in a 
# sheet is equal to the previous rows are appended to the open csv file. If they are not
# a new csv file is opened beginning with the header of the new sheet.
# 
# Johan Arvelius, SMHI, 2016-11-21

import argparse
import codecs
import openpyxl

parser = argparse.ArgumentParser()
parser.add_argument('--infile', help="input xlsx file")
parser.add_argument('--outfilebase', help="output base filename")
parser.add_argument('--rows_per_file', help="number of rows per output file", 
                    default = 2500)
parser.add_argument('--add_to_rows',help="extra text to put after last comma on each row",
                    default = None)

su = parser.parse_args()

#import pdb; pdb.set_trace()

wb = openpyxl.load_workbook(filename=su.infile, read_only=True)

class LimitedFileLimitError(Exception):
    pass

class Writer:
    def __init__(self,su):
        self.filenumber = 0

    def writerow(self,row):
        try:
            self.outfile.writerow(row)
        except LimitedFileLimitError:
            self.newfile(self.header)
            self.outfile.writerow(row)

    def newfile(self,header):
        self.header = header
        if hasattr(self,'outfile'):
            self.outfile.close()
        self.filenumber += 1
        self.outfile = LimitedFile(self.filenumber,su)
        self.outfile.writerow(header)
        

class LimitedFile:
    def __init__(self, number, su):
        self.rowcounter = 0
        self.fh = codecs.open(su.outfilebase+unicode(number)+'.csv','w','utf8')
        print "opening file "+su.outfilebase+unicode(number)+'.csv'

    def close(self):
        self.fh.close()

    def writerow(self,row):
        if self.rowcounter < su.rows_per_file:
            for cell in row:
                self.fh.write(cell.value)
                self.fh.write(",")
            if su.add_to_rows:
                self.fh.write(su.add_to_rows)
            self.fh.write("\n")
            self.rowcounter += 1
        else:
            raise LimitedFileLimitError("maximum number of rows reached")

def comprows(r1,r2):
    citer = iter(r2)
    if len(r1) != len(r2):
        return False
    for cell in r1:
        if cell.value != next(citer).value:
            return False
    return True

wtr = Writer(su)
header = list()
for ws in wb.worksheets:
    rowiter = iter(ws.rows)
    tmpheader = next(rowiter)
    if len(header)==0 or not comprows(header,tmpheader):
        header = tmpheader
        wtr.newfile(header)
        
    for row in rowiter:
        if row[0].value:
            wtr.writerow(row)
