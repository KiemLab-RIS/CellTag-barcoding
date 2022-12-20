#
#  look through R1 and R1 of 10x fastq files for 
#    R1: [cell barcode][umi]
#    R2:[expBarcode gfp seq][expBarcode][anchor]
#
#  R2 will be reverse complement so:
#
#  ACGAGCTGTACAAGTAAACCGGT........GAA
#
#
# example fastq
#
# CGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTaaaccggtGATGCATCG
# GATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCAATCTCGGCATGGACGAGCTGTACAAGTaaaccggtTTTTTTTTTTGA
# GGATCACTCTCGGCATGGACGAGCTGTACAAGTaaaccggtATTAACTCgaattCGATGACAGGCGCAGCTTCCGAGGGATTTGAGATCC
# GGACGAGCTGTACAAGTaaaccggtAGTACAAGgaattCGATGACAGGCGCAGCTTCCGAGGGATTTGAGATCCAGACATGATAAGATAC
# GCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTaaaccggtGTCGTATGgaattCGATGACAGGCGC
#
import sys
import os
from collections import defaultdict
import argparse
#-----------------------------------------------------------------------------------------------
# main start
#-----------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Build script files to run gene editing pipeline')
parser.add_argument("-r1", "--read1", type=str, help="fastq read 1",required=True)
parser.add_argument("-r2", "--read2", type=str, help="fastq read 2",required=True)
parser.add_argument("-o", "--outFile", type=str, help="output fil,prefix",required=True)
args = parser.parse_args()
#
# 
#
file1 = args.read1
file2 = args.read2
out   = args.outFile
print(f'expBarcodeAn3.py {file1} {file2} {out}')
#
# make sure output directory exists
#
path = './'
if not os.path.isdir('filtered'):
  os.makedirs('filtered')
#--------------------------------------------------------------------------------
#
# checkQuality: how many bases are below Q
#
#--------------------------------------------------------------------------------
def checkQuality(qualSeq):
  #
  # count how many bases are below low quality
  #
  fail = 0
  i = 0
  for char in qualSeq:
    q = ord(char) - ord('!')
    if (q < 15):
      fail += 1
  return fail

#--------------------------------------------------------------------------------
#
# bioFindStart  - find sub-seq in sequence 
#
#-------------------------------------------------------------------------------


def bioFindStart(seq):
  gp1 = 'ACGAGCTGTACAAGTAAACCGGT'
  anchor = 'GAA'
  bc = None
  loc = seq.find(gp1)
  # is there an exact match
  if loc > -1:
    #
    # yes... does location allow for 8 barcode bases  (and thee anchor bases)?
    #
    if len(seq) - (loc + len(gp1)) >= (8 + len(anchor)):
      
      p1 = loc + len(gp1)
      p2 = p1 + 8
      p3 = p2 + 3
      anchorSeq = seq[p2:p3]
      expbc = seq[p1:p2]
      if anchorSeq == anchor:
        #print()
        #print(f'{seq}')
        #print(f'{" " * loc}{gp1}')
        #print(f'{" " * (loc + len(gp1))}{expbc}')
        return(expbc)
  return('')

#--------------------------------------------------------------------------------
# read main
#--------------------------------------------------------------------------------

q1 = (line.strip() for line in open(file1,'r'))
q2 = (line.strip() for line in open(file2,'r'))
#
# cell barcode dict...just how many times was each cell barcode is seen
#
cellDict = defaultdict(int)
#
# cbc is a dictionary where the key is the cell barcode.
# the value is a list of (expBarcode,umi)
#
cbc = {}
#
totalExpBarcodes = 0
good = 0
bad = 0
index = 0

while True:
  index += 1
  try:
    id1  = next(q1)
    seq1 = next(q1)
    next(q1)
    qual1 = next(q1)

    id2  = next(q2)
    seq2 = next(q2)
    next(q2)
    qual2 = next(q2)

    #print(qual)
    if checkQuality(qual1) > 4 or checkQuality(qual2) > 4:
      bad += 1
    else:
      good+= 1

      cellBC = seq1[0:16]
      umi    = seq1[16:]

      cellDict[cellBC] += 1
      #
      # is there an expression barcode?
      #
      expBC = bioFindStart(seq2)
      if expBC != '':
        totalExpBarcodes += 1
        # add this expBarcode to cell
        # not currently using UMI
        try:
          tagList = cbc[cellBC]
          tagList.append((expBC,umi))
        except KeyError:
          cbc[cellBC] = [(expBC,umi)]
    
  except StopIteration:
    break

  if index % 1000000 == 0:
    sd = sorted(cellDict.items(),key=lambda x:x[1],reverse=True)
    ssd = sorted(cbc.items(),key=lambda x:len(x[1]),reverse=True) 
    print(f'{index}   goodQ = {good}   badQ = {bad} most frequent cell bc = {sd[0][0]} count = {sd[0][1]}    total cell bc = {len(cellDict)}')
    print(f'   length of cell barcode dictionary = {len(cbc)}')
    if len(ssd) > 0:
      k = ssd[0][0]
      l = ssd[0][1]
      print(f'cell bc with most tags = {k} length = {len(l)}')

  #if index > 1000000:
  #  break


sd = sorted(cellDict.items(),key=lambda x:x[1],reverse=True)
ssd = sorted(cbc.items(),key=lambda x:len(x[1]),reverse=True) 

out1 = 'filtered/' + out + '_cellDict.log'
out2 = 'filtered/' + out + '_expDict.tsv'
out3 = 'filtered/' + out + '_dist.log'
#
# cell barcode counts
#
with open (out1,'w') as fh:
  fh.write(f'#CellBarcode Count     : total exp barcodes = {totalExpBarcodes}\n')
  for k,v in sd:
    fh.write(f'{k} {v}\n')
#
# cell barcode with expBarcode list
#
with open (out2,'w') as fh:
  fh.write('#cell barcode \t exp tag \t umi\n')
  for k,v in ssd:
    for tag,umi in v:
      fh.write(f'{k}')
      fh.write(f'\t{tag}\t{umi}\n')
#
# cell exp barcode list
#
with open (out3,'w') as fh:
  fh.write(f'#  cell, total unique tags\n')
  for k,v in ssd:
    tags = defaultdict(int)
    total = 0
    for tag,umi in v:
      tags[tag] += 1
    fh.write(f'{len(tags)}\t{total}\n')


