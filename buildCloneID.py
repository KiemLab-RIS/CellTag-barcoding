#
# buildCloneID [inputFiles] --o outfile
#
# read summary files from fastq of cell expression barcodes
#
# format of input files is tab separateed: [cellBarcode] ([tag umi] * N)
#
#
# jaccard merge does not sort and renumber output, maybe change back to do
# this?
#
from subprocess import call
from sys import argv
import csv
import os
import time
import argparse
from collections import defaultdict
#
# command strings
#

#-----------------------------------------------------------------------------------------------
# main start
#-----------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Build script files to run gene editing pipeline')
parser.add_argument('inputFiles',type=str,nargs='*',help='control file(s) specifying input sequencing files and editing detials')
parser.add_argument("-o", "--outFile", type=str, help="output file prefix",required=True)
parser.add_argument("-b", "--barcode", type=str, nargs='*',help="10X cell barcode file",required=True)
parser.add_argument("-p", "--prefix", type=str, nargs='*',help="10X cell barcode prefix")
parser.add_argument("-w", "--whitelist", type=str, help="Expressed Barcode whitelist",required=True)
args = parser.parse_args()

files = args.inputFiles
out = args.outFile
barfiles = args.barcode
prefix = args.prefix
whiteFile = args.whitelist

print(files,out,barfiles,prefix,whiteFile)
#
# add placeholder if no prefix given
#
if (prefix is None) or len(prefix) == 0:
  prefix = []
  for index in range(0,len(barfiles)):
    prefix.append('')



if len(barfiles) != len(prefix):
  print('error, barcode files and prefixes must be same size')
  quit()
#
# globals tag whitelist and 10x valid barcode list
#
whitelist = []
crBarcodeList = []
#------------------------------------------------------------------------------------
#
# read whitelist
#
#------------------------------------------------------------------------------------
def readWhiteList(whiteF):
    white = []
    print(f'read WhiteList{white}')
    wl = (line.strip().split(',') for line in open(whiteF,'r'))
    next(wl)  # skip header
    while True:
        try:
            white.append(next(wl)[0])
        except StopIteration:
            break
    print(f'whitelist complete: length = {len(white)}')
    return(white)
#------------------------------------------------------------------------------------
#
#
# compare to whitelist
#
#
#------------------------------------------------------------------------------------
def compToWhite(tag,whitelist):
  minError = 8
  minIndex = -1
  for wi,w in enumerate(whitelist):
    error = 0
    for i in range(0,8):
      if tag[i] != w[i]:
        error += 1
    if error < minError:
      minError = error
      minIndex = wi
  return((minError,minIndex))
#------------------------------------------------------------------------------------
#
# read barcode list
# return a list of all barcodes
# and a list[list] of barcodes with prefix attached
#
#------------------------------------------------------------------------------------
def readBarcodeList(barFiles,prefix):
    bl = []
    bl_pre_list = []
    for barFile,px in zip(barFiles,prefix):
      bl_pre = []
      bl_pre_list.append(bl_pre)
      print(f'read 10x barcode file {barFile}')
      bcl = (line.strip() for line in open(barFile,'r'))
      while True:
        try:
            cellbc = next(bcl)[:-2]
            if len(cellbc) != 16:
              print(f'format error in barcode {cellbc} from line:')
              print(bcl)
              quit()
            bl.append(cellbc)
            bl_pre.append(px + cellbc)
        except StopIteration:
            break

      print(f'finished reading {len(bl)} 10x barcodes')
    return(bl,bl_pre_list)
#------------------------------------------------------------------------------------
#
#  read expressed barcode summary file
#       validate each cell barcode with 10x list
#       validate each tag [expressed barcode] with whitelist
#
#  return dictionary of {cellBarcode: {tag:umiList}}
#
#
#------------------------------------------------------------------------------------
def readInputFiles(files,prefix,crBarcodeList,whitelist):
    bcDict = {}
    goodCells=0
    badCells=0
    #
    # input file format is
    #
    # cellBarcode tag umi
    #
    # read all input files and put them into a dictionary based on cell barcode
    #   { cellBarcpde : {tag : [umi]} }
    #
    # cellBarcode :
    #
    #
    for file1,px in zip(files,prefix):
        print(f' Read input file {file1} and add prefix {px}')
        df = (line.strip().split('\t') for line in open(file1,'r'))
        next(df) # skip header
        index = 0
        for nl in df:
            cellBC = nl[0]
            #
            # is cell barcode in 10x list?
            #
            if cellBC in crBarcodeList:
                goodCell = True
                goodCells += 1
                # add file prefix once library match is complete
                cellBC = px + cellBC
            else:
                badCells += 1
                goodCell = False
                #print(f'ERROR    {cellBC} not found  in 10x list')
            if goodCell:
              #
              # read expression barcode pairs of (tag,umi)
              #
              k   = nl[1]
              umi = nl[2]
              # check tag against whitelist
              if k in whitelist:
                pass
              else:
                # see if there is a matching tag with one substitution
                error,whIndex = compToWhite(k,whitelist)
                #print(error,index,k,whitelist[index])
                if error < 2:
                  #print(f'error correct {k} to {whitelist[index]}')
                  k = whitelist[whIndex]
                else:
                  print(f'ERROR: Could not correct tag {k}, min error = {error}')
                  k=''
              #
              # did tag pass whitelist check?
              #
              if k != '':
                #
                # tag to the cell
                #
                # is cell already in bcDict?
                #
                if cellBC in bcDict.keys():
                  #
                  # must combine tag : umi
                  #
                  tagDict = bcDict[cellBC]

                  try:
                    # already seen this tag, add umi
                    umiDict = tagDict[k]
                    umiDict[umi] += 1
                  except KeyError:
                    # new tag, creade umiDict and add umi
                    umiDict = defaultdict(int)
                    umiDict[umi] += 1
                    tagDict[k] = umiDict
                else:
                  # new cell BC
                  umiDict = defaultdict(int)
                  umiDict[umi] += 1
                  tagDict = {k:umiDict}
                  bcDict[cellBC] = tagDict

            index += 1
            if index % 1000 == 0:
                print(f'index = {index}, lenght of cell barcode dict = {len(bcDict)}  good cells = {goodCells} bad cells = {badCells}')
        print(f'Final index = {index}, lenght of cell barcode dict = {len(bcDict)}  good cells = {goodCells} bad cells = {badCells}')
    return(bcDict)
#-----------------------------------------------------------------------------------------
#
#
# add 10x cells not found in express barcode data
#
#
#-----------------------------------------------------------------------------------------
def add10x(cbc,bc_pre_list):
    print(f'add cells where no tags were found, initial dict size = {len(cbc)}')
    for bcl in bc_pre_list:
      print(f'barcode list size = {len(bcl)}')
      for cellBC in bcl:
        if not cellBC in cbc.keys():
            cbc[cellBC] = {}
    print(f'final cbc size = {len(cbc)}')
    return(cbc)

#-----------------------------------------------------------------------------------------
#
#
#   for each list of UMIs, compress to just the number of unique UMIS
#  
#       input is {cellBarcode : {tag:UMI}}
# 
#       output is {cellBarcode : {tag:numUniqueUMI}}
#
#-----------------------------------------------------------------------------------------
def compressUMI(d):
    ret = {}

    for cellBC,tagDict in d.items():
      td = {}
      for tag,umiDict in tagDict.items():
        td[tag] = len(umiDict)
      ret[cellBC] = td

    return(ret)
#-----------------------------------------------------------------------------------------
#
#
#   initCloneDict - make a dictionay for clones as { cloneID : [cells] }
#   starting with a unique clone for each cell
#
#
#  input is cell dictionary of {cellBC : {tag : umiCount}}
#
#-----------------------------------------------------------------------------------------
def initCloneDict(cbc):
    clones = {}
    cloneID = 0
    for cellBC, td in cbc.items():
        clones[cloneID] = [cellBC]
        cloneID += 1
    return clones
#-----------------------------------------------------------------------------------------
#
#
#   JaccardScore
#
#   calc jaccard scores baesed on umi counts of same and different tags per clone
#
#  cbc is the initail {cellBC : {tag : umiCount} }
#  clones is the initial clone dict  {cloneID : [cellList]}
#
#-----------------------------------------------------------------------------------------
def jaccardScore(cbc,cloneInput):
  with open('fred.log','w') as fh:
    for cloneID1,cellList1 in cloneInput.items():
        j_scores = []
        cell1_0 = cellList1[0]
        tags1 = list(cbc[cell1_0].keys())
        umis1 = []
        for tag in tags1:
            umis1.append(cbc[cell1_0][tag])
        for cloneID2, cellList2 in cloneInput.items():
            if cloneID1 == cloneID2:
                continue
            cell2_0 = cellList2[0]
            tags2 = list(cbc[cell2_0].keys())
            umis2 = []
            for tag in tags2:
                umis2.append(cbc[cell2_0][tag])
            # jaccard denominator
            jd = 0
            jn = 0
            for u in umis1:
                jd += u
            for u in umis2:
                jd += u
            # jaccard numerator for tags  in cell1 found in cell2
            for index,tag in enumerate(tags1):
                if tag in tags2:
                    jn += umis1[index]
            # jaccard numerator for tags in cell2 found in cell1
            for index,tag in enumerate(tags2):
                if tag in tags1:
                    jn += umis2[index]
            score = 1.0
            if jd > 0:
                score = jn / jd
            j_scores.append(score)
        for s in j_scores:
            fh.write(f'{s:2.2f}\n')
#-----------------------------------------------------------------------------------------
#
# get clone umi score
#
#
#-----------------------------------------------------------------------------------------
def getCloneUmiScore(cbc,cellList):
    umiCount = defaultdict(int)
    totalUmiCount = 0
    for cell in cellList:
        tagUMI = cbc[cell]
        for tag,umi in tagUMI.items():
            umiCount[tag] += umi
            totalUmiCount += umi
    tagList = {}
    for tag,umi in umiCount.items():
        umi = umi/totalUmiCount
        tagList[tag] = umi
    return(tagList)



#-----------------------------------------------------------------------------------------
#
#
#   JaccardScore
#
#   calc jaccard scores baesed on umi counts of same and different tags per clone
#
#  cbc is the initail {cellBC : {tag : umiCount} }
#  clones is the initial clone dict  {cloneID : [cellList]}
#
#-----------------------------------------------------------------------------------------
def jaccardMerge(cbc,clones,minScore,oPrefix):
  oFile = 'output/' + oPrefix + '_jaccard_merge.tsv'
  with open(oFile,'w') as fh:
    cloneInput = clones
    newClones = {}
    merges = 0
    index = 0
    cloneIDs = sorted(cloneInput.keys())
    for cloneID1 in cloneIDs:
        cellList1 = cloneInput[cloneID1]
        if len(cellList1) == 0:
            continue
        newClones[cloneID1] = cellList1
        tagD1 = getCloneUmiScore(cbc,cellList1)
        if len(tagD1.values()) == 0:
                continue
        #print(f'----{cloneID1}------------------')
        #print(tagD1)
        for cloneID2 in cloneIDs:
            #print(cloneID1,cloneID2)
            cellList2 = cloneInput[cloneID2]
            # only merge smaller clones into larger
            if cloneID2 <= cloneID1:
                continue
            if len(cellList2) == 0:
                continue
            tagD2 = getCloneUmiScore(cbc,cellList2)
            #print(f'+++++++++++++++{cloneID2}++++++++')
            #print(f'    {tagD2}')
            # jaccard denominator
            jd = 0
            jn = 0
            tags1 = tagD1.keys()
            tags2 = tagD2.keys()
            for t in tags1:
                if t in tags2:
                    jn += tagD2[t]
                    jn += tagD1[t]
                    #print(f'tag match {t}')
            # jaccard numerator
            for v in tagD1.values():
                jd += v
            for v in tagD1.values():
                jd += v
            score = 1.0
            if jd > 0:
                score = jn / jd
            #print(jn,jd,score)
            # merge
            if score >= minScore:
                # combine clones
                newClones[cloneID1].extend(cellList2)
                cloneInput[cloneID2] = []
                merges += 1
                fh.write(f'MERGE   {score}  {cloneID1}   {cloneID2}\n')
                for cell in cellList1:
                    fh.write(f'  clone1  {cbc[cell]}\n')
                fh.write('------\n')
                for cell in cellList2:
                    fh.write(f'  clone2  {cbc[cell]}\n')
                fh.write('\n')
        #index += 1
        #if index > 2: break
    print(f'JaccardMerge: merged {merges} clones : from {len(clones)} to {len(newClones)}')
    tc = 0
    for k,v in newClones.items():
        tc = tc + len(v)
    print(f'JaccardMerge: final cells = {tc}')
    ## finally sort by clone size
    sd = sorted(newClones.items(),key=lambda x:len(x[1]),reverse=True)
    rd = {}
    newCloneIndex = 0
    for k,v in sd:
        rd[newCloneIndex] = v
        newCloneIndex += 1
    #return(rd)
    return(newClones)
#-----------------------------------------------------------------------------------------
#
#
#   combineExactTag
#
#   combine clones where all cells have same tag(s)
#
#  cbc is the initail {cellBC : {tag : umiCount} }
#  clones is the initial clone dict  {cloneID : [cellList]}
#
#  returns {cloneID:[cells]}
#
#-----------------------------------------------------------------------------------------
def combineExactTag(cbc,cloneInput):
    merges = 0
    clones = cloneInput  # do not destroy
    newClones = {}
    for cloneID1,cellList1 in clones.items():
        # clone not already merged
        if len(cellList1) == 0:
            continue
        newClones[cloneID1] = cellList1
        # only identical tags are looked at now so first cell can provice
        tags1 = cbc[cellList1[0]].keys()

        for cloneID2, cellList2 in clones.items():
            if cloneID1 != cloneID2 and len(cellList2) > 0:
                tags2 = cbc[cellList2[0]].keys()
                if tags1 == tags2:
                    # combine clones
                    newClones[cloneID1].extend(cellList2)
                    clones[cloneID2] = []
                    merges += 1
    print(f'combineExactTag: merged {merges} clones : from {len(cloneInput)} to {len(newClones)}')
    print(f'combineExactTag: total cells at start = {len(cloneInput)}')
    tc = 0
    for k,v in newClones.items():
        tc = tc + len(v)
    print(f'combineExactTag: final cells = {tc}')
    # finally sort by clone size
    sd = sorted(newClones.items(),key=lambda x:len(x[1]),reverse=True)
    rd = {}
    newCloneIndex = 0
    for k,v in sd:
        rd[newCloneIndex] = v
        newCloneIndex += 1
    return(rd)

#-----------------------------------------------------------------------------------------
#
#
# helper function...how many strings to l1 and l2 have in common
#
#
#-----------------------------------------------------------------------------------------
def tagsInCommon(l1,l2):
  count = 0
  for t in l1:
    if t in l2:
      count += 1
  return(count)
#-----------------------------------------------------------------------------------------
#
#
# helper function...get list of tags with umi > 1
#
# tagD {tag:umi}
#-----------------------------------------------------------------------------------------
def getLargeTags(tagD):
    ret = []
    for tag,umi in tagD.items():
        if umi > 1:
            ret.append(tag)
    return(ret)



def getLargeAveTags(cbc,cellList):
  tagSum = defaultdict(int)
  for cell in cellList:
    #print(f'  {cell}')
    tagD = cbc[cell]
    for tag,umi in tagD.items():
      #print(f'    tag = {tag}, umi = {umi}')
      tagSum[tag] += umi
  tagAve = defaultdict(int)
  for tag,umiCount in tagSum.items():
    tagAve[tag] = umiCount / len(cellList)

  ret = []
  rumi = []
  for tag,umiAve in tagAve.items():
    if umiAve > 1.0:
      ret.append(tag)
      rumi.append(umiAve)

  #print(f'tagSum: {tagSum}')
  #print(f'tagAve: {tagAve}')
  #print(f'ret:    {ret}')
  #value = input("Please enter a string:\n")
  return((ret,rumi))

#-----------------------------------------------------------------------------------------
#
#
#   merge clones that all have three matching tags
#
#
#-----------------------------------------------------------------------------------------
def combineMatching(cbc,cloneInput,matchCount):
    print(f'combineMatching {matchCount}: ')
    merges = 0
    clones = cloneInput  # do not destroy
    newClones = {}
    for cloneID1,cellList1 in clones.items():
        # clone not already merged
        if len(cellList1) == 0:
            continue
        newClones[cloneID1] = cellList1
        # get tags > 1 umi from first cell
        #tags1 = getLargeTags(cbc[cellList1[0]])
        tags1,umi1 = getLargeAveTags(cbc,cellList1)

        for cloneID2, cellList2 in clones.items():
            # clone not already merged
            if cloneID1 != cloneID2 and len(cellList2) > 0:
                # get tag list from first cell
                #tags2 = getLargeTags(cbc[cellList2[0]])
                tags2,umi2 = getLargeAveTags(cbc,cellList2)
                if tagsInCommon(tags1,tags2) >= matchCount:
                    # combine clones
                    if False:
                      print(f'    Merge')
                      print(f'        {tags1}   {len(cellList1)}')
                      print(f'        {umi1}    {len(cellList1)}')
                      print(f'        {tags2}   {len(cellList2)}')
                      print(f'        {umi2}    {len(cellList2)}')
                      for cell in cellList1:
                        print(f'        cell1 {cell}')
                        for tag,umi in cbc[cell].items():
                          print(f'            {tag}:{umi}')


                      for cell in cellList2:
                        print(f'        cell2 {cell}')
                        for tag,umi in cbc[cell].items():
                          print(f'            {tag}:{umi}')

                      value = input("Please enter a string:\n")
                    # end of debug
                    newClones[cloneID1].extend(cellList2)
                    clones[cloneID2] = []
                    merges += 1
    print(f'combineMatching {matchCount}: merged {merges} clones : from {len(cloneInput)} to {len(newClones)}')
    tc = 0
    for k,v in newClones.items():
        tc = tc + len(v)
    print(f'combineMatching {matchCount}: final cells = {tc}')
    # finally sort by clone size
    sd = sorted(newClones.items(),key=lambda x:len(x[1]),reverse=True)
    rd = {}
    newCloneIndex = 0
    for k,v in sd:
        rd[newCloneIndex] = v
        newCloneIndex += 1
    return(rd)
#-----------------------------------------------------------------------------------------
#
#
#   generate clone-cell-tag report
#
#     cbc is the initial data structure will cell - tag info
#     cloneDict is the current merged clone - cell list
#
#
#-----------------------------------------------------------------------------------------
def generateCellTagReport(cbc,cloneDict,oPrefix):
    oFile = 'output/' + oPrefix + '_cell_tag_report.tsv'
    with open(oFile,'w') as fh:
        fh.write(f'CloneID\tCellBarcode\tExpBarcode\tUMICount\n')
        for cloneID,cellList in cloneDict.items():
            fh.write(f'# Clone ID {cloneID} number of cells = {len(cellList)}\n')
            for cell in cellList:
              #fh.write(f'{cell}: ')
                tags = sorted(list(cbc[cell].keys()))
                if len(tags) == 0:
                  fh.write(f'{cloneID}\t{cell}\tNone\t1\n')
                else:
                  for tag in tags:
                    umi = cbc[cell][tag]
                    fh.write(f'{cloneID}\t{cell}\t{tag}\t{umi:5d}\n')



#-----------------------------------------------------------------------------------------
#
#
#   generate human-readable clone-cell-tag report
#
#     cbc is the initial data structure will cell - tag info
#     cloneDict is the current merged clone - cell list
#
#
#  CLONE:
#    CELL1     tagA:count     tagB:count  
#    CELL2                    tagB:count    tagC:count
#
#-----------------------------------------------------------------------------------------
#
# helper----get sorted list of all tags clone
#
def getTagList(cbc,cellList):
  allTags = defaultdict(int)
  for cell in cellList:
    tags = sorted(list(cbc[cell].keys()))
    for tag in tags:
      allTags[tag] += 1
  #
  # sort to find tags with most counts
  #   return tags in order of most common
  #
  sd = sorted(allTags.items(),key=lambda x:x[1],reverse=True)
  rt = []
  for t,c in sd:
    rt.append(t)

  return(rt)


def generateHumanCellTagReport(cbc,cloneDict,oPrefix):
    oFile = 'output/' + oPrefix + '_human_cell_tag_report.tsv'
    with open(oFile,'w') as fh:
        fh.write(f'CloneID\tCellBarcode\tExpBarcode\tUMICount\n')
        for cloneID,cellList in cloneDict.items():
            fh.write(f'# Clone ID {cloneID} number of cells = ,{len(cellList)}\n')
            tagList = getTagList(cbc,cellList)
            for cell in cellList:
              fh.write(f'{cell}   ')
              tags = cbc[cell].keys()
              for t in tagList:
                if t in tags:
                  umi = cbc[cell][t]
                  fh.write(f'{t},{umi:5d}, ')
                else:
                  fh.write(f'________,..... ,')
              fh.write('\n')


#-----------------------------------------------------------------------------------------
#
#
# generate tag report
#    for each tag find all cell-clone it is used in
#   cbc is:
#     { cellBarcode : {tag : [umi]} }
#
#   cloneDict is:
#     { cloneID : [cells] }
#-----------------------------------------------------------------------------------------

def generateTagReport(cbc,cloneDict,oPrefix):
  #
  # make tag dict as:
  #
  #   {tag : [cellBarcodes]
  #
  tg = {}
  for cellBarcode,tagDict in cbc.items():
    for tag in tagDict.keys():
      try:
        tg[tag].append(cellBarcode)
      except KeyError:
        tg[tag] = [cellBarcode]


  print(f'tag dir has {len(tg)} tags')
  #
  # sort for largest tags fist
  #
  sd = sorted(tg.items(),key = lambda x:len(x[1]),reverse=True)
  #with open('fred.txt','w') as fh:
  #  for k,v in sd:
  #    fh.write(f'{k},{v}\n')
  #
  # from cloneDict make dict of cellBarcode:cloneID
  #
  cd = defaultdict(int)
  for cloneID,cellList in cloneDict.items():
    for cell in cellList:
      cd[cell] = cloneID
  print(f'cell list lenght = {len(cd)}')
  #
  #scd = sorted(cd.items(),key=lambda x:x[1])
  #with open('barney.txt','w') as fh:
  #  for k,v in scd:
  #    fh.write(f'{k},{v}\n')
  #
  # add clone ID to cell barcodes and write
  #
  #
  oFile = 'output/' + oPrefix + '_tag_cell_clone_report2.tsv'
  with open(oFile,'w') as fh:
    fh.write('Tag\tnClones\tnCells\n')
    for tag,cellList in sd:
      # find out how many cells and how many unique clones
      uClone = defaultdict(int)
      for cell in cellList:
        cell_clone_id = cd[cell]
        uClone[cell_clone_id] += 1
      lenCell = len(cellList)
      lenClone = len(uClone)
      fh.write(f'{tag}\t{lenClone}\t{lenCell}\n')

#-----------------------------------------------------------------------------------------
#
#
#   generate clone output
#
#   cell barcodes are column labels
#   row 1 is cloneID for each barcode
#
#
#-----------------------------------------------------------------------------------------
def generateOutputFile(cbc,cloneDict,oPrefix):
    oFile = 'output/' + oPrefix + '_10x_clone_barcodes.csv'
    with open(oFile,'w') as fh:
        bca = []
        ida = []
        for cloneID,cellList in cloneDict.items():
          for cell in cellList:
            bca.append(cell)
            ida.append(cloneID)

        fh.write(f'{bca[0]}')
        for bc in bca[1:]:
          fh.write(f',{bc}')
        fh.write('\n')

        fh.write('CloneID')
        for cid in ida:
          fh.write(f',{cid}')
        fh.write('\n')

#def generateOutputFile(cloneDict,barcodeList,oPrefix):
#    oFile = 'output/' + oPrefix + '_10x_clone_barcodes.csv'
#    with open(oFile,'w') as fh:
#        first = True
#        for bc in barcodeList:
#            if first:
#                fh.write(f'"{bc}"')
#                first = False
#            else:
#                fh.write(f',"{bc}"')
#        fh.write('\n')
#
#        fh.write('CloneID')
#        for bc in barcodeList:
#            # find id for cell
#            found = False
#            for cloneID,cellList in cloneDict.items():
#                if bc in cellList:
#                    fh.write(f',{cloneID}')
#                    found = True
#                    break
#            if found == False:
#                print(f"ERROR, could not find {bc} in clone dict")
#        fh.write('\n')
#-----------------------------------------------------------------------------------------
#
#
#   ErorrCorrectTags
#
#
#-----------------------------------------------------------------------------------------
#
# error-correct tags   off-by-one assign to other tag?
#
def tagDiff(seq1,seq2):
  diff = 0
  for i in range(0,8):
    if seq1[i] != seq2[i]:
      diff += 1
  return diff

def errorCorrectTags(cbc):
    errorCount = 0
    for cell,tagD in cbc.items():
        tp = []
        for tag,umi in tagD.items():
            tp.append((tag,umi))
        for tag,umi in tp:
            if umi <= 1:
                # try to correct
                for tag2,umi2 in tp:
                    if  tag2 != tag and umi2 > 1:
                        diff = tagDiff(tag,tag2)
                        if diff <= 1:
                            # correct tag
                            tagD[tag2] += 1
                            tagD[tag] = 0
                            errorCount += 1
                            #print(f' {cell}  repair {tag} {umi} to {tag2} {umi2}   {tagD[tag]} {tagD[tag2]}')
                            break
    #
    # remove empty tags
    #
    for cell,tagD in cbc.items():
        newTagD = {}
        for tag,umi in tagD.items():
            if umi > 0:
                newTagD[tag] = umi
        cbc[cell] = newTagD
    print(f'Error Correction fixed {errorCount} 1 base pair errors')
#-----------------------------------------------------------------------------------------
#
#
#   Main Processing
#
#
#-----------------------------------------------------------------------------------------
#
# read expressed barcode whitelist and
# read 10x valid barcode lists (from 10x outs directory)
#

whitelist = readWhiteList(whiteFile)
crBarcodeList,bc_pre_list = readBarcodeList(barfiles,prefix)
cellBcDict = readInputFiles(files,prefix,crBarcodeList,whitelist)
#
# compress unique UMI in cellBcDict [one tag for all umi]
#
cellBcDict = compressUMI(cellBcDict)
#
# add cell barcodes with no tags detected
#
cellBcDict = add10x(cellBcDict,bc_pre_list)
#
# error-correct tags
#
errorCorrectTags(cellBcDict)
#
# initialize a clone dictionary
#
# the format for the clone dictionary is { cloneID : [cells] }
#
# the starting format of the clone dictionary is each clone is one cell
#
cloneDict1 = initCloneDict(cellBcDict)
print(f'Initial Clone dictionary length = {len(cloneDict1)}')
#
#
# now combine clones with cells with exact tag matches
#
#
cloneDict2 = combineExactTag(cellBcDict,cloneDict1)
print(f'Exact match Clone dictionary length = {len(cloneDict2)}')
cloneDict3 = combineMatching(cellBcDict,cloneDict2,3)
print(f'Combine matching tags: Clone dictionary length = {len(cloneDict3)}')
cloneDict4 = combineMatching(cellBcDict,cloneDict3,5)
print(f'Combine matching tags: Clone dictionary length = {len(cloneDict4)}')
cloneDict5 = combineMatching(cellBcDict,cloneDict4,4)
print(f'Combine matching tags: Clone dictionary length = {len(cloneDict5)}')

cloneDict6 = jaccardMerge(cellBcDict,cloneDict5,0.95,out)
print(f'Jaccard 0.95 marge: Clone dictionary length = {len(cloneDict6)}')
#
#
# generate cloneID output
#
#
generateOutputFile(cellBcDict,cloneDict6,out)
generateCellTagReport(cellBcDict,cloneDict6,out)
generateHumanCellTagReport(cellBcDict,cloneDict6,out)
#
# for each tag, find all cell-cloneID it is in
#
generateTagReport(cellBcDict,cloneDict6,out)







