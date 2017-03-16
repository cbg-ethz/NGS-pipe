import io, sys, os

sam = sys.stdin
countAll = 0
countPaired = 0
countMapped = 0
countUnmapped = 0
countReadMappedInProperPair = 0
countMateUnmapped = 0
countReadReverse = 0
countMateReverse = 0
countFirstInPair = 0
countSecondInPair = 0
countSecondaryAlignment = 0
countFailedQualityCheck = 0
countPCROrOpticalDuplicated = 0
countSupplementaryAlignment = 0


for line in sam:
    split = line.split('\t')
    flag = int(split[1].strip())
    countAll += 1
    if flag & 1:
        #print 'read paired'
        countPaired += 1
    if flag & 4:
        #print 'read unmapped'
        countUnmapped += 1
    else:
        countMapped += 1
    if flag & 2:
        #print 'read mapped in proper pair'
        countReadMappedInProperPair += 1
    if flag & 8:
        #print 'mate unmapped'
        countMateUnmapped += 1
    if flag & 16:
        #print 'read reverse'
        countReadReverse += 1
    if flag & 32:
        #print 'mate reverse'
        countMateReverse += 1
    if flag & 64:
        #print 'first in pair'
        countFirstInPair += 1
    if flag & 128:
        #print 'second in pair'
        countSecondInPair += 1
    if flag & 256:
        #print 'secondary alignment'
        countSecondaryAlignment += 1
    if flag & 512:
        #print 'failed quality check'
        countFailedQualityCheck += 1
    if flag & 1024:
        #print 'pcr or optical duplicate'
        countPCROrOpticalDuplicated += 1
    if flag & 2048:
        #print 'supplementary alignment'
        countSupplementaryAlignment += 1
print 'Alignments\t' + str(countAll)
print 'Paired reads\t' + str(countPaired)
print 'Aligned reads\t' + str(countMapped)
print 'Unaligned reads\t' + str(countUnmapped)
print 'Read mapped in proper pair\t' + str(countReadMappedInProperPair)
print 'Unaligned Mate\t' + str(countMateUnmapped)
print 'Read reverse strand\t' + str(countReadReverse)
print 'Mate reverse strand\t' + str(countMateReverse)
print 'First in Pair\t' + str(countFirstInPair)
print 'Second in Pair\t' + str(countSecondInPair)
print 'Primary Alignments\t' + str(countAll - countUnmapped - countSecondaryAlignment - countFailedQualityCheck)
print 'Secondary Alignments\t' + str(countSecondaryAlignment)
print 'Failed Quality Check\t' + str(countFailedQualityCheck)
print 'PCR or Optical Duplicate\t' + str(countPCROrOpticalDuplicated)
print 'Supplementary Alignments\t' + str(countSupplementaryAlignment)
