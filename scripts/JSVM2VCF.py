import glob,sys
from math import log
'''

'''
infile = sys.argv[1]
outfile = sys.argv[2]
f = open(infile, 'rU')
o = open(outfile, 'w')
o.write("##fileformat=VCFv4.0\n")
o.write("##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads supporting reference in normal\">\n")
o.write("##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in normal\">\n")
o.write("##INFO=<ID=TR,Number=1,Type=Integer,Description=\"Number of reads supporting reference in tumor\">\n")
o.write("##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in tumor\">\n")
o.write("##INFO=<ID=FREQ,Number=.,Type=Float,Description=\"Variant frequency\">\n")
o.write("##INFO=<ID=RPS,Number=.,Type=Float,Description=\"Raw probability of somatic mutation P=p_AA_AB+p_AA_BB\">\n")
o.write("##INFO=<ID=PS,Number=.,Type=Float,Description=\"Post processed probability of somatic mutation\">\n")
o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")

f.readline() # skip the first line

for line in f:
    splits = line.split()
    o.write(splits[0] + '\t')    # chr
    o.write(splits[1] + '\t')    # pos
    o.write('.\t')               # identifier
    o.write(splits[2] + '\t')    # ref
    o.write(splits[3] + '\t')    # alt

    pSomatic = float(splits[17])
    if (pSomatic != 0):
        if (pSomatic<1.0):
            pSomatic = -10.0 * log(1.0 - pSomatic)
        else:
            pSomatic = 255.0
        if (pSomatic > 255.0):
            pSomatic = 255.0

    o.write(str(pSomatic) + '\t')       # qual
    o.write('.\t')                      # filter
    o.write('NR=' + splits[4])          # NR
    o.write(';NV=' + splits[5])         # NV
    o.write(';TR=' + splits[6])         # AN
    o.write(';TV=' + splits[7])         # AV
    
    freq = float(splits[7]) / (float(splits[6]) + float(splits[7]))

    o.write(';FREQ=' + str(freq))       # FREQ

    rsp = float(splits[9]) + float(splits[10])

    o.write(';RSP=' + str(rsp))   # RSP
    o.write(';SP=' + splits[17])  # SP
    o.write('\t.\t.\n')              # FORMAT field

f.close()
o.close()
