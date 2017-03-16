#Take mapping file with filename from htseq files -> condition as input
#Take output filename as input
#Create outfile matrix that can be used for DeSeq2
import sys

if len(sys.argv) != 3:
    print('usage: python DeSeq2Preparer mappingfile outfileDir')
    print('Where mappingfile is a file that maps the htseq input files to a condition in a tab separated file')
    exit()



mappingfile = sys.argv[1]
outfileDir = sys.argv[2]


#read mapping file
mfile = open(mappingfile, 'rU')
folder = ''
files = []
compares = []
while True:
    line = mfile.readline()
    if not line:
        break
    if line.startswith('Folder:'):
        folder = mfile.readline().strip()
        line = folder
    if line.startswith('Files:'):
        line = mfile.readline().strip()
        while not line.startswith('DeSeq2Calls:'):
            files.append(line.strip().split()[0])
            line = mfile.readline().strip()
    if line.startswith('DeSeq2Calls:'):
        line = mfile.readline().strip()
        while True:
            if not line:
                break
            compares.append((line.split()[0].strip().split(','),line.split()[1].strip().split(','), line.split()[2].strip()))
            line = mfile.readline().strip()

#done reading the mapping file


all = {}
pos = 0
for file in files:
    f = open(folder + '/' + file, 'rU')
    for line in f:
        name = line.split('\t')[0].strip()
        count = int(line.split('\t')[1].strip())
        if count == 0:
            continue
        if name not in all:
            all[name] = [0] * len(files)
        all[name][pos] = int(count)
    pos = pos + 1
    f.close()



for compare in compares:
    a = compare[0]
    b = compare[1]
    name = compare[2]
    outfile = open(outfileDir + '/' + name + '.txt', 'w')
    outfile.write('Gene\t')
    for x in a:
        outfile.write('A\t')
    for x in b:
        outfile.write('B\t')
    outfile.write('\n')

    for gene,counts in all.iteritems():
        newcounts = [0]*(len(a) + len(b))
        pos = 0
        for x in a:
            newcounts[pos] = counts[files.index(x)]
            pos = pos + 1
        for x in b:
            newcounts[pos] = counts[files.index(x)]
            pos = pos + 1
        if sum(newcounts) == 0:
            continue
        outfile.write(gene + '\t')
        for value in newcounts:
            outfile.write(str(value) + '\t')
        outfile.write('\n')
    outfile.close()
