import glob,sys
'''
takes the filtered bicseq file and translates it to a format that is readable by annovar
'''
infile = sys.argv[1]
outfile = sys.argv[2]
f = open(infile, 'rU')
o = open(outfile, 'w')
for line in f:
    if line.startswith('#'):
        continue
        o.write(line.strip() + '\n')
    else:
        splits = line.split()
        o.write(splits[0] + '\t' + splits[1] + '\t' + splits[2] + '\t0\t0\tCOMMENTS:[' + line.strip() + ']\n')
f.close()
o.close()
