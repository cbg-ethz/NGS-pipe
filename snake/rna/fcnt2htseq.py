import sys

fcfile = open(sys.argv[1],'rU')
htsfile = open(sys.argv[2],'w')


for line in fcfile:
    if line.startswith('#') or line.startswith('Gene'):
        continue
    split = line.split()
    ens = split[0].strip()
    cnt = split[6]
    htsfile.write(ens + ' ' + str(int(round(float(cnt), 0))) + '\n')
fcfile.close()
htsfile.close()

