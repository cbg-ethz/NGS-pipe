import gzip,sys


deseqfile = open(sys.argv[1],'rU')#open('/Users/hansr/Downloads/xeno/new/xeno_liver.results.tsv', 'rU')
ensemblidfile = gzip.open(sys.argv[2], 'rU')#gzip.open('/Users/hansr/Desktop/EnsembleId2GeneName.csv.gz', 'rU')
outfile = open(sys.argv[3],'w')#open('/Users/hansr/Downloads/xeno/new/xeno_liver.results_geneNames.tsv', 'w')

id2name = {}
for line in ensemblidfile:
    if line.startswith('#'):
        continue
    id2name[line.split('\t')[0].strip()] = line.split('\t')[0].strip() + '__' + line.split('\t')[1].strip() + '__' + line.split('\t')[2].strip()
ensemblidfile.close()

for line in deseqfile:
    if line.startswith('base'):
        outfile.write(line)
        continue
    ensembleid = line.split('\t')[6].strip()
    gene_name = ensembleid
    if ensembleid in id2name:
        gene_name = id2name[ensembleid]
    line = line.replace(ensembleid, gene_name)
    outfile.write(line)
deseqfile.close()
outfile.close()
