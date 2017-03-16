
file = open('/Users/hansr/Desktop/Homo_sapiens.GRCh37.75.gtf', 'rU')
hgncfile = open('/Users/hansr/Desktop/hgnc_complete_set.txt', 'rU')
outfile = open('/Users/hansr/Desktop/EnsembleId2GeneName.csv', 'w')

id2symbol = {}

id2name = {}
for line in hgncfile:
    if line.startswith('hgnc_id'):
        continue
    splits = line.split('\t')
    symbol = splits[1].strip()
    name = splits[2].strip()
    id2name[symbol] = name





for line in file:
    if line.startswith('#'):
        continue
    splits = line.split()
    ensemblid = splits[splits.index('gene_id')+1].strip()
    gene_name = splits[splits.index('gene_name')+1].strip()
    ensemblid = ensemblid[1:len(ensemblid)-2]
    gene_name = gene_name[1:len(gene_name)-2]
    if ensemblid not in id2symbol:
        id2symbol[ensemblid] = set()
    id2symbol[ensemblid].add(gene_name)
outfile.write('#generated from Homo_sapiens.GRCh37.75.gtf\n#Ensembleid\tgene_name\n')




for id, gene_names in id2symbol.iteritems():
    symbol = gene_names.pop()
    name = 'No HGNC Name'
    if symbol in id2name:
        name = id2name[symbol]
    outfile.write(id + '\t' + symbol + '\t' + name + '\n')
outfile.close()
