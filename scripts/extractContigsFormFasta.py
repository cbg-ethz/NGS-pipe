#!/cluster/work/bewi/modules/biopython/1.66/bin/python3
import sys, os

refFasta = sys.argv[1]

from Bio import SeqIO
handle = open(refFasta, "r")
if not os.path.exists(refFasta + '_contigs/'):
    os.makedirs(refFasta + '_contigs/')
for record in SeqIO.parse(handle, "fasta") :
    out = open(refFasta + "_contigs/" + record.id + ".fasta", "w")
    SeqIO.write(record, out, "fasta")
handle.close()
