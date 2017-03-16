#!/bin/bash
echo -e '##fileformat=VCFv4.1'
echo -e '##INFO=<ID=NR,Number=1,Type=Integer,Description="Number of reads supporting reference in normal">'
echo -e '##INFO=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting variant in normal">'
echo -e '##INFO=<ID=TR,Number=1,Type=Integer,Description="Number of reads supporting reference in tumor">'
echo -e '##INFO=<ID=TV,Number=1,Type=Integer,Description="Number of reads supporting variant in tumor">'
echo -e '##INFO=<ID=RPS,Number=.,Type=Float,Description="Raw probability of somatic mutation P=p_AA_AB+p_AA_BB">'
echo -e '##INFO=<ID=PS,Number=.,Type=Float,Description="Post processed probability of somatic mutation">'
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

fprintf(fd_out, ".\t%s\t%s\t%f\t.\tNR=%s;NV=%s;TR=%s;TV=%s;RPS=%f;PS=%f\t.\n",   ref, var, p_somatic,  normal_a,normal_b, tumor_a, tumor_b,PS,p_somatic);

cat ${jsm} | awk -F "\t" 'NR!=1' | \
    awk -F "\t" '{print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\t.\tAAAB=" $10 ";AABB=" $11 "\tRD:AD\t" $5 ":" $6 "\t" $7 ":" $8}' | \
    sort -k1,1V -k2,2n
