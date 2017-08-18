# NGS-Pipe

## Description
NGS-Pipe provides analyses for large scale DNA and RNA sequencing experiments.
The scope of pre-implemented functions spans the detection of germline variants,
somatic single nucleotide variants (SNV) and insertion and deletion (InDel)
identification, copy number event detection, and differential expression
analyses. Further, it provides pre-configured workflows, such that the final
mutational information as well as quality reports and all intermediate results
can be generated quickly, also by inexperienced users. In addition, the pipeline
can be used on a single computer or in a cluster environment where independent
steps are executed in parallel. If one of the steps of the pipeline fails and
produces incomplete or no results, the computation of all depending steps is
halted and an error message is shown. However, after the issue is resolved the
pipeline independently resumes the analyses at the appropriate point,
eliminating the need to rerun the complete analysis or manual deletion of
erroneous files.

See also the wiki pages of this repository for more information about NGS-pipe.

## Workflows for WES, WGS, and RNA-seq data
We have implemented and tested predefined workflows for the automated analysis
of WES, WGS, and RNA-seq data (Fig. 1).

<h1 align="center">
<img src="https://github.com/cbg-ethz/NGS-pipe/blob/master/docs/pipe_step_by_step_compact.png?raw=true" alt="Workflows"/></h1>

The primary data analysis steps include Trimmomatic (Bolger, 2014) to process
raw files, BWA (Li, 2009) or STAR (Dobin, 2013) to align reads, and Picard
tools (http://broadinstitute.github.io/picard), SAMtools (Li, 2009Samtools) and
GATK (McKenna, 2010) to process the aligned reads. 

Detecting genomic variants is highly dependent on properties of the input data,
such as variant frequency, coverage, or contamination (Cai, 2016; Hofmann,
2017). For this reason, we included several variant callers in NGS-pipe, viz.
Mutect (Cibulskis, 2013), JointSNVMix2 (Roth, 2012), VarScan2 (Koboldt, 2012),
VarDict (Lai, 2016), SomaticSniper (Larson, 2011), Strelka (Saunder, 2012), and
deepSNV (Gerstung, 2012). Further, we included SomaticSeq (Fang, 2015), which
combines the results of multiple variant callers and ranked high in the
ICGC-TCGA DREAM Somatic Mutation Calling Challenge (Ewing, 2015), and the rank
aggregation scheme introduced in (Hofmann, 2017).

Copy number events are detected by FACETS (Shen, 2016), or BIC-seq2 (Xi, 2016),
which has been designed specifically for whole genome data.

The results of the experiments can be annotated and manipulated using SnpEff
(Cingolani 2012), SnpSift (Cingolani, 2012) and ANNOVAR (Wang, 2010).

RNA-seq data is analyzed to quantify gene expression levels. We include quality
control, alignment, and gene counting using the SubRead (Liao, 2014) package.
Output files are reformatted to serve as direct input to tools that perform
differential gene expression analysis.

## Example
The directory *examples/wes/* contains a ready to go example for the analysis
of three leukemia patients (Cifola, 2015). This example downloads 
tumor-control matched exome data sets from the Sequence Read Archive, installs
the required programs, downloads the necessary reference files and builds the
essentials indices. Afterwards, an analysis starting with the mapping of the
reads via BWA (Li 2009) all the way to the somatic variant calling with
VarScan2 (Koboldt 2012). After the [installation of all tools via conda](https://github.com/cbg-ethz/NGS-pipe/wiki/Installation-&-Testing)
you can proceed like:
```sh
#1. Go to examples folder:
cd examples/dna
#2. Download test data: We provide an additional snakemake pipeline to 
#   download test sequences, databases and adapter files:
./run_prepare_data_locally.sh
# This will download 6 test data sets, the adapters, regions file,
# the human reference and build the BWA database index
#3. Execute the DNA Pipeline:
./run_analysis_locally.sh
# This will execute: RAW --> QC(Trimmomatic) --> Mapping(BWA) --> Sort(Picard)
# --> Merge(Picard) --> Remove Secondary Alignments(Samtools) --> MarkDuplicates(Picard)
# --> RemoveDuplicates(Samtools) --> SNV Calling (VarScan2)
```

An example for RNA-seq data analysis can be found in *examples/rna/* and [here](https://github.com/cbg-ethz/NGS-pipe/wiki/Installation-&-Testing).

## References
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer
for Illumina sequence data. Bioinformatics, 30(4), 2114-2120.

Cai, L., Yuan, W., Zhou Zhang, L. H., & Chou, K. C. (2016). In-depth comparison
of somatic point mutation callers based on different tumor next-generation
sequencing depth data. Scientific reports, 6.

Cibulskis, K., Lawrence, M. S., Carter, S. L., Sivachenko, A., Jaffe, D.,
Sougnez, C., ... & Getz, G. (2013). Sensitive detection of somatic point
mutations in impure and heterogeneous cancer samples. Nature biotechnology,
31(3), 213-219.  ISO 690	

Cifola, I., Lionetti, M., Pinatel, E., Todoerti, K., Mangano, E., Pietrelli. A.
et al. (2015). Whole-exome sequencing of primary plasma cell leukemia discloses
heterogeneous mutational patterns. Oncotarget, 6(19), 17543-17558.

Cingolani, P., Patel, V. M., Coon, M., Nguyen, T., Land, S. J., Ruden, D. M., &
Lu, X. (2012). Using Drosophila melanogaster as a model for genotoxic chemical
mutational studies with a new program, SnpSift. Toxicogenomics in non-mammalian
species, 3, 35.

Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., ... &
Ruden, D. M. (2012). A program for annotating and predicting the effects of
single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila
melanogaster strain w1118; iso-2; iso-3. Fly, 6(2), 80-92.

Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S.,
... & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner.
Bioinformatics, 29(1), 15-21.

Ewing, A. D., Houlahan, K. E., Hu, Y., Ellrott, K., Caloian, C., Yamaguchi, T.
N., ... & Calling, I. T. D. S. M. (2015). Combining tumor genome simulation
with crowdsourcing to benchmark somatic single-nucleotide-variant detection.
Nature methods, 12(7), 623-630.

Fang, L. T., Afshar, P. T., Chhibber, A., Mohiyuddin, M., Fan, Y., Mu, J. C.,
... & Koboldt, D. C. (2015). An ensemble approach to accurately detect somatic
mutations using SomaticSeq. Genome biology, 16(1), 197.

Gerstung, M., Beisel, C., Rechsteiner, M., Wild, P., Schraml, P., Moch, H., &
Beerenwinkel, N. (2012). Reliable detection of subclonal single-nucleotide
variants in tumour cell populations. Nature communications, 3, 811.

Hofmann, A. L., Behr, J., Singer, J., Kuipers, J., Beisel, C., Schraml, P., ...
& Beerenwinkel, N. (2017). Detailed simulation of cancer exome sequencing data
reveals differences and common limitations of variant callers. BMC
Bioinformatics, 18(1), 8.

Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller,
C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and
copy number alteration discovery in cancer by exome sequencing. Genome
Research, 22(3), 568--576.

Koster, J. and Rahmann, S. (2012). Snakemake–a scalable bioinformatics workflow
engine. Bioinformatics, 28(19), 2520–2522.

Lai, Z., Markovets, A., Ahdesmaki, M., Chapman, B., Hofmann, O., McEwen, R.,
... & Dry, J. R. (2016). VarDict: a novel and versatile variant caller for
next-generation sequencing in cancer research. Nucleic acids research, 44(11),
e108-e108.

Larson, D. E., Harris, C. C., Chen, K., Koboldt, D. C., Abbott, T. E., Dooling,
D. J., ... & Ding, L. (2011). SomaticSniper: identification of somatic point
mutations in whole genome sequencing data. Bioinformatics, 28(3), 311-317.

Li H. and Durbin R. (2009). Fast and accurate short read alignment with
Burrows-Wheeler Transform. Bioinformatics, 25(14), 1754-1760.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... &
Durbin, R. (2009). The sequence alignment/map format and SAMtools.
Bioinformatics, 25(16), 2078-2079.

Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general
purpose program for assigning sequence reads to genomic features.
Bioinformatics, 30(7), 923-930.

McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky,
A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce
framework for analyzing next-generation DNA sequencing data. Genome research,
20(9), 1297-1303.

Roth, A., Ding, J., Morin, R., Crisan, A., Ha, G., Giuliany, R., ... & Marra,
M. A. (2012). JointSNVMix: a probabilistic model for accurate detection of
somatic mutations in normal/tumour paired next-generation sequencing data.
Bioinformatics, 28(7), 907-913.

Saunders, C. T., Wong, W. S., Swamy, S., Becq, J., Murray, L. J., & Cheetham,
R. K. (2012). Strelka: accurate somatic small-variant calling from sequenced
tumor–normal sample pairs. Bioinformatics, 28(14), 1811-1817.

Shen, R., & Seshan, V. E. (2016). FACETS: allele-specific copy number and
clonal heterogeneity analysis tool for high-throughput DNA sequencing. Nucleic
acids research, 44(16), e131-e131.

Wang, K., Li, M., & Hakonarson, H. (2010). ANNOVAR: functional annotation of
genetic variants from high-throughput sequencing data. Nucleic acids research,
38(16), e164-e164.

Xi, R., Lee, S., Xia, Y., Kim, T. M., & Park, P. J. (2016). Copy number
analysis of whole-genome data using BIC-seq2 and its application to detection
of cancer susceptibility variants. Nucleic acids research, gkw491.
