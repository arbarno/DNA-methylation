# DNA-methylation
Analyzing DNA methylation data from nanopore sequencing reads involves several steps to process, align, and extract methylation information from the raw data. Below is a general outline of a DNA methylation analysis pipeline for nanopore reads data:

(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02510-z)  

- Another useful tuterial: https://timkahlke.github.io/LongRead_tutorials/

## Quality Control:

Assess the quality of the raw nanopore reads using NanoPlot to identify any issues that may need attention.

https://github.com/wdecoster/NanoPlot

````bash
NanoPlot -t 16 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots dot --legacy hex
````
  
## Adapter Removal using PoreChop:

````bash
porechop –-input reads.fastq -o reads_porechopped.fastq --discard_middle
````
## Kmer profiling:
Usually, this is for short-reads or high-accurate long reads as "HiFi technology" but we could try
### genomescope2.0
https://github.com/tbenavi1/genomescope2.0

I ran the kmer profiling using KMC, then inputted the graph into the Genome-scope website: 

http://genomescope.org/genomescope2.0/

### Smudgeplot
I tried smudgeplot, but it kept crashing (memory).

https://github.com/KamilSJaron/smudgeplot

````bash
kmc -k21 -t50 -m128 -ci1 -cs10000 reads_porechopped.fastq kmcdb output/
kmc_tools transform kmcdb histogram kmcdb_21.histo -cx10000
````

## Genome assembly
Tested using using Flye (with the nano_corr, nano_hq, and nano_raw setting), NECAT, and CANU

The best assembly was Flye nano_corr, so this was used for all samples.

### 1a. Flye (nano_corr)

````bash
flye --nano-corr [input.fastq] -g 421m -o [output_directory] --scaffold -t 32
````
### 1b. Flye (nano_hq)

````bash
flye --nano-hq [input.fastq] -g 421m -o [output_directory] --scaffold -t 32
````
### 1c. Flye (nano_raw)

````bash
flye --nano-raw [input.fastq] -g 421m -o [output_directory] --scaffold -t 32
````
### 2.NECAT

````bash
necat.pl correct acro_config.txt
necat.pl assemble acro_config.txt
necat.pl bridge acro_config.txt
````
- Configuration file for NECAT

````bash
PROJECT=acro_necat
ONT_READ_LIST=/ibex/project/c2208/nanopore/tests/read_list.txt
GENOME_SIZE=421000000
THREADS=1
MIN_READ_LENGTH=1000
PREP_OUTPUT_COVERAGE=20
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=10
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true
````
### 3. CANU

````bash
canu -p [output_prefix] -d [output_directory] genomeSize=421m -nanopore -trimmed -correct -assemble [input.fastq]
````
## Medaka for polishing the contigs

### Medaka

````bash
mini_align -r ${f} -i ../porechop/${name}_porechopped.fastq.gz -m -p ../medaka/${name}_calls_to_draft -t 32

samtools view -H ${f} | grep "^@SQ" | cut -f 2,3 | sed 's/SN://; s/\tLN//' > ${name}_regions.txt
sed -i 's/:/:0-/g' ${name}_regions.txt
split -d -n l/20 ${name}_regions.txt ${name}_split_regions- --additional-suffix .txt

xargs -n 1 < ${name}_split_regions-00.txt | xargs medaka consensus ${f} /ibex/project/c2208/nanopore/medaka/temp/${name}_split_regions-00.hdf --model r1041_e82_400bps_sup_v4.2.0 --batch 200 --threads 8 --region

medaka stitch temp/${name}*.hdf ../flye/${name}_assembly.fasta results/${name}.polished.assembly.fasta
````

## Scaffolding
LINKS was used first for scaffolding, then the gaps were filled from this output using TGS-GapCloser

### LINKS

````bash
LINKS -f ${f} -s ${name}_porechop.fof -k 21 -t 10 -b ../../links/${name}
````
### TGS-GapCloser

````bash
tgsgapcloser --scaff ${f} \
  --reads ../porechop/fasta/${name}_porechopped.fasta \
	--output ../tgs_gapcloser/${name} \
	--racon /ibex/sw/rl9c/racon/1.5.0/rl9_conda3/env/bin/racon \
	--thread 32
````

### Redundans (WAS NOT USED!)

````bash
redundans.py -f <assembly.fasta> -o <output_directory>

-f <assembly.fasta>: Replace <assembly.fasta> with the path to your input nanopore assembly in FASTA format.

-o <output_directory>: Specify the directory where Redundans should save the output files, including the improved assembly.
````

## Genome assembly assessment 
### BlobToolKit (BTK)
- You can find information about the tool in this course https://www.futurelearn.com/courses/eukaryotic-genome-assembly-how-to-use-blobtoolkit-for-quality-assessment
- Here is a handbook for installation and running the tool https://github.com/blobtoolkit/tutorials/tree/main/futurelearn 

## Read Alignment:

Map the nanopore reads to a reference genome using a suitable aligner like Minimap2 or GraphMap.

````bash
minimap2 -ax map-ont -t 40 Reference.fasta reads.fastq |samtools sort -@40 -O BAM -o mapped.bam -
````

## DNA Methylation Calling:

Detect methylated bases in the aligned reads. For nanopore sequencing data, methylation is typically detected through the modification of the electrical signal by DNA methyltransferases.
You can use tools like Tombo, Nanopolish, or DeepSignal to perform methylation calling.

## Methylation Data Analysis:

Convert the methylation calls into a more interpretable format (e.g., BED or BEDGraph) for downstream analysis.
Perform differential methylation analysis, if applicable, to compare methylation patterns between different samples or conditions.

## Visualization:

Visualize the methylation patterns using genome browsers or specialized tools like Integrative Genomics Viewer (IGV) or UCSC Genome Browser.

## Annotation:

Annotate the differentially methylated regions (DMRs) with genomic features (e.g., genes, CpG islands) to gain biological insights.

## Functional Enrichment Analysis:

Perform functional enrichment analysis on genes associated with DMRs to understand the potential functional implications of DNA methylation changes.

## Integration with Other Data:

Integrate DNA methylation data with other omics data, such as gene expression or chromatin accessibility data, for a more comprehensive understanding of gene regulation.

## Validation:

If possible, validate the DNA methylation changes using an independent method, such as bisulfite sequencing.

Please note that specific tools and parameters used in the pipeline might vary depending on the characteristics of the nanopore data and the research question you are addressing. It is crucial to refer to the respective software documentation and literature for the most up-to-date protocols and guidelines for DNA methylation analysis with nanopore sequencing data.
