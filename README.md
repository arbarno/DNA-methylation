# DNA-methylation
Analyzing DNA methylation data from nanopore sequencing reads involves several steps to process, align, and extract methylation information from the raw data. Below is a general outline of a DNA methylation analysis pipeline for nanopore reads data:

(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02510-z)  

- Another useful tuterial: https://timkahlke.github.io/LongRead_tutorials/


  
## Adapter Removal using PoreChop:

````bash
porechop –-input reads.fastq -o reads_porechopped.fastq --discard_middle
````
## Quality Control:

Assess the quality of the raw nanopore reads using tools like FastQC or NanoPlot to identify any issues that may need attention.

https://github.com/wdecoster/NanoPlot

````bash
NanoPlot -t 2 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots dot --legacy hex
````
- Not very informative but we can give it a try

## Genome assembly
using Flye, CANU, and MaSuRCA .....

### 1. Flye

````bash
flye --nano-hq -g 0.421g --input [input.fastq] --out-dir [output_directory] --scaffold -t 50
````
### 2. CANU

````bash
canu -p [output_prefix] -d [output_directory] genomeSize=0.421g stopOnLowCoverage=5 -nanopore-raw [input.fastq]
````
### 3.MaSuRCA

````bash
runCanu.sh nanopore-[read_type] [config_file]
````
- Configuration file for MaSuRCA

````bash
# Configuration file for MaSuRCA

DATA
  PE = 
  JUMP = 
  OTHER = nanopore-[read_type] raw_reads.fastq
  # Add additional libraries as needed

PARAMETERS
  # Specify assembly parameters here, such as genome size estimate, k-mer size, etc.

````
## Kmer profiling
Usually, this is for short-reads or high-accurate long reads as "HiFi technology" but we could try
### genomescope2.0
https://github.com/tbenavi1/genomescope2.0

````bash
ls ../*.fastq > FILES
./bin/kmc -k21 -t50 -m64 -ci1 -cs10000 @FILES reads tmp/
./bin/kmc_tools transform reads histogram reads.histo -cx10000
genomescope.R -i ../reads.histo -o output_dir -k 21
````

### Smudgeplot

https://github.com/KamilSJaron/smudgeplot


## Scaffolding
Use any of the following:

### LINKS

````bash
links -f <contigs.fasta> -s <scaffolds.fasta> -k <k-mer_size> -l <library_name> -d <library_mean> -o <output_directory>

- -f <contigs.fasta>: This should be replaced with the path to your input contig assembly in FASTA format.

- -s <scaffolds.fasta>: Specify the path where you want the scaffolded output to be saved in FASTA format.

- -k <k-mer_size>: Set the k-mer size used for building the contig graph. The choice of k-mer size depends on your data, but typical values range from 17 to 21.

- -l <library_name>: Assign a name or identifier to your sequencing library.

- -d <library_mean>: Set the mean insert size of your library. This parameter depends on the specific library preparation method and data you have. You may need to calculate this value based on your data or refer to the documentation of your sequencing library.

- -o <output_directory>: Specify the directory where LINKS should save the output files.
````
### TGS-GapCloser

````bash
TGS-GapCloser -l <reads.fastq> -s <input_assembly.fasta> -o <output_directory>

- -l <reads.fastq>: Specify the path to the nanopore reads in FASTQ format. These reads will be used to close gaps and improve the existing assembly.

- -s <input_assembly.fasta>: Provide the path to the existing nanopore assembly that you want to improve.

- -o <output_directory>: Specify the directory where TGS-GapCloser should save the improved assembly and other output files.
````

### Redundans

````bash
redundans.py -f <assembly.fasta> -o <output_directory>

-f <assembly.fasta>: Replace <assembly.fasta> with the path to your input nanopore assembly in FASTA format.

-o <output_directory>: Specify the directory where Redundans should save the output files, including the improved assembly.
````


## Genome assembly assessment 
-using BTK,.....
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
