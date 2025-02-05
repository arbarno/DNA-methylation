# DNA-methylation
Analyzing DNA methylation data from nanopore sequencing reads involves several steps to process, align, and extract methylation information from the raw data.
Below is a general outline of a DNA methylation analysis pipeline for nanopore reads data:

(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02510-z)  

- Another useful tuterial: https://timkahlke.github.io/LongRead_tutorials/

## Quality Control:

### NanoPlot
https://github.com/wdecoster/NanoPlot

NanoPlot v1.42.0<br/>
input files: ONT long reads

````bash
NanoPlot -t 16 --fastq long_reads.fastq.gz --maxlength 40000 --plots dot --legacy hex
````
  
## Kmer profiling:
Usually, this is for short-reads or high-accurate long reads as "HiFi technology" but we could try

### genomescope2.0
https://github.com/tbenavi1/genomescope2.0 
http://genomescope.org/genomescope2.0/
### Smudgeplot
https://github.com/KamilSJaron/smudgeplot

## Genome assembly

### MaSuRCA
https://github.com/alekseyzimin/masurca

MaSuRCA v4.1.0<br/>
input files: Ultima short reads (SE 300bp) + ONT long reads<br/>
Species: _Acropora glandularis_

````bash
masurca -i short_reads.fastq.gz -r long_reads.fastq.gz -t 32
````

## Reference-base scaffolding

### CSAR
https://github.com/ablab-nthu/CSAR

CSAR v1.1.1<br/>
input files: MaSuRCA assembly + reference genome (in this case _Acropora millepora_ because it was the closest relative from BLAST)<br/>
_A. millepora_ reference: GCA_013753865.1

````bash
conda activate csar

php csar.php -t masurca_assembly.fasta -r GCA_013753865.1_Amil_v2.1_genomic.fna --nuc -o output_folder

conda deactivate
````

## Polishing the assembly

### Racon
https://github.com/isovic/racon

minimap2 v2.24<br/>
racon v1.5.0<br/>
input files: CSAR scaffolds + ONT long reads

````bash
minimap2 -t 32 -ax map-ont csar_scaffolds.fna long_reads.fastq.gz -o mapped_long_reads.sam

racon -u -m 3 -x -5 -g -4 -w 500 -t 32 long_reads.fastq.gz mapped_long_reads.sam csar_scaffolds.fna > racon_polished_scaffold.fasta
````

### Pilon
https://github.com/broadinstitute/pilon/wiki

java v19.0.1<br/>
minimap2 v2.24<br/>
samtools v1.16.1<br/>
pilon v1.24<br/>
input files: racon-polished scaffolds + Ultima short reads

````bash
minimap2 -d  racon_polished_scaffold.mmi racon_polished_scaffold.fasta
minimap2 -ax sr racon_polished_scaffold.mmi short_reads.fastq.gz > mapped_short_reads.sam
samtools view -bS mapped_short_reads.sam | samtools sort -o mapped_short_reads.bam
samtools index mapped_short_reads.bam

java -Xms64G -Xmx128G -jar $PILON \
	--genome racon_polished_scaffold.fasta \
	--unpaired mapped_short_reads.bam \
	--output output_prefix \
	--changes \
	--threads 48
````

## Creating homozygous genome assembly

### Redundans
https://github.com/Gabaldonlab/redundans

redundans v2.0.1
input files: pilon/racon-polished scaffolds + Ultima short reads + ONT long reads

````bash
redundans.py -v -i short_reads.fastq.gz \
	-l long_reads.fastq.gz \
	-f pilon_racon_polished_scaffold.fasta \
	-o output_folder \
	-t 32 \
	--log log.txt
````

## Genome assembly assessment 
### BlobToolKit (BTK)
https://www.futurelearn.com/courses/eukaryotic-genome-assembly-how-to-use-blobtoolkit-for-quality-assessment
https://github.com/blobtoolkit/tutorials/tree/main/futurelearn

Create the blobtools input files (blastn, diamond, busco, and coverage)

#### blastn v2.13.0
````bash
blastn -db NCBI/nt \
	-query redundans_reduced_scaffolds.fa \
	-outfmt "6 qseqid staxids bitscore std" \
	-max_target_seqs 10 \
	-max_hsps 1 \
	-evalue 1e-25 \
	-num_threads 30 \
	-out scaffolds_nt.blastn
````
#### diamond v2.1.6
````bash
diamond blastx --query redundans_reduced_scaffolds.fa \
	--db uniprot/reference_proteomes.dmnd \
	--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
	--max-target-seqs 1 \
	--evalue 1e-25 \
	--fast \
	--verbose \
	--threads 48 > scaffolds_uniprot_diamond.out
````
#### busco v5.7.1
````bash
busco -i redundans_reduced_scaffolds.fa \
	-l busco/eukaryota_odb10 \
	-o scaffolds_output_busco_folder \
	-m genome \
	--cpu 30
````
#### coverage (minimap2 v2.24, samtools v1.16.1)
````bash
minimap2 -ax sr -t 30 redundans_reduced_scaffolds.fa \
	short_reads.fastq.gz | samtools sort -@30 -O BAM -o scaffolds_coverage.bam
samtools index -c scaffolds_coverage.bam
````

Combine files into a single blobtoolkit folder for visualization

#### blobtoolkit v4.4.0
````bash
conda activate btk

blobtools create --fasta redundans_reduced_scaffolds.fa \
	--meta acropora.yaml \
	--taxdump ncbi_taxdump/ \
	btk_folder
blobtools add --hits scaffolds_nt.blastn \
	--hits scaffolds_uniprot_diamond.out \
	--taxrule bestsumorder \
	--taxdump ncbi_taxdump/ \
	btk_folder
blobtools add --cov scaffolds_coverage.bam \
	--threads 30 \
	btk_folder
blobtools add --busco scaffolds_output_busco_folder/run_eukaryota_odb10/full_table.tsv \
	btk_folder

conda deactivate
````

Filter the btk files to keep only the contigs that hit to phylum = Cnidaria, genus = Acropora, coverage = 0.01

````bash
conda activate btk

blobtools filter \
	--param bestsumorder_phylum--Inv=Cnidaria \
	--param bestsumorder_genus--Inv=Acropora \
	--param scaffolds_coverage--Min=0.01 \
	--fasta redundans_reduced_scaffolds.fa \
	--summary filter_summary \
	btk_folder

conda deactivate
````

## Read Alignment:

### Adapter Removal using PoreChop:

````bash
porechop –-input reads.fastq -o reads_porechopped.fastq --discard_middle
````

Map the nanopore reads to the consensus genomes using Minimap2, then filter the original bam files that contain the methylation data.

````bash
minimap2 -ax map-ont -t 40 ${f} ../raw/${name}.fastq.gz | samtools sort -@40 -O BAM -o ../mapped_bam/${name}.mapped.bam
samtools view -b -F 4 ../mapped_bam/${name}.mapped.bam > ../mapped_bam/${name}.mapped_filtered.bam
samtools view ../mapped_bam/${name}.mapped_filtered.bam | cut -f 1 > ../mapped_bam/${name}.mapped_read_names.txt
samtools view -h ../bam/${name}.bam | grep -F -w -f ../mapped_bam/${name}.mapped_read_names.txt | samtools view -Sb - > ../filtered_bam/${name}.filtered.bam
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
