# Analyzing coral DNA methylation in parents and gametes using ONT
Analyzing DNA methylation data from nanopore sequencing reads involves several steps to process, align, and extract methylation information from the raw data.
Below is a general outline of a DNA methylation analysis pipeline for nanopore reads data:

(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02510-z)  

- Another useful tuterial: https://timkahlke.github.io/LongRead_tutorials/

# Assembling an _Acropora_ reference genomes

We must first make a reference assembly using one of the samples to perform downstream methylation analysis.

## Quality Control:

### NanoPlot
https://github.com/wdecoster/NanoPlot

NanoPlot v1.42.0  
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

MaSuRCA v4.1.0  
input files: Ultima short reads (SE 300bp) + ONT long reads  
Species: _Acropora glandularis_

````bash
masurca -i short_reads.fastq.gz -r long_reads.fastq.gz -t 32
````

## Reference-base scaffolding

### CSAR
https://github.com/ablab-nthu/CSAR

CSAR v1.1.1  
input files: MaSuRCA assembly + reference genome (in this case _Acropora millepora_ because it was the closest relative from BLAST)  
_A. millepora_ reference: GCA_013753865.1

````bash
conda activate csar

php csar.php -t masurca_assembly.fasta -r GCA_013753865.1_Amil_v2.1_genomic.fna --nuc -o output_folder

conda deactivate
````

## Polishing the assembly

### Racon
https://github.com/isovic/racon

minimap2 v2.24  
racon v1.5.0  
input files: CSAR scaffolds + ONT long reads

````bash
minimap2 -t 32 -ax map-ont csar_scaffolds.fna long_reads.fastq.gz -o mapped_long_reads.sam

racon -u -m 3 -x -5 -g -4 -w 500 -t 32 long_reads.fastq.gz mapped_long_reads.sam csar_scaffolds.fna > racon_polished_scaffold.fasta
````

### Pilon
https://github.com/broadinstitute/pilon/wiki

java v19.0.1  
minimap2 v2.24  
samtools v1.16.1  
pilon v1.24  
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

## Gene prediction

Using EDTA and RepeatMasker, DNA repeat regions were first identified and softmasked.  
Then using BRAKER to generate full gene structure annotations.

### The Extensive de novo TE Annotator (EDTA)
https://github.com/oushujun/EDTA  

edta v2.0.1  
input files: btk-filtered _A. gladularis_ assembly genome; related genomes (in this case, 3 NCBI RefSeq _Acropora_ genomes - GCF_000222465.1, GCF_013753865.1, GCF_036669905.1)

````bash
EDTA.pl --genome genome.fa --sensitive 1 --threads 40
````

Concatenate the outputs into a single repeat database and get statistics 

````bash
cat *.TElib.fa > acropora_repeat_library.fa
seqkit rmdup -s -i acropora_repeat_library.fa > acropora_repeat_library_uniq.fa
grep '>' acropora_repeat_library_uniq.fa | sed -r 's/.+#//' | sed -r 's/\s+.+//' | sort | uniq -c > acropora_repeat_stats.txt
````

### RepeatMasker
https://github.com/Dfam-consortium/TETools  

repeatmasker v4.1.7  
input files: btk-filtered _A. gladularis_ assembly genome + _Acropora_ repeat library

````bash
RepeatMasker genome.fa -lib acropora_repeat_library.fa -pa 1 -norna -xsmall -dir repeat_masked
````

### BRAKER
https://github.com/Gaius-Augustus/BRAKER  

singularity v3.9.7  
braker v3.0.8  
input files: btk-filtered _A. gladularis_ assembly genome + OrthoDB12 (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Metazoa.fa.gz) concatenated with related proteomes (NCBI RefSeq _Acropora_ genomes - GCF_000222465.1, GCF_013753865.1, GCF_036669905.1)

````bash
singularity exec braker3.sif braker.pl --species=A_glandularis --threads 40 \
	--genome=repeat_masked/genome.fa.masked \
	--prot_seq=orthodb12_metazoa_ncbi_acropora_proteins.fa
````

### tRNA prediction

````bash
tRNAscan-SE -E -I -H --detail --thread 50 -o trnascan-se.out -f trnascan-se.tbl -m trnascan-se.log  repeat_masked/genome.fa.masked

EukHighConfidenceFilter -i trnascan-se.out -s trnascan-se.tbl -o eukconf -p filt
````

- Implement the tRNA-Scan results.

````bash
#covert tRNA to gff after removing non-high confident
perl convert_tRNAScanSE_to_gff3.pl --input=trna_data.txt > trna_annotation.gff
#merge gff files
agat_sp_merge_annotations.pl --gff braker.gff --gff trna_annotation.gff --out merged.gff
#export protein sequences to proceed with functional annotation
gffread merged.gff -g genome.fa -y Acropora.braker.prot.fasta
````
# Functional annotation

- transmembrane topology and signal peptide predictor, which requires install phobius
````bash
./phobius/phobius.pl -short Acropora.braker.prot.fasta > phobius.results.txt
````
- hint: remove the '*' at the end of each protein sequence at Acropora.braker.prot.fasta

- IterProScan (i'm using the funnannotate command)

````bash
funannotate iprscan -i Acropora.braker.prot.fasta -m docker -c 50
````
- eggnog-mapper
````bash

emapper.py --cpu 50 -m mmseqs --data_dir funannotate_DB  -i Acropora.braker.prot.fasta -o Acropora_eggnog
````

- Implement annotation using funannotate

````bash
funannotate annotate --gff merged.gff --fasta genome.fa -s "YOU SPECIES NAME" --busco_db  metazoa --eggnog Acropora_eggnog.emapper.annotations --iprscan Acropora_iprscan.xml --cpus 100 -o Acropora_anno
````

# DNA methylation analysis
Now that the reference genome is assembled, we can analyze the DNA methylation of individual samples

## Read Alignment:

### Map the nanopore bam files containing methylation data to the consensus genomes using Dorado.  

dorado v0.9.1 
samtools v1.16.1  
input files: _A. gladularis_ assembly genome + raw ONT .bam files (containing methylation information)

````bash
dorado aligner aglandularis_ref_genome.fa reads.bam --output-dir aligned_bam
````

### Determining _Acropora_ total methylation per sample

modkit v0.2.2  
input files: filtered .bam file that was mapped to the reference genome assembly

````bash
modkit summary --no-sampling --threads 32 aligned_reads.bam > modkit_summary.txt
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
