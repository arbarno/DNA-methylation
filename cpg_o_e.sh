## ADAPTED FROM https://github.com/jldimond/Coral-CpG/blob/master/ipynb/Amil_exp_CpG_ratio.ipynb ##

## get just the headers from the positive strands and make txt file using excel (filter function) from the gff file

## filter out the positive strand mRNA transcripts from the mrna-transcripts file obtained from funannotate
# I first had to make the individual sequences be on one line each under the headers
awk '/^>/ {if (seq) print seq; print; seq=""} /^>/ {next} {seq=seq $0} END {if (seq) print seq}' Acropora_gladularis.mrna-transcripts.fa > single_line_fasta.fa
# the filter out just the positive strand mRNA transcripts
awk '{print $1}' mrna_pos_strand_headers.txt | grep -A1 -F -f - single_line_fasta.fa | sed '/-$/d' > Acropora_gladularis.mrna-pos_strand-transcripts.fa

# make the positive strand mRNA transcripts into tab format
perl -e '$count=0; $len=0; while(<>) {s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) {print "\n"} s/ |$/\t/; $count++; $_ .= "\t";} else {s/ //g; $len += length($_)} print $_;} print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n";' \
Acropora_gladularis.mrna-pos_strand-transcripts.fa > a_gland_fasta2tab

# add length column for each sequence
perl -e '$col = 2;' -e 'while (<>) { s/\r?\n//; @F = split /\t/, $_; $len = length($F[$col]); print "$_\t$len\n" } warn "\nAdded column with length of column $col for $. lines.\n\n";' a_gland_fasta2tab > tab1

# filter out just the column with sequence data (column 2) to do the CpG/C/G calculations
awk '{print $2}' tab1 > tab2

# get CpG counts
echo "CG" | awk -F\[Cc][Gg] '{print NF-1}' tab2 > CG

# get C counts
echo "C" | awk -F\[Cc] '{print NF-1}' tab2 > C

# get G counts
echo "G" | awk -F\[Gg] '{print NF-1}' tab2 > G

# combine the count data with the tabular formatted mRNA transcript fasta
paste tab1 CG C G > comb

# perform calculation per Gavery & Roberts (2010) - https://link.springer.com/article/10.1186/1471-2164-11-483#Sec8
# make sure C > 0, G > 0, and length > 1, because this will result in dividing by zero
awk '{if ($5 != 0 && $6 != 0 && $3 != 1) print $1, "\t", (($4)/($5*$6))*(($3^2)/($3-1))}' comb > ID_CpG

# sort the CpG file
sort ID_CpG > ID_CpG.sorted
