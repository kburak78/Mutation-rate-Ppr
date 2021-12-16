# Where to find what 


## Table of Contents  
[Headers](#headers)  
[Emphasis](#emphasis)  
...snip...    
<a name="headers"/>
# Headers

all data in /RAID/Data/liinda/all_data/

mapped_data/ 
-> sorted bam files (*.sort.bam) and indexes (*sort.bam.bai)
-> sorted bam files + duplicates removed (*.sort.de.bam) and index (*sort.de.bam.bai) 

mapped_data/genome_coverage/ 
-> samtools flagstat (*.flagstat)
-> bedtools gencov 
-> List of the average genome coverage per sample (average.genomecov.all)

-> shell for calculating the average coverage of each sample (genome_coverage.sh) and it's log file (genome_coverage.log) 
-> test file for mean_cov.sh (test_genome_coverage) and it's test file (test.average)

-> Genome Size of each sample (chr.bed)
