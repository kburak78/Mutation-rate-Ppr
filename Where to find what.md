# Where to find what 

all data in /RAID/Data/linda/all_data/

## mapped_data/ 
-> mapped and sorted bam files (*.sort.bam) and indexes (*sort.bam.bai)
-> mapped and sorted bam files + duplicates removed (*.sort.de.bam) and index (*sort.de.bam.bai) 

### mapped_data/genome_coverage/ 
-> samtools flagstat (*.flagstat)
-> bedtools gencov -bga outputfile (*_genome_coverage)
-> List of the average genome coverage per sample (average.genomecov.all)
-> Genome Size of each sample (chr.bed)
-> shell for calculating genome coverage with bedtools genomecov -bga and samtools flagstat (genome_coverage.sh) 
-> shell for calculating the coverage of each sample (genome_coverage.sh) and it's log file (genome_coverage.log) 
-> shell for calculating mean coverage from each sample (mean_cov.sh) and it's log-file (mean_cov.log). Result in average.genomecov.all) 
-> test file for mean_cov.sh (test_genome_coverage) and it's test file (test.average)



