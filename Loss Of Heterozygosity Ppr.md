# Estimation of Mutation Rate in Plathynothrus peltifer
- - - -

This file encompasses all scripts of how I calculated the mutation rate in Plathynothrus peltifer.

The file is structured as follows:
1. [Sampling and Sequencing](#Sampling)
    1. [Sampling](#Sampling1)
    2. [Sequncing](#Seq)
    3. [Umi tools, Trimming and concatenating read data](#Umi)
2. [Alignment & Post-Alignment Processing](#Alignment)
3. [Variant Calling with GATK](#GATK)
4. [Detect sites with LoH](#Mutations)
    1. [Filter DP](#DP)
    2. [Filter GT](#GT)
    3. [Filter GQ & PL](#GQ)
    4. [Filter AD](#AD)
5. [Callable genome estimation](#CG)
    1. [Filter DP](#DP1)
    2. [Filter GT](#GT1)
    3. [Filter GQ & PL](#GQ1)
    4. [Filter AD](#AD1)
6. [Calculating the LoH rate](#Mutrate)

You will also find two other files in the repository: 
'R Graphs' - the scripts to visualize the results 
'Ppr Eggs' - what I did for the additional data of two eggs 

- - - -

## 1. Sampling and Sequencing <a name="Sampling"></a>
### a. Sampling
The oribatid mites were sampled early June at the “Moorpfad” in Dahlem, Germany (50.389204, 6.568780).
- M6 -> S1 (Mother 1) 
- E6.2 -> S2 (Egg from Mother 1)
- M30 -> S3 (Mother 2)
- E30.7 -> S4 (Egg from Mother 2)
- M31 -> S5 (Mother 3)
- E31.6 -> S6 (First egg from Mother 3)
- E31.2 -> S7 (Second egg from Mother 3)
- E31.5 -> S8 (Third egg from Mother 3) 
    
### b. Sequncing 

The DNA was send to be sequenced with Illumina Next Generation Sequencing (NovaSeq6000) with 2 additional kits -the 'NEBNext Ultra II FS DNA Library Prep Kit for Illumina' and 'NEBNext Multiplex Oligos for Illumina (Unique Dual Index UMI Adaptors DNA Set 1)'.
FastaQ files provided by Cologne's Center for Genomics.

The average reading lenght was 151. 

Raw reads in four folders in /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads.
In jbast_JB03_September6/ and jbast_JB03_Oktober3/ are the reads of samples S1-S6.
In jbast_JB03_Jan1/ and jbast_JB03_Februar1/ are the reads of samples S7 and S8.

Sequences provided by Cologne Center for Genomics (CCG). September ones went already through umi_tools and trimming in https://github.com/wintergoldcrest/umi_tools.

### c. Umi tools, Trimming and concatenating all Data
Trimming Illumina adapters and Umi adapters from raw reads.
To look into the script: cat /RAID/Data/linda/october/test.sh

Before running on motoko do >conda activate uni_tools and >mkdir clean_data

    
    umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin= (R2.fastq.gz-file) --read2-in= (R1.fastq.gz-file) --stdout= (R1.fastq.gz-outputfile) --read2-stdout
    umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin= (R2.fastq.gz-file) --read2-in= (R3.fastq.gz-file) --stdout= (R3.fastq.gz-outputfile) --read2-stdout
    /NVME/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --fastqc --paired --output_dir clean_data (R1.fastq.gz-outputfile) (R3.fastq.gz-outputfile)
    done
    

Concatenated previous September reads with October reads and January reads with february reads with cat. 

- - - -
## 2. Alignment & Post-Alignment Processing <a name="Alignment"></a>
Mapping trimmed reads to the reference genome.
To look into script: cat /RAID/Data/linda/all_data/mapping.sh 
Before, runnning the command: 

bwa index (reference)
conda activate uni_tools
mkdir tmp
         
    bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" (reference) (trimmed R1 output) (trimmed R3 output)> (samfile)
                #Run Aligntment with masked genome
    samtools view -bS (samfile) > (bamfile)
                #compress .sam into .bam
    samtools sort (bamfile) -o (sorted bamfile)
                #sort mapped reads according to position in the genome
    samtools index (sorted bamfile)
                #index sorted bams, creates .sort.bam.bai-file
    umi_tools dedup -I (sorted bamfile) --output-stats=deduplicated --paired -S (sorted deduplicated bamfile) --temp-dir=tmp
                #remove all duplications 
    samtools index (sorted deduplicated bamfile)
                #index deduplicated bam
    done
- - - -
## 3. Variant Calling with GATK <a name="GATK"></a>

Creating a GVCF and VCF-File with GATK HaplotypeCaller and GenotypeGVFCs.

    #software
    gatk=(GATK location)
    #data
    ref=(reference)
    #file names
    bam1=$1 bam2=$2 gvcf1=$3 gvcf2=$4 mergedgvcf=$5 vcf=$6
    #build dict for genome
    $gatk CreateSequenceDictionary -R $ref -O (outout-dict)
    #call gvcf
    $gatk HaplotypeCaller -R $ref --emit-ref-confidence GVCF -I $bam1 -O $gvcf1
    #call gvcf
    $gatk HaplotypeCaller -R $ref --emit-ref-confidence GVCF -I $bam2 -O $gvcf2
    #merge them 
    $gatk CombineGVCFs -R $ref -V $gvcf1 -V $gvcf2 -O merged_gvcf/$mergedgvc
    #detect SNPs
    $gatk GenotypeGVCFs -R $ref -V $mergedgvcf -O $vcf
    #compress
    bgzip -f $vcf
    tabix -p vcf ${vcf}.gz    
    # Filter out only SNPs from VCF 
    $gatk SelectVariants -select-type SNP -V ${vcf}.gz -O ${vcf}.snp.gz   
    # filter SNPs by parameters
    $gatk VariantFiltration -V ${vcf}.snp.gz --filter-expression "QD <2.0 || MQ <40.0 || FS >60.0 || SOR >3.0 || ReadPosRankSum < -8.0" --filter-name "PASS" -O ${vcf}.snp.f.gz
    
Executed the script with this: 

    sh VCFcalling.sh ../masked_mapping_data/153621_S1.sort.de.bam ../masked_mapping_data/153622_S2.sort.de.bam gvcf/s1.g.vcf gvcf/s2.g.vcf merged_vcf/s12.g.vcf merged_vcf/s12.vcf 
    sh VCFcalling.sh ../masked_mapping_data/153623_S3.sort.de.bam ../masked_mapping_data/153624_S4.sort.de.bam gvcf/s3.g.vcf gvcf/s4.g.vcf merged_vcf/s34.g.vcf merged_vcf/s34.vcf
    sh VCFcalling.sh ../masked_mapping_data/153625_S5.sort.de.bam ../masked_mapping_data/15366_S6.sort.de.bam gvcf/s5.g.vcf gvcf/s6.g.vcf merged_vcf/s56.g.vcf merged_vcf/s56.vcf
- - - -
## 4. Detect sites with LoH <a name="Mutations"> </a>

Following sites where kept in the VCF:
i) Filter DP: read coverage between 0,5x and 2x of the average genome coverage of all samples
ii) Filter GT: SNPs that are homozygous in the Mother and heterozygous in the Daughter
iii) Filter GQ and PL: GQ = 99 and PL difference 120 in Mother sites and 100 in Daughter sites 
iv) Filter AD:  


### i. Filter DP <a name="DP"> </a>

Calculated the average genome coverage with bedtools

Script in RAID/Data/linda/all_data/mapped_data/genome_coverage/genome_coverage.sh

    bedtools genomecov -ibam (bam-file) -bga > (output-file)

Then I calculated the average with awk and created a flagfile for each sample. (Script in RAID/Data/linda/all_data/mapped_data/genome_coverage/mean_cov.sh)

Average genome coverage of each sample in average-genomecov.all

    Average coverage of S1= 66.297
    Average coverage of S2= 65.3743
    Average coverage of S3= 66.4849
    Average coverage of S4= 71.9713
    Average coverage of S5= 64.3218
    Average coverage of S6= 65.4336
    
Average genome coverage for all samples is 66.64715. 

### ii. Filter GT <a name="GT"> </a>
### iii. Filter GQ and PL <a name="GQ"> </a>
### iv. Filter AD <a name="AD"> </a>

- - - -
## 5. Callable Genome estimation <a name="CG"> </a>

adding header to hf and gcov filtered data:

cat /RAID/Data/linda/Mother_Egg_Pairs_Ppr/vcf/merged_gvcf/filter/SNPs_after_gcov_filter/addheader.sh

### i. Filter DP <a name="DP1"> </a>
### ii. Filter GT <a name="GT1"> </a>
### iii. Filter GQ and PL <a name="GQ1"> </a>
### iv. Filter AD <a name="AD1"> </a>

- - - -
## 6. Calculating the LoH rate <a name="Mutrate"> </a> 

The Mutation rate can be calculated as: 

Mutation counts / (2 x callable sites) = Mutation rate 

With the first calculations something between 1-2x10^-6 comes out which would be higher then for example spider mites. If Ppr has a lower mutation rate we would expect something under around 1*10^-9. 

### Identifying Coding and Non-coding regions with existing gtf-file

extract coding regions from gtf-file.

    awk -v OFS='\t' '$3=="gene"{print$1,$4,$5}' /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf > coding_area_Ppr
    
    coding_area.sh

/RAID/Data/linda/all_data/mapped_data/genome_coverage/coding_area_gcov 
