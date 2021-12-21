# Mutation-rate-Ppr

Sequences provided by CCG in September and October. September ones went already through umi_tools and trimming in https://github.com/wintergoldcrest/umi_tools. 

## 1. Umi tools and Trimming of October Data
Trimming Illumina adapters and Umi adapters from raw reads.
To look into the script: cat /RAID/Data/linda/october/test.sh

    for i in A006200184_153621_S1 \
    A006200184_153622_S2 \
    A006200184_153623_S3 \
    A006200184_153624_S4 \
    A006200184_153625_S5 \
    A006200184_153626_S6
    do
    #umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_Oktober3/${i}_L002_R2_001.fastq.gz --read2-in=/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_Oktober3/${i}_L002_R1_001.fastq.gz --stdout=${i}_add_barcode_R1.fastq.gz --read2-stdout
    #umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_Oktober3/${i}_L002_R2_001.fastq.gz --read2-in=/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_Oktober3/${i}_L002_R3_001.fastq.gz --stdout=${i}_add_barcode_R3.fastq.gz --read2-stdout
    /NVME/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --fastqc --paired --output_dir clean_data ${i}_add_barcode_R1.fastq.gz ${i}_add_barcode_R3.fastq.gz
    done
    

## 2. Concatenated previous September data and new coverage data from October

Result in /home/linda/Data/all_data/

    i in 153621_S1 \
    153622_S2 \
    153623_S3 \
    153624_S4 \
    153625_S5 \
    153626_S6
    do
    cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/A006200178_${i}*add_barcode_R1_val_1.fq.gz ../october/clean_data/A006200184*${i}*add_barcode_R1_val_1.fq.gz > $i.R1_val_1.fq.gz
    cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/A006200178*${i}*add_barcode_R3_val_2.fq.gz ../october/clean_data/A006200184*${i}_add_barcode_R3_val_2.fq.gz > $i.R3_val_2.fq.gz
    done



## 3. Mapping 
Mapping trimmed reads to the reference genome.
To look into script: cat /RAID/Data/linda/all_data/mapping.sh 

bwa index ref.fe
conda activate umi_tools
bash nohup mapping.sh > mapping.log &

    for i in 153621_S1 \
    153622_S2 \
    153623_S3 \
    153624_S4 \
    153625_S5 \
    153626_S6
    do
    bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" /RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa ${i$
    samtools view -bS mapped_data//${i}.sam > mapped_data//${i}.bam
    samtools sort mapped_data/${i}.bam -o mapped_data/${i}.sort.bam
    samtools index mapped_data/${i}.sort.bam
    umi_tools dedup -I mapped_data/${i}.sort.bam --output-stats=deduplicated --paired -S mapped_data/${i}.sort.de.bam --temp-dir=tmp 
    done


## 4. SNP and Indel calling 

### with GATK 

Calling Variants with HaplotypeCaller and GenotypeGVCF.

Script in /RAID/Data/linda/all_data/sim_SNP_calling_w_GATK.sh

    #software
    gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk
    #reference_data
    ref=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa
    bam=$1
    out1=$2
    out2=$3
    #build dict for genome, this time not because I already had build one
    #/NVME/Software/popgen/gatk-4.1.9.0/gatk CreateSequenceDictionary -R Ppr.FINAL.fa -O Ppr.FINAL.dict
    #call gvcf
    $gatk HaplotypeCaller \
    -R $ref \
    --emit-ref-confidence GVCF \
    -I $bam \
    -O $out1

    #detect SNPs
    $gatk GenotypeGVCFs \
    -R $ref \
    -V $out1 \
    -O $out2

    #compress
    bgzip -f $out2
    tabix -p vcf $out2.gz
    
To bash everything I used 'echo'. Script in /RAID/Data/linda/all_data/mapping.sh

    for i in 153621_S1 \
    153622_S2 \
    153623_S3 \
    153624_S4 \
    153625_S5 \
    153626_S6
    do

    echo "sh sim_SNP_calling_w_GATK.sh mapped_data/${i}.sort.de.bam vcf/${i}.g.vcf vcf/${i}.vcf " > shell/${i}.sh
    nohup sh shell/${i}.sh/) > shell/$i.log &
    done

Commands and log in /RAID/Data/linda/all_data/shell, compressed results in /RAID/Data/linda/all_data/vcf.
  

## 5. Hardfiltering the Variants

For Filters I used this, which was also proposed for 
QD <2.0 || MQ <40.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0"
QD <2.0 || FS >100.0 || SOR >5.0 || ReadPosRankSum < -8.0


Script in /RAID/Data/linda/all_data/vcf/filter.sh

    ###software
    gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk
    ###data
    ref=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa
    vcf=$1
    snpvcf=$2
    indelvcf=$3
    filterSNP=$4
    filterINDEL=$5
    finalvcf=$6
    ###SelectVariants SNP
    $gatk SelectVariants \
    -select-type SNP \
    -V $vcf \
    -O $snpvcf
    ###SelectVariants INDEL
    $gatk SelectVariants \
    -select-type INDEL \
    -V $vcf \
    -O $indelvcf
    ###filter SNP
    $gatk VariantFiltration \
    -V $snpvcf \
    --filter-expression "QD <2.0 || MQ <40.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O $filterSNP
    ###filter INDEL
    $gatk VariantFiltration \
    -V $indelvcf \
    --filter-expression "QD <2.0 || FS >100.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O $filterINDEL
    ###merge SNP INDEL
    $gatk MergeVcfs \
    -I $filterSNP \
    -I $filterINDEL \
    -O $finalvcf
    ###delete temp
    rm -f $snpvcf $indelvcf $filterSNP $filterINDEL
        
        
To run everything simultaneously and I used 'echo'.
Script in /RAID/Data/linda/all_data/vcf/bash_filter_sim.sh

    for i in A006200178_153621_S1 \
    A006200178_153622_S2 \
    A006200178_153623_S3 \
    A006200178_153624_S4 \
    A006200178_153625_S5 \
    A006200178_153626_S6
    do
    echo "sh filter.sh ${i}.vcf.gz ${i}.snp.vcf.gz ${i}.indel.vcf.gz ${i}.f.snp.vcf.gz ${i}.f.indel.vcf.gz ${i}.f.vcf.gz" >filter_shell/${i}.f.sh
    nohup sh filter_shell/${i}.f.sh > filter_shell/${i}.f.sh &
    done

## 6. Calculated the average genome coverage with bedtools

Sorted BAM-Files in /RAID/Data/linda/all-files/mapped_data
Script in cat RAID/Data/linda/all_data/mapped_data/genome_coverage/genome_coverage.sh


    for i in 153621_S1 \
    153622_S2 \
    153623_S3 \
    153624_S4 \
    153625_S5 \
    153626_S6
    do
    bedtools genomecov -ibam ../${i}.sort.de.bam -bga > ${i}_genome_coverage
    done

Then I calculated the average with awk.
Script in RAID/Data/linda/all_data/mapped_data/genome_coverage/mean_cov.sh 
And created a samtools flagfile for each sample.

Average genome coverage of each sample in average-genomecov.all

    Average coverage of S1= 66.297
    Average coverage of S2= 65.3743
    Average coverage of S3= 66.4849
    Average coverage of S4= 71.9713
    Average coverage of S5= 64.3218
    Average coverage of S6= 65.4336

## 7. Identifying Coding and Non-coding regions with existing gtf-file

extract coding regions from gtf-file.

    awk -v OFS='\t' '$3=="gene"{print$1,$4,$5}' /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf > coding_area_Ppr
    
    
coding_area.sh

/RAID/Data/linda/all_data/mapped_data/genome_coverage/coding_area_gcov 

## 8. Selected Biallelic SNPs
In /RAID/Data/linda/all_data/vcf/f.biallelic

First I used GATK SelectVariants to omit any multiallelic SNPs ( cat /RAID/Data/linda/all_data/vcf/f.biallelic/f.biallelic.sh ).

    for i in 153621_S1 \
    153622_S2 \
    153623_S3 \
    153624_S4 \
    153625_S5 \
    153626_S6
    do
    /NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
    -V ../${i}.f.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    -O ${i}.f.bi.vcf.gz
    done

And then selected all SNPs to omit the Indels ( cat /RAID/Data/linda/all_data/vcf/f.biallelic/select_SNPs.sh ).

    for i in 153621_S1 \
    153622_S2 \
    153623_S3 \
    153624_S4 \
    153625_S5 \
    153626_S6
    do
    /NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
    -select-type SNP \
    -V ${i}.f.bi.vcf.gz \
    -O ${i}.f.bi.snp.vcf.gz
    done

