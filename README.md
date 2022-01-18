# Mutation-rate-Ppr

Sequences provided by CCG in September and October. September ones went already through umi_tools and trimming in https://github.com/wintergoldcrest/umi_tools. 

## Samples

FastaQ files provided by Cologne's Center for Genomics.

S1 -> M6 
S2 -> E6.2
S3 -> M30
S4 -> E30.7
S5 -> M31 
S6 -> E31.6

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



## 3.1 Calculated the average genome coverage with bedtools

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

## 3.2 Identifying Coding and Non-coding regions with existing gtf-file

extract coding regions from gtf-file.

    awk -v OFS='\t' '$3=="gene"{print$1,$4,$5}' /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf > coding_area_Ppr
    
    
coding_area.sh

/RAID/Data/linda/all_data/mapped_data/genome_coverage/coding_area_gcov 

## 4. Merging Mother and Daughter & SNP Calling

Combining Mother-daughter pairs with CombineGVCFs. Converting GVCF to VCF with GenotypeGVCFs. Filtering out SNPs with SelectVariant and filtering with VariantFiltration. And then omitting multiallelic SNPs. The following code is an example for the first mother and daughter pair.  

Script in: /RAID/Data/linda/all_data/vcf/s12.filter.sh

    gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk
    ref=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa

    $gatk CombineGVCFs \
        -R $ref \
        -V 153621_S1.g.vcf \
        -V 153622_S2.g.vcf \
        -O merged_gvcf/s12.g.vcf

    $gatk GenotypeGVCFs \
            -R $ref \
            -V merged_gvcf/s12.g.vcf \
            -O merged_gvcf/s12.vcf.gz

    # Filter out only SNPs from VCF 
    $gatk SelectVariants \
    -select-type SNP \
    -V merged_gvcf/s12.vcf.gz \
    -O merged_gvcf/s12.snp.vcf.gz

    # filter SNPs by parameters
    $gatk VariantFiltration \
    -V merged_gvcf/s12.snp.vcf.gz \
    --filter-expression "QD <2.0 || MQ <40.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O merged_gvcf/s12.snp.f.vcf.gz

    # filter for biallelic snps
    $gatk SelectVariants \
        -V merged_gvcf/s12.snp.f.vcf.gz \
        --restrict-alleles-to BIALLELIC \
        -O merged_gvcf/s12.all.vcf.gz



## 5. Filtering SNPs
Removed TE-Area: 

    for i in s12 s34 s56
    do
    #bedtools intersect -a ../$i.all.vcf.gz -b /home/shangao/Data/EDTA/Ppr/Ppr_instagrall/Ppr_instagrall.polished.fa.mod.EDTA.TEanno.gff3 -v > $i.all.removedTE.vcf
    python /home/shangao/script/python/vcf_filter-same.3.py -s $i.all.removedTE.vcf -o $i.all.removedTE_fDP_hap.vcf
    done

Filtered out 
i) heterozygous mother SNPs.
ii) all homozygous mother SNPs which had an allele depth under 60 the respective allele 
iii) all homozygous daugther SNPs where either the reference or alternate allele depth is under 60
iv) GQ under 60. Left were only Homozygous to Heterozygous SNPs with an allele depth over 60 for the supported allele and GQ over 60. 

cat /home/shangao/script/python/vcf_filter-same.3.py

    import sys
    import getopt
    import os
    def usage():
        print('''Useage: python script.py [option] [parameter]
        -s/input_file           input the lift file 
        -t/temp                 blasr_result
        -l/length                novel_seq
        -o/--output              the output results file
        -h/--help                show possible options''')
    #######################default
    opts, args = getopt.getopt(sys.argv[1:], "hs:t:o:l:",["help","sequence_file=","temp=","length","output="])  
    for op, value in opts:
        if op == "-s" or op=="--sequence_file":
            sequence_file = value
        elif op == "-o" or op =="--output": 
            output = value
        elif op == "-l" or op =="--length": 
            length = value
        elif op == "-t" or op =="--temp": 
            temp = value
        elif op == "-h" or op == "--help":
            usage()
            sys.exit(1)
    f1=open(sequence_file)
    #f2=open(temp)
    #f4=open(length,'w')
    f3=open(output,'w')
    total={}
    for l in f1.readlines():
            i=l.strip().split()
            if i[9].split(':')[0] == '0/0' or i[9].split(':')[0] == '0|0' or i[9].split(':')[0] == '1/1' or i[9].split(':')[0] == '1|1':
                if i[9].split(':')[0] == '0/0' or i[9].split(':')[0] == '0|0':
                    if int(i[9].split(':')[1].split(',')[0])<=60:
                        pass
                    elif int(i[9].split(':')[1].split(',')[1])!=0:
                        pass
                if i[9].split(':')[0] == '1/1' or i[9].split(':')[0] == '1|1':
                    if int(i[9].split(':')[1].split(',')[0])!=0:
                        pass
                    elif int(i[9].split(':')[1].split(',')[1])<=60:
                        pass
                elif i[10].split(':')[0] == '0/1' or i[10].split(':')[0] == '0|1' or i[10].split(':')[0] == '1|0':
                    if int(i[10].split(':')[1].split(',')[0])<=60:
                        pass
                    elif int(i[10].split(':')[1].split(',')[1])<=60:
                        pass
                    elif int(i[9].split(':')[3])<=90:
                        pass
                    elif int(i[10].split(':')[3])<=90:
                        pass
                    else:
                        f3.write(l)
                else:pass
    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

