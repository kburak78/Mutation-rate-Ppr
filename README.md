# Mutation-rate-Ppr

## 1. Trimming
Trimming Illumina adapters and Umi adapters from raw reads.
To look into the script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/uni_tools.sh

    for i in A006200178_153621_S1
    A006200178_153622_S2
    A006200178_153623_S3
    A006200178_153624_S4
    A006200178_153625_S5
    A006200178_153626_S6 
    do umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-#in=${i}_L002_R1_001.fastq.gz --stdout=${i}_add_barcode_R1.fastq.gz --read2-stdout #umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-#in=${i}_L002_R3_001.fastq.gz --stdout=${i}_add_barcode_R3.fastq.gz --read2-stdout 
    done
    
## 2. Mapping 
Mapping trimmed reads to the reference genome.
To look into script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_linda.sh 

    for i in A006200178_153621_S1
    A006200178_153622_S2
    A006200178_153623_S3
    A006200178_153624_S4
    A006200178_153625_S5
    A006200178_153626_S6 
    do 
    /NVME/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAA$ bwa mem -t 40 /home/shangao/Data/juicer/Ppr/Ppr_withoutchange/review/Ppr.FINAL.sort.fasta /home/linda/Data/map$ 
    python /home/shangao/script/python/umi_tools_change_reads_title.py -s ./${i}.sam -o ./${i}.replace.sam 
    LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH 
    export LD_LIBRARY_PATH 
    bowtie --threads 40 -v 2 -m 10 -a /home/shangao/Data/juicer/Ppr/Ppr_withoutchange/review/Ppr.FINAL.sort -1 ./$ 
    samtools view -bS ./${i}.sam > ./${i}.bam 
    samtools sort ./${i}.bam -o ./${i}.sort.bam 
    samtools index ./${i}.sort.bam 
    umi_tools dedup -I ./${i}.sort.bam --output-stats=deduplicated --paired -S ./${i}.sort.de.bam --temp-dir=tmp done
  command: sh clean_linda.sh
  
  
 ## 3. SNP calling with GATK4
 Genotype each line seperately with HaplotypeCaller.
 To look into the script: cat/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/gatk.sh 
 
 Software:  
 gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk 
 data 
 ref=/home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa #bam=$1 #out1=$2 #out2=$3 ###build dict for genome #/NVME/Software/popgen/gatk-4.1.9.0/gatk 
 CreateSequenceDictionary -R Ppr.FINAL.fa -O Ppr.FINAL.dict 
 call gvcf 
 $gatk HaplotypeCaller \-R $ref \--emit-ref-confidence GVCF \-I $bam \-O $out1
 
 detect SNPs 
 $gatk GenotypeGVCFs \-R $ref \-V $out1 \-O $out2
 
 Compress 
 bgzip -f $out2 #tabix -p vcf $out2.gz
 
 command:
 
    sh /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/gatk.sh
 
 ## 4. Setting hard filters with parameters as follows: 
QD (variant confidence standardized by depth) <2 
MQ (mapping quality of a SNP) <30 
FS (strand bias in support for REF versus ALT allele calls) >60 
SOR (sequencing bias in which one DNA strand is favored over the other) >5 
ReadPosRankSumTest < -8 
also Indel: 
QD <2 
FS >100 
SOR >5 
ReadPosRankSumTest < -8

To look into the script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/merged_gvcf.sh 

/NVME/Software/popgen/gatk-4.1.9.0/gatk GenotypeGVCFs \-R /home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa \
-V merged_gvcf/merge.g.vcf.gz \
-O merged_gvcf/merge.vcf.gz

/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
-select-type SNP \
-V merged_gvcf/merge.vcf.gz \
-O merged_gvcf/merge.snp.vcf.gz

/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants \
-select-type INDEL \
-V merged_gvcf/merge.vcf.gz \
-O merged_gvcf/merge.indel.vcf.gz

/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantFiltration \
-V merged_gvcf/merge.snp.vcf.gz \
--filter-expression "QD <2.0 || MQ <30.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
--filter-name "PASS" \
-O merged_gvcf/merge.snp.f.vcf.gz

###filter INDEL #/NVME/Software/popgen/gatk-4.1.9.0/gatk VariantFiltration \
-V merged_gvcf/merge.indel.vcf.gz \
--filter-expression "QD <2.0 || FS >100.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
--filter-name "PASS" \
-O merged_gvcf/merge.indel.f.vcf.gz

merge SNP INDEL #/NVME/Software/popgen/gatk-4.1.9.0/gatk MergeVcfs \
-I merged_gvcf/merge.snp.f.vcf.gz \
-I merged_gvcf/merge.indel.f.vcf.gz \
-O merged_gvcf/merge.f.vcf.gz

Command: sh merged-gvcf.sh

## 5. Selecting of the variants 
To look into script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/merged_gvcf/filter.BIALLELIC.sh

/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants -V merge.snp.f.vcf.gz --restrict-alleles-to BIALLELIC -O merge.snp.f.bi.vcf.gz

command:  sh filter.BIALLELIC.sh

## VCFtools for mutation rates: comparing mothers and daughters to get the mutation rate
to look into script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/merged_gvcf/paired/extract_SNP.sh 

vcftools --gzvcf ../merge.snp.f.bi.vcf.gz \--recode-INFO-all \
--maxDP 30 \
--minDP 10 \
--minQ 30 \
--recode \
--stdout \
--maf 0.05 \
--min-meanDP 20 \
--max-missing 0.95 \
--indv A006200178_153625_S5 \
--indv A006200178_153626_S6 \
--out s56.vcf > s56.vcf

grep -v '|' s56.vcf > s56.removedshuxian.vcf

#bedtools intersect -a s56.removedshuxian.vcf -b /home/shangao/Data/EDTA/Ppr/Ppr/Ppr.FINAL.fa.mod.EDTA.TEanno.gff3 -v > s56.removedshuxian.afterbed.vcf

#python /home/shangao/script/python/vcf_filter-same.py -s s56.removedshuxian.afterbed.vcf -o s56.removedshuxian.afterbed.removesame.vcf

#python /home/shangao/script/python/vcf_filter-same.1.py -s s56.removedshuxian.afterbed.removesame.vcf -o s56.removedshuxian.afterbed.removesame.filterread.vcf

#command: sh extract_SNP.sh #result: final VCF files

## 7. Principal Component analysis (pca) to compare the samples to each other
To  look into the script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/merged_gvcf/paired/pca/test.sh

vcftools --gzvcf ../../merge.snp.f.bi.vcf.gz --plink --out merge.snp.f.bi.remove.repeat.vcf #/home/shangao/Software/plink --noweb --file merge.snp.f.bi.remove.repeat.vcf --make-bed --out merge.snp.f.bi.remove.repeat.vcf_bfile #/home/shangao/Software/plink --threads 16 --bfile merge.snp.f.bi.remove.repeat.vcf_bfile --pca 3 --out merge.snp.f.bi.remove.repeat.vcf_pca3_bfile #library(ggplot2) #library(ggrepel) #data<-read.table("merge.snp.f.bi.remove.repeat.vcf_pca3_bfile.eigenvec",header=T) #> pdf('pca12.pdf') #> ggplot(data,aes(x=pc1,y=pc2, colour = col,shape=col))+ geom_point()+ geom_text_repel(aes(x=pc1,y=pc2,label=name))+scale_x_continuous(limits=c(-1, 1)) #> dev.off() #pdf
2

#> pdf('pca13.pdf') #> ggplot(data,aes(x=pc1,y=pc3, colour = col,shape=col))+ geom_point()+ geom_text_repel(aes(x=pc1,y=pc3,label=name))+scale_x_continuous(limits=c(-1, 1)) #> dev.off() #pdf
2

#> pdf('pca23.pdf') #> ggplot(data,aes(x=pc2,y=pc3, colour = col,shape=col))+ geom_point()+ geom_text_repel(aes(x=pc2,y=pc3,label=name))+scale_x_continuous(limits=c(-1, 1)) #> dev.off()

#command: sh test.sh #result: pca image

## 8. Landscape of the mutation to visualize the number of mutations across the genome
To look into the script: cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/merged_gvcf/paired/plot_dot/test1.sh

python /home/shangao/script/python/calulate_SNp_num.py -s ../s56.removedshuxian.afterbed.removesame.filterread.vcf -o s56.slide

#pdf('All.pdf',width=20) #> ggplot(data,aes(x=pos,y=num, colour = chr,shape=as.factor(shape)))+ geom_point(size=3)+ scale_y_continuous(limits=c(0, 20),expand=c(0,0))+scale_x_continuous(expand=c(0,0))+theme(panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5),axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"))+ggtitle("Number of mutation") #dev.off()

#command: sh test1.sh


