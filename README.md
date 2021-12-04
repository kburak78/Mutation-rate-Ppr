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
    
## 2. Cat previous data and new data

Location: /home/linda/Data/all_data/mapping.sh 

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
  
  
 
