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
4. [Detect de novo mutations](#Mutations)
    1. [Filter DP](#DP)
    2. [Filter GT](#GT)
    3. [Filter GQ & PL](#GQ)
    4. [Filter AD](#AD)
5. [Callable genome estimation](#CG)
    1. [Filter DP](#DP1)
    2. [Filter GT](#GT1)
    3. [Filter GQ & PL](#GQ1)
    4. [Filter AD](#AD1)
    5. [Count positions](#count)
6. [Calculating the mutation rate](#Mutrate)

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
## 4. Detect de novo mutations <a name="Mutations"> </a>

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
            #Mother must be homozygous, Daughter must be heterozygous. Everything else will be passed.
            if i[9].split(':')[0] == '0/0' or i[9].split(':')[0] == '0|0' or i[9].split(':')[0] == '1/1' or i[9].split(':')[0] == '1|1':
                if i[10].split(':')[0] == '0/1' or i[10].split(':')[0] == '0|1' or i[10].split(':')[0] == '1|0' or i[10].split(':')[0] == '2|1' or i[10].split(':')[0] == '1|2' or i[10].split(':')[0] == '1/2' or i[10].split(':')[0] == '0|2' or i[10].split(':')[0] == '2|0':
                    f3.write(l)
                else:pass 
            else:pass

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

### iii. Filter GQ and PL <a name="GQ"> </a>

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

            #GQ equal or under 98 will be passed
            if int(i[9].split(':')[3]) <= 98 or int(i[10].split(':')[3]) <= 98:
                pass

            #if PL is in the 7th position in FORMAT
            if 'PGT' in i[8]:
                #PL values are sorted in increasing order for both samples
                #e.g. a PL like this [100, 0, 500], will be sorted to [0, 100, 500]
                l2 = i[9].split(':')[6].split(',')
                if '.' in l2: pass
                else: PL2 = sorted([int(x) for x in l2]) 

                l3 = i[10].split(':')[6].split(',')
                if '.' in l3: pass
                else: PL3 = sorted([int(x) for x in l3])      

                #if second likeliest PL is 0, the difference between the third and second likeliest PL is calculated. If the second likeliest PL is not zero the difference between second likeliest and most likely PL is calculated. For the Mother the must be a difference over 200 and for the daughter a difference of over 100. 
                if int(PL2[1]) == 0 and  int(PL2[2]) - int(PL2[1]) < 120:
                        pass
                elif int(PL3[1]) == 0 and int(PL3[2]) - int(PL3[1]) < 100:
                        pass
                elif int(PL2[1]) != 0 and int(PL2[1]) - int(PL2[0]) < 120:
                    pass
                elif int(PL3[1]) != 0 and int(PL3[1]) - int(PL3[0]) < 100:
                    pass
                else: f3.write(l)
            else: 
                l0 = i[9].split(':')[4].split(',')
                if '.' in l0: pass
                else: PL0 = sorted([int(x) for x in l0])

                l1 = i[10].split(':')[4].split(',')
                if '.' in l1: pass
                else: PL1 = sorted([int(x) for x in l1])

                if int(PL0[1]) == 0 and int(PL0[2]) - int(PL0[1]) < 120:
                        pass
                elif int(PL1[1]) == 0 and int(PL1[2]) - int(PL1[1]) < 100:
                        pass
                elif int(PL0[1]) != 0 and int(PL0[1]) - int(PL0[0]) < 120:
                    pass
                elif int(PL1[1]) != 0 and int(PL1[1]) - int(PL1[0]) < 100: 
                    pass
                else: f3.write(l)

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()



### iv. Filter AD <a name="AD"> </a>

cat filter.AD.Mom.py

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
            i=l.strip().split()[9].split(':')[1].split(',')
            # if there are 3 AD values for the Mother, one should be over 59   
            if len(i) == 3:
                if i[0] == '0' and i[1] == '0' :
                    f3.write(l)
                elif i[0] == '0'  and i[2] == '0':
                    f3.write(l)
                elif  i[1] == '0' and i[2] == '0':
                    f3.write(l)
                else: pass
            elif len(i) == 2:
                if i[0] == '0':
                    f3.write(l)
                elif i[1] == '0' :
                    f3.write(l)
                else: pass
            else: print(l)

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()


cat filter.AD.Daughter1.py

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
            # Sum of all AD reads must equal the number of total reads (DP)
            if len(i[10].split(':')[1].split(',')) == 2:
                if int(i[10].split(':')[1].split(',')[0]) + int(i[10].split(':')[1].split(',')[1]) == int(i[10].split(':')[2]):
                    f3.write(l)
                else: pass
            elif len(i[10].split(':')[1].split(',')) == 3:
                if int(i[10].split(':')[1].split(',')[0]) + int(i[10].split(':')[1].split(',')[1]) + int(i[10].split(':')[1].split(',')[2]) == int(i[10].split(':')[2]):
                    f3.write(l)
                else: pass
            else: print(l)

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()
    
  
  cat filter.AD.Daughter2.py 
  
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
            i=l.strip().split()[10].split(':')[1].split(',')
            if len(i) == 2:
                if int(i[0]) > 27 and int(i[1]) > 27:
                    f3.write(l)
                else: pass
            elif len(i) == 3:
                if int(i[0]) > 27 and int(i[1]) > 27:
                    f3.write(l)
                elif int(i[0]) > 27 and int(i[2]) > 27:
                    f3.write(l)
                elif int(i[2]) > 27 and int(i[1]) > 27:
                    f3.write(l)
                else: pass
            else: print(l)

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

- - - -
## 5. Callable Genome estimation <a name="CG"> </a>

adding header to hf and gcov filtered data:

cat /RAID/Data/linda/Mother_Egg_Pairs_Ppr/vcf/merged_gvcf/filter/SNPs_after_gcov_filter/addheader.sh

### i. Filter DP <a name="DP1"> </a>

cat filter.hap.DP.py

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
        if '#' in l:
            pass
        else:
            i=l.strip().split()
            if i[9].split(':')[0] == '0/0' or i[9].split(':')[0] == '0|0' or i[9].split(':')[0] == '1/1' or i[9].split(':')[0] == '1|1':
                if 'AD' not in i[8]:
                    if 31<=int(i[9].split(':')[1])<=124:
                        f3.write(l)
                    else:pass
                elif 31<=int(i[9].split(':')[2])<=124:
                        f3.write(l)
                else:pass
            else:pass
    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

### ii. Filter GT <a name="GT1"> </a>

see i. filter DP 

### iii. Filter GQ and PL <a name="GQ1"> </a>

cat filter.GQ.py
 
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
            if 'AD' not in i[8]:
                if int(i[9].split(':')[2]) <= 98:
                    pass
                else: f3.write(l)
            elif 'AD' in i[8]:
                if int(i[9].split(':')[3]) <= 98:
                    pass
                else: f3.write(l)

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()
    
cat filter.PL.py
    
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
                if 'PGT' in i[8]:
                    l0 = i[9].split(':')[6].split(',')
                    if '.' in l0: pass
                    else: PL0 = sorted([int(x) for x in l0])

                    if int(PL0[1]) == 0 and  int(PL0[0]) == 0:
                        if int(PL0[2]) - int(PL0[1]) < 120:
                            pass
                        else: f3.write(l)
                    elif int(PL0[1]) - int(PL0[0]) < 120:
                        pass
                    else: f3.write(l)

                else:
                    l1 = i[9].split(':')[4].split(',')
                    if '.' in l1: pass
                    else: PL1 = sorted([int(x) for x in l1])

                    if int(PL1[1]) == 0 and  int(PL1[0]) == 0:
                        if int(PL1[2]) - int(PL1[1]) < 120:
                            pass
                        else: f3.write(l)
                    elif int(PL1[1]) - int(PL1[0]) < 120:
                        pass
                    else: f3.write(l)

        f1.close()
        #f2.close()
        f3.close()
        #f4.close()


### iv. Filter AD <a name="AD1"> </a>

 cat filter.AD.py
 
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
                i=l.strip().split()[9].split(':')[1].split(',')

                if 'AD' not in l.strip().split()[8]:
                    f3.write(l)
                else:
                    if 'AD'  in l.strip().split()[8]:
                    # if content is 3 possible outcomes are
                        if len(i) == 3:
                            if i[0] == '0' and i[1] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[2] == '0':
                                f3.write(l)
                            elif i[1] == '0' and i[2] == '0':
                                f3.write(l)
                            else: pass
                        elif len(i) == 2:
                            if i[0] == '0':
                                f3.write(l)
                            elif i[1] == '0' :
                                f3.write(l)
                            else: pass
                        elif len(i) == 4:
                            if i[0] == '0' and i[1] == '0' and i[3] == '0':
                                f3.write(l)
                            elif i[0] == '0'  and i[2] == '0' and i[3] == '0':
                                f3.write(l)
                            elif  i[1] == '0' and i[2] == '0' and i[3] == '0':
                                f3.write(l)
                            elif  i[1] == '0' and i[2] == '0' and i[0] == '0':
                                f3.write(l)
                            else: pass
                        elif len(i) == 5:
                            if i[0] == '0' and i[1] == '0' and i[3] == '0' and i[4] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[2] == '0' and i[3] == '0' and i[4] == '0':
                                f3.write(l)
                            elif i[1] == '0' and i[2] == '0' and i[3] == '0' and i[4] == '0':
                                f3.write(l)
                            elif i[1] == '0' and i[2] == '0' and i[0] == '0'and i[3] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[1] == '0' and i[2] == '0'and i[4] == '0':
                                f3.write(l)
                            else: pass
                        elif len(i) == 6:
                            if i[0] == '0' and i[1] == '0' and i[2] == '0' and i[3] == '0' and i[4] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[2] == '0' and i[3] == '0' and i[4] == '0' and i[5] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[1] == '0' and i[3] == '0' and i[4] == '0' and i[5] == '0':
                                f3.write(l)
                            elif i[1] == '0' and i[2] == '0' and i[3] == '0' and i[4] == '0' and i[5] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[1] == '0' and i[2] == '0'and i[3] == '0'and i[5] == '0':
                                f3.write(l)
                            elif i[0] == '0' and i[1] == '0' and i[2] == '0'and i[4] == '0' and i[5] == '0':
                                f3.write(l)
                            else: pass

        f1.close()
        #f2.close()
        f3.close()
        #f4.close()


### v. Count positions <a name="count"> </a>

cat count_positions.sh

    for i in s1 \
        s3 \
        s5
    do

    python count_positions.py -s ../filter.DP.GT/${i}.DP.GT.g.vcf -o ${i}.DP.GT.txt
    awk '{sum+= $1} END {print sum}' ${i}.DP.GT.txt > ${i}.DP.GT.positions

    python count_positions.py -s ../filter.DP.GT.GQ/${i}.DP.GT.GQ.g.vcf -o ${i}.DP.GT.GQ.txt
    awk '{sum+= $1} END {print sum}' ${i}.DP.GT.GQ.txt > ${i}.DP.GT.GQ.positions

    python count_positions.py -s ../filter.DP.GT.GQ.PL/${i}.DP.GT.GQ.PL.g.vcf -o ${i}.DP.GT.GQ.PL.txt
    awk '{sum+= $1} END {print sum}' ${i}.DP.GT.GQ.PL.txt > ${i}.DP.GT.GQ.PL.positions

    python count_positions.py -s ../filter.DP.GT.GQ.PL.AD/${i}.DP.GT.GQ.PL.AD.g.vcf -o ${i}.DP.GT.GQ.PL.AD.txt
    awk '{sum+= $1} END {print sum}' ${i}.DP.GT.GQ.PL.AD.txt > ${i}.DP.GT.GQ.PL.AD.positions

    cat ${i}.DP.GT.positions ${i}.DP.GT.GQ.positions ${i}.DP.GT.GQ.PL.positions ${i}.DP.GT.GQ.PL.AD.positions > ${i}.all.positions

    rm ${i}.DP.GT.txt
    rm ${i}.DP.GT.GQ.txt
    rm ${i}.DP.GT.GQ.PL.txt
    rm ${i}.DP.GT.GQ.PL.AD.txt

cat count_positions.py

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
            if '.' in i[7]:
                f3.write('1\n')
            elif int(i[7].split('=')[1]) == int(i[1]):
                f3.write('1\n')
            else: f3.write(str(int(i[7].split('=')[1])-int(i[1]))+'\n')

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()



- - - -
## 6. Calculating the mutation rate <a name="Mutrate"> </a> 

The Mutation rate can be calculated as: 

Mutation counts / (2 x callable sites) = Mutation rate 

With the first calculations something between 1-2x10^-6 comes out which would be higher then for example spider mites. If Ppr has a lower mutation rate we would expect something under around 1*10^-9. 

### Identifying Coding and Non-coding regions with existing gtf-file

extract coding regions from gtf-file.

    awk -v OFS='\t' '$3=="gene"{print$1,$4,$5}' /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf > coding_area_Ppr
    
    coding_area.sh

/RAID/Data/linda/all_data/mapped_data/genome_coverage/coding_area_gcov 
