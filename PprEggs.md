# Ppr Eggs 

# 1. Samples 

  Eggs from Mother 31
  
# 2. Making the VCF-files 

  ## 2.1 Trimming and Extracting UMI-Adapters 
  
  
  ## 2.2 Mapping
  
  ## 2.2.1 Calculating mean coverage 
  
  ## 2.3 Creating GVCF and VCF-Files 
  
  ## 2.4 Filtering VCF 
  
  Filtered DP, GT, GQ, PL and AD. 
  
  ### Filter DP
  
  On Motoko in /Scratch/linda/Kat_MutRate/Eggs_from_M31_Ppr/vcf/Comparison_of_Daughters/filter.eggs.DP.py
  
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

    ##  pass every site in which the DP of the samples is '.' or not between 50% and 200% of the average DP 
    for l in f1.readlines():
        if '#' in l:
            pass
        else:
            i=l.strip().split()
            if i[10].split(':')[0][0] == '.':
                pass
            elif i[9].split(':')[0][0] == '.':
                pass
            elif 29 <= int(i[9].split(':')[2]) <= 116 and 29 <= int(i[10].split(':')[2]) <= 116 and 29 <= int(i[11].split(':')[2]) <=116:
                f3.write(l)
            else: pass

    f1.close()
    #f2.close()
    f3.close()
    #f4.close()


  ### Filter GT 
  
  On Motoko in /Scratch/linda/Kat_MutRate/Eggs_from_M31_Ppr/vcf/Comparison_of_Daughters/filter.eggs.hap.py
  

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
            # pass every site that where all three samples have the same Genotype
            if i[9].split(':')[0] == i[10].split(':')[0] and i[9].split(':')[0] == i[11].split(':')[0]:
                pass
            # pass every site where all samples are homozygous for the alternative allele
            elif i[9].split(':')[0] == '1|1' or i[9].split(':')[0] == '1/1':
                if i[10].split(':')[0] == '1|1' or i[10].split(':')[0] == '1/1':
                    if i[11].split(':')[0] == '1|1' or i[11].split(':')[0] == '1/1':
                        pass
                    else: f3.write(l)
                else: f3.write(l)
            # pass every site where all samples are heterozygous for alternative and reference
            elif i[9].split(':')[0] == '0|1' or i[9].split(':')[0] == '0/1':
                if i[10].split(':')[0] == '0|1' or i[10].split(':')[0] == '0/1':
                    if i[11].split(':')[0] == '0|1' or i[11].split(':')[0] == '0/1':
                        pass
                    else: f3.write(l)
                else: f3.write(l)
            #pass every site where all samples are heterozygous  for the second alternative allele and reference
            elif i[9].split(':')[0] == '0|2' or i[9].split(':')[0] == '0/2':
                if i[10].split(':')[0] == '0|2' or i[10].split(':')[0] == '0/2':
                    if i[11].split(':')[0] == '0|2' or i[11].split(':')[0] == '0/2':
                        pass
                    else: f3.write(l)
                else: f3.write(l)
            #everything else will be written
            else: f3.write(l)
    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

    
  ### Filter GQ 
  
  On Motoko in /Scratch/linda/Kat_MutRate/Eggs_from_M31_Ppr/vcf/Comparison_of_Daughters/filter.eggs.GQ.py
  
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

            #if GQ is '.' pass
            if str(i[9].split(':')[3]) == '.' or str(i[10].split(':')[3]) == '.' or str(i[11].split(':')[3]) == '.' :
                pass
            #if GQ is not 99 or 0 pass. Why did I left 0? Because sometimes the most likely PL and second likely PL are 0, which menas both GT have the same high likelihood. If that happens the VCF-file will show a GQ of 0.
            elif int(i[9].split(':')[3]) != 99 and int(i[9].split(':')[3]) != 0:
                pass
            elif int(i[10].split(':')[3]) != 99 and int(i[10].split(':')[3]) != 0:
                pass
            elif int(i[11].split(':')[3]) != 99 and int(i[11].split(':')[3]) != 0:
                pass
            #everything else will be written 
            else: f3.write(l)


    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

    
  ### Filter PL
  
  On Motoko in /Scratch/linda/Kat_MutRate/Eggs_from_M31_Ppr/vcf/Comparison_of_Daughters/filter.eggs.GQ.PL.py
  
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

            #if PL is in the 7th position in FORMAT
            if 'PGT' in i[8]:
                # listing contents of PL by increasing order for all three samples
                # for example if PL is '100, 500, 0'; it will be listed as '0, 100, 500' 
                l0 = i[9].split(':')[6].split(',')
                if '.' in l0: pass
                else: PL0 = sorted([int(x) for x in l0])

                l1 = i[10].split(':')[6].split(',')
                if '.' in l1: pass
                else: PL1 = sorted([int(x) for x in l1])

                l2 = i[11].split(':')[6].split(',')
                if '.' in l2: pass
                else: PL2 = sorted([int(x) for x in l2])

                # if the second PL is not 0, calculate the difference between first and second likely PL. A PL difference of a 100 or less will be passed.
                if int(PL0[1]) != 0 and (int(PL0[1]) - int(PL0[0])) =< 100:
                    pass
                # if the second PL is 0, calculate the difference between second and third likely PL. A PL difference of a 100 or less will be passed.
                elif int(PL0[1]) == 0 and int(PL0[2]) - int(PL0[1]) =< 100:
                    pass
                #same as above for second sample
                elif int(PL1[1]) != 0 and (int(PL1[1]) - int(PL1[0])) =< 100:
                    pass
                elif int(PL1[1]) == 0 and int(PL1[2]) - int(PL1[1]) =< 100:
                    pass
                #same as above for thrid sample
                elif int(PL2[1]) != 0 and (int(PL2[1]) - int(PL2[0])) =< 100:
                    pass
                elif int(PL2[1]) == 0 and int(PL2[2]) - int(PL2[1]) =< 100:
                    pass
                #everything else will be writtem
                else: f3.write(l)

            # same as above, for if PL is in 5th position in FORMAT
            else: 
                l3 = i[9].split(':')[4].split(',')
                if '.' in l3: pass
                else: PL3 = sorted([int(x) for x in l3])

                l4 = i[10].split(':')[4].split(',')
                if '.' in l4: pass
                else: PL4 = sorted([int(x) for x in l4])

                l5 = i[11].split(':')[4].split(',')
                if '.' in l5: pass
                else: PL5 = sorted([int(x) for x in l5])

                if int(PL3[1]) != 0 and (int(PL3[1]) - int(PL3[0])) =< 100:
                    pass
                elif int(PL3[1]) == 0 and int(PL3[2]) - int(PL3[1]) =< 100:
                    pass
                elif int(PL4[1]) != 0 and (int(PL4[1]) - int(PL4[0])) =< 100:    
                    pass
                elif int(PL4[1]) == 0 and int(PL4[2]) - int(PL4[1]) =< 100:
                    pass
                elif int(PL5[1]) != 0 and (int(PL5[1]) - int(PL5[0])) =< 100:
                    pass
                elif int(PL5[1]) == 0 and int(PL5[2]) - int(PL5[1]) =< 100:
                    pass
                else: f3.write(l)



    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

  ### Filter AD 

  On Motoko in /Scratch/linda/Kat_MutRate/Eggs_from_M31_Ppr/vcf/Comparison_of_Daughters/filter.eggs.GQ.PL.py

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
        a = l.strip().split()[9].split(':')[1].split(',')
        b = l.strip().split()[10].split(':')[1].split(',')
        c = l.strip().split()[11].split(':')[1].split(',')
        #number of reads supporting non-GT in homozygous
        x = 6
        #number of reads supporting GT in homozygous
        y = 48

        ####################### First Sample #################################

        #if the Genotype is homozygous
        if '0|0' or  '1|1' or '2|2' or '0/0' or '1/1' or '2/2' in l.strip().split()[9].split(':')[0]:
            #if there are three AD values, filter so one is supported with at most with y reads & the other two with at least y reads
            if len(a) == 3:
                if int(a[0]) > x and int(a[1]) > x and int(a[2]) < y:
                        pass
                elif int(a[0]) > x and int(a[1]) < y and int(a[2]) > x:
                        pass
                elif int(a[0]) < y and int(a[1]) > x and int(a[2]) > x:
                        pass
                #if there are two AD values, pass according to parameters x and y
            elif len(a) == 2:
                    if int(a[0]) > x and  int(a[1]) < y:
                        pass
                    elif int(a[1]) > x  and  int(a[0]) < y:
                        pass
                    else: f3.write(l)
        #if the GT is heterozygous
        else: 
            # if there are 2 AD values, both should be at least above 28 reads 
            if len(a) == 2:
                if int(a[0]) < 27 and int(a[1]) < 27:
                    pass
                else: f3.write(l)
            # if there are 3 AD values, two should be at least above 28 and one should be not equal and higher then x
            elif len(a) == 3:
                if int(a[0]) < 27 and int(a[1]) < 27 and int(a[2]) > x:
                    pass
                elif int(a[0]) < 27 and int(a[1]) > x and int(a[2]) < 27:
                    pass
                elif int(a[0]) > x and int(a[1]) < 27 and int(a[2]) < 27:
                    pass
                else: f3.write(l)
            else: pass

       ######################## Second Sample #################################
       #same as above just for the second sample
        if '0|0' or  '1|1' or '2|2' or '0/0' or '1/1' or '2/2' in l.strip().split()[10].split(':')[0]:
            if len(b) == 3:
                if int(b[0]) > x and int(b[1]) > x and int(b[2]) < y:
                        pass
                elif int(b[0]) > x and int(b[1]) < y and int(b[2]) > x:
                        pass
                elif int(b[0]) < y and int(b[1]) > x and int(b[2]) > x:
                        pass
            elif len(b) == 2:
                    if int(b[0]) > x and  int(b[1]) < y:
                        pass
                    elif int(b[1]) > x  and  int(b[0]) < y:
                        pass
                    else: f3.write(l)
        else:
            if len(b) == 2:
                if int(b[0]) < 27 and int(b[1]) < 27:
                    pass
                else: f3.write(l)
            elif len(b) == 3:
                if int(b[0]) < 27 and int(b[1]) < 27 and int(b[2]) > x:
                    pass
                elif int(b[0]) < 27 and int(b[1]) > x and int(b[2]) < 27:
                    pass
                elif int(b[0]) > x and int(b[1]) < 27 and int(b[2]) < 27:
                    pass
                else: f3.write(l)
            else: pass

        ########################### Third sample #####################################

        if '0|0' or  '1|1' or '2|2' or '0/0' or '1/1' or '2/2' in l.strip().split()[11].split(':')[0]:
            if len(c) == 3:
                if int(c[0]) > x and int(c[1]) > x and int(c[2]) < y:
                        pass
                elif int(c[0]) > x and int(c[1]) < y and int(c[2]) > x:
                        pass
                elif int(c[0]) < y and int(c[1]) > x and int(c[2]) > x:
                        pass
            elif len(c) == 2:
                    if int(c[0]) > x and  int(c[1]) < y:
                        pass
                    elif int(c[1]) > x  and  int(c[0]) < y:
                        pass
                    else: f3.write(l)
        else:
            if len(c) == 2:
                if int(c[0]) < 27 and int(c[1]) < 27:
                    pass
                else: f3.write(l)
            elif len(c) == 3:
                if int(c[0]) < 27 and int(c[1]) < 27 and int(c[2]) > x:
                    pass
                elif int(c[0]) < 27 and int(c[1]) > x and int(c[2]) < 27:
                    pass
                elif int(c[0]) > x and int(c[1]) < 27 and int(c[2]) < 27:
                    pass
                else: f3.write(l)
            else: pass


    f1.close()
    #f2.close()
    f3.close()
    #f4.close()

  ## 2.5 Result 
  
# 3. Loss of heterogozity 

  ## 3.1 Creating VCF-File of Mother 31 with each egg 
  
  ## 3.2 Filtering VCFs 
  
  Hardfiltering
  
  filtering DP 10=<x =< 32 (50% and 150% of the average coverage 21,6154)
  
  filtering from Mom (hetero) to Daughter (homo) -> loss of heterogozity
  
  filtering GQ
