## 8. Creating a Ideogram with mutations 

Created an chromosome Ideogramm with S12, S34 and S56 mutations. I had to be on shan's account to plot the diagramm, because otherwise the R package wouldn't load or the scriopt wouldn't work. 

Data in /Scratch/gaoshan/breaker/05RIdeogram_plot/potato/
    
    conda activate BRAKER
    R
    > require(RIdeogram)
    Loading required package: RIdeogram
    > mut <- read.table("mutations", sep = "\t", header = T, stringsAsFactors = F)
    > karyo <- read.table("karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
    > ideogram(karyotype = karyo, label = mut, label_type = "marker")
    > convertSVG("chromosome.svg", device = "png")
    I had to errors, which one I solved by switching to Shan's account and activating BRAKER. The other one was solves like this: 

    Error: C stack usage is too close to the limit 
    $ ulimit -s # print default
    8192
    $ R --slave -e 'Cstack_info()["size"]'
       size 
    8388608
    $ ulimit -s 16384 # enlarge stack limit to 16 megs
    $ R --slave -e 'Cstack_info()["size"]'
        size 
    16777216 
    

