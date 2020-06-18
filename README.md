# ances-al 
R based for assigning ancestral alleles and comparison of alleles frequency group/population of interest to defined ancestral alleles based on WGS data. \n
##INPUTS 
Main inputs are from VCFTools frequency for SNP variants. Example linux command for preparing inputs as follow, for example: \n
for i in {1..29}; \n
do vcftools --vcf ~/data/Cattle/filteredSNP.vcf --chr $i --freq --keep taurus_list.txt --out ~/data/Cattle/taurus_Chr_$i; \n
done\n
Given vcf file with name filteredSNP.vcf in the folder of ~/data/Cattle, taurus_list.txt is simple text file containing individual ids from group/population of interest in vcf file. Above command will output frequency spectrum of intended group/population per chromosome 1 to 29 in format of "taurus_Chr_$i", where group is the name of group or population of interest, Chr literally written as it is, and i is the number of the chromosome. \n
\n
##Function explanations
Workflows are starting from AncesAll function to assign the ancestral alleles with default taking major allele with frequency of 1. Then R-functions work accordingly to the explanation below.\n \n
#Function ancealls - calling putative ancestral alleles from freq files. Inputs for running this function are: \n
gr: group/population name of respective freq files \n
n: number of individuals in the group/population respective to freq files \n
cr: chromosome number to be called respective to freq files \n
fr: allele frequency to be called as ancestral (default is 1) \n
dip: Ploidy (default is 2) \n
This function will output a list containing two dataframes, i.e. Ances_Allele and summary. Ances_Allele comprises of  5 columns and each row is for physical location of SNP, i.e.  Chr:Chromosome number, Pos: Position, Alleles_n: number of alleles, AA: Defined ancestral allele, Freq: Frequency of alleles. Summary dataframe containing columns for each A, C, G, T allele and a single row for sum of defined ancestral allele for respective bases. \n \n 
#Function persist.ancesA - counting sum of ancestral alleles in  a group/population of interest within fixed size of scanning window. Inputs for running this function are: \n
AncesAA: First dataframe from ancealls function \n
gr: name of group/population respective to freq files \n
n: number of individuals in the group/population respective to freq files \n
cr: chromosome number to be called respective to freq files \n
dip: Ploidy (default is 2) \n
scan: size of canning window (default:10000) \n
o: First top treshold percentage (default: 0.1%) \n
p: Second top treshold percentage (default: 0.01%) \n
q: Third top treshold percentage (default: 0.001%) \n
This function will output a list with two dataframes, i.e. AAcounts and treshold. AAcounts is the details of conserved ancestral alleles in population of interest. This dataframe contains rows according to number windows with ancestral alleles and 6 columns, i.e. Chr:Chromosome number, Start:Start position of scanning window, End:End position of scanning window, Window:Consecutive number of windows with ancestral alleles, Score: 1 if current allel match with ancestral allel otherwise 0, Ancestral_count: sum of sites with ancestral alleles within current window. Treshold is numeric treshold for ancestral allele counts with columns bar1, bar2, bar3, for respective o, p, q parameter for each chromosome. \n \n  
#Function of annv.ances - separate windows from persist.ancesA output based on ancestral alleles count. Inputs for running this function are: \n
AAcount: First dataframe from ancealls function \n
treshold: Second dataframe from ancealls function \n
bar: treshold for extracting windows with high ancestral allele (default: bar1) \n
This function will output three dataframes, i.e. Above, Zero, and Ratio. Above and zero dataframes are complied for input of annovar for annotation of regions of interest in regards its high ancestral count above the treshold and regions without ancestral allele. Ratio dataframe is ratio of windows without ancestral alleles to total scanning windows per chromosome wise. \n \n 
#Function of plt.ances - creating a simple manhattan plot for windows where ancestral allele exists in the group/population of interest. Inputs for running fuction are: \n
AAcount: First dataframe from ancealls function \n
treshold: Second dataframe from ancealls function \n
cr: chromosome number to be plotted \n
This function will create a manhattan plot based on umber of ancestral alleles found wtin scanning windows and adding horizontal bars corresponding to treshold of bar1, bar2, and bar3. \n \n 
#Function AverageAA - as it name, getting mean of ancestral allele counts for scanning windows. Inputs are: \n
AAcount: First dataframe from ancealls function \n
treshold: Second dataframe from ancealls function \n
Note: windows without ancestral allele are inclusive in mean calculation \n \n
#Function CheckNullRegion - checking whether scanning windows without ancestral allele are due to change of ancestral allele or there were no ancestral alleles defined within the scanning windows. Inputs for running fuction are: \n
AncesAA: First dataframe from ancealls function \n
Region_wo_AA: Dataframe zero from the output of annv.ances function \n
cr: chromosome number \n
This function will output a dataframe with 5 columns, i.e. Chr:Chromosome number, Start:Start position of scanning window, End:End position of scanning window, AA_count: Must be 0 (No ancestral allele detected), Actual_AA_Sites: Actual number of sites with defined ancestral allele. \n 
If 0 is found in column 5 meaning within respective scanning window, there were no ancestral allele defined. Else, if any positive number is found, then mutation is likely happening changing the allele from the defined ancestral alleles. \n

