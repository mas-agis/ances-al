# ances-al *UNDER DEVELOPMENT*
R based for assigning ancestral alleles and comparison based on WGS data. 
Main inputs are from VCFTools frequency for SNP variants. Example linux command for preparing inputs as follow, for example: 
for i in {1..29};
do vcftools --vcf ~/data/Cattle/filteredSNP.vcf --chr $i --freq --keep taurus_list.txt --out ~/data/Cattle/taurus_Chr_$i;
done
Given vcf file with name filteredSNP.vcf in the folder of ~/data/Cattle, taurus_list.txt is simple text file containing individual ids from group/population of interest in vcf file. Above command will output frequency spectrum of intended group/population per chromosome 1 to 29 in format of "taurus_Chr_$i", where group is the name of group or population of interest, Chr literally written as it is, and i is the number of the chromosome. 
Workflows are starting from AncesAll function to assign the ancestral alleles with default taking major allele with frequency of 1. 
