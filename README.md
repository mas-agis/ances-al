# ances-al
R-based for assigning ancestral alleles and breed comparison based on allele frequency of WGS data. \n
Inputs are frequency summary from vcftools function --freq with specific --chr to be assigned. \n
For preparing inputs, one need calling SNP variants and vcf ready-file. In linux environemnt, these commands should be executed using vcftools: \n
  for i in {1..(number of chromosome)};
  do vcftools --vcf (path-to-vcf) --chr $i --freq --keep (list-of-individuals-in-a-population) --out (path-to-output/breeds__Chr_$i);
  done
Workflows are starting from AncesAll function to assign the ancestral alleles with default \n
taking major allele with frequency of 1. 

