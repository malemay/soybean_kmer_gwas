#!/usr/bin/gawk -f

# This script takes a VCF in which all alleles are explicitly
#  coded (the sequence of each allele is written) as output
#  by BayesTyper and outputs a new vcf which contains only
#  alleles considered as SVs (here 50 nucleotides or more
#  but should be customizable by an argument eventually)
#  This version of the scripts also keeps any variants
#  for which the alternative allele is >= 50 nt.

# The output is printed to console and should be handled accordingly

# Lines beginning with a # (header) are output no matter what
/^#/ { print $0 }

!/^#/ {
  # Creating an array of alternative alleles
  n_alleles = split($5, alleles, ",")
  
  # This variable will be used as a flag to indicate whether this line should be output
  # We want to output any line for which at least one allele is >= 50 or <= -50
  output = 0

  # Then we iterate over all alleles
  for(i in alleles) {
    allele_length = length(alleles[i]) - length($4)
    if (allele_length >= 50 || allele_length <= -50 || length(alleles[i]) >= 50) output = 1
  }

  # We print the full line if the condition is met
  if(output == 1) print $0
}

