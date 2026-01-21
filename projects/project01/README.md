# BINF6250
# Project01
# Introduction
VCF file parsing

# Pseudocode
input: vcf file
output: dictionary containing frequency of disease types

Open vcf file and read line by line 
for each line determine rarity of variant using allele frequencies from the reported AF_EXAC value in the file
If variant is not rare, do not record variant
If the variant is considered rare, record associated disease
If no disease is associated, do not record
start a running count that tracks the frequency of diseases associated with rare variants
Each time a disease associated with a rare variant is recorded, increase associated disease's tally by 1
When entire file is read, report frequency each disease associate with a rare variant appeared in the file
# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus

