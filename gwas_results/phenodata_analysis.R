# This script has a detailed look at the unexpected haplotype/phenotype
# associations based on k-mers to try and find explanations

# Reading the table of sample names and their SRR accesion numbers
srr_ids <- read.table("utilities/srr_id_correspondence.txt")
srr <- srr_ids[[2]]
names(srr) <- srr_ids[[1]]

# Flower color
flower_color <- readRDS("gwas_results/kmers/flower_color_W1_phenodata.rds")
table(flower_color$phenotype, flower_color$haplotype_id)
flower_color[flower_color$phenotype == "White" & flower_color$haplotype_id == 1, -2]
#         sample haplotype_id phenotype
# NA        <NA>           NA      <NA>
# 83  SRR1533389            1     White
# 275    USB-480            1     White
srr["SRR1533389"]
#   SRR1533389 
# "SRR1533389" 
srr["USB-480"]
#                   USB-480 
# "SRR12366485;SRR12366488" 

# the fastqc data for SRR1533389 has a weird GC content, which suggests contamination 
# the two SRA accesions for USB-480 did not suffer from that problem


# Pubescence color (T locus)
pubescence_color_all <- readRDS("gwas_results/kmers/pubescence_color_all_T_phenodata.rds")
table(pubescence_color_all$phenotype, pubescence_color_all$haplotype_id)
pubescence_color_all[pubescence_color_all$phenotype == "Tawny" & pubescence_color_all$haplotype_id == 2, -2]
#    sample haplotype_id phenotype
# 38  HN057            2     Tawny
srr["HN057"]
#                                 HN057 
# "SRR12190489;SRR12190490;SRR12190491" 
# There does not seem to be any problem with the SRR IDs associated with this sample

# Pubescence color (Td locus)
pubescence_color_nogray <- readRDS("gwas_results/kmers/pubescence_color_nogray_Td_phenodata.rds")
table(pubescence_color_nogray$phenotype, pubescence_color_nogray$haplotype_id)
pubescence_color_nogray[pubescence_color_nogray$phenotype == "Light tawny" & pubescence_color_nogray$haplotype_id == 1, -2]
#      sample haplotype_id   phenotype
# NA     <NA>           NA        <NA>
# 6     HN018            1 Light tawny
# 23  USB-029            1 Light tawny
# 51  USB-146            1 Light tawny
# 67  USB-323            1 Light tawny
# 102 USB-530            1 Light tawny
# 134 USB-763            1 Light tawny
# 137 USB-780            1 Light tawny
srr["HN018"]
#        HN018 
# "SRR2163310" 
srr["USB-029"]
#                                                       USB-029 
# "SRR12191316;SRR12191317;SRR12191318;SRR12191319;SRR12191321" 
srr["USB-146"]
#       USB-146 
# "SRR12335713" 
srr["USB-323"]
#                   USB-323 
# "SRR12335744;SRR12335745" 
srr["USB-530"]
#       USB-530 
# "SRR12366168" 
srr["USB-763"]
#       USB-763 
# "SRR12366514" 
srr["USB-780"]
#       USB-780 
# "SRR12366208" 

# there does not seem to be problems with any of those samples

# Seed coat color (green-yellow)
seed_coat_color_greenyellow <- readRDS("gwas_results/kmers/seed_coat_color_greenyellow_G_phenodata.rds")
table(seed_coat_color_greenyellow$phenotype, seed_coat_color_greenyellow$haplotype_id)
seed_coat_color_greenyellow[seed_coat_color_greenyellow$phenotype == "Green" & seed_coat_color_greenyellow$haplotype_id == 1, -2]
#     sample haplotype_id phenotype
# NA    <NA>           NA      <NA>
# 40 USB-419            1     Green
# 95 USB-541            1     Green
srr["USB-419"]
#                   USB-419 
# "SRR12366370;SRR12366371" 
srr["USB-541"]
#       USB-541 
# "SRR12366341" 
seed_coat_color_greenyellow[seed_coat_color_greenyellow$phenotype == "Yellow" & seed_coat_color_greenyellow$haplotype_id == 2, -2]
#      sample haplotype_id phenotype
# 136 USB-762            2    Yellow
srr["USB-762"]
#       USB-762 
# "SRR12366515" 

# Stem termination
stem_termination_all_Dt1 <- readRDS("gwas_results/kmers/stem_termination_all_Dt1_phenodata.rds")
table(stem_termination_all_Dt1$phenotype, stem_termination_all_Dt1$haplotype_id)
stem_termination_all_Dt1[stem_termination_all_Dt1$phenotype %in% c("Indeterminate", "Semi-determinate") & stem_termination_all_Dt1$haplotype_id == 2, -2]
#       sample haplotype_id        phenotype
# NA      <NA>           NA             <NA>
# NA.1    <NA>           NA             <NA>
# 151  USB-156            2    Indeterminate
# 154  USB-162            2 Semi-determinate
# 170  USB-228            2 Semi-determinate
# 200  USB-318            2 Semi-determinate
# 287  USB-504            2 Semi-determinate
# 341  USB-697            2 Semi-determinate
# NA.2    <NA>           NA             <NA>
srr["USB-156"]
#       USB-156 
# "SRR12336102" 
srr["USB-162"]
#                               USB-162 
# "SRR12336091;SRR12336092;SRR12336093" 
srr["USB-228"]
#                   USB-228 
# "SRR12336152;SRR12336153" 
srr["USB-318"]
#                   USB-318 
# "SRR12335672;SRR12335673" 
srr["USB-504"]
#                   USB-504 
# "SRR12366298;SRR12366299" 
srr["USB-697"]
#       USB-697 
# "SRR12366399" 

# Hilum color
hilum_color_all_T <- readRDS("gwas_results/kmers/hilum_color_all_T_phenodata.rds")
table(hilum_color_all_T$phenotype, hilum_color_all_T$haplotype_id)
hilum_color_all_T[hilum_color_all_T$phenotype == "Yellow" & hilum_color_all_T$haplotype_id == 1, -2]
#      sample haplotype_id phenotype
# NA     <NA>           NA      <NA>
# 25    HN038            1    Yellow
# 151 USB-217            1    Yellow
# 237 USB-449            1    Yellow
# 296 USB-642            1    Yellow
srr["HN038"]
#                                 HN038 
# "SRR12191331;SRR12190723;SRR12190734" 
srr["USB-217"]
#                   USB-217 
# "SRR12336056;SRR12336057" 
srr["USB-449"]
#       USB-449 
# "SRR12366234" 
srr["USB-642"]
#       USB-642 
# "SRR12366139" 
hilum_color_all_T[hilum_color_all_T$phenotype == "Brown" & hilum_color_all_T$haplotype_id == 2, -2]
#    sample haplotype_id phenotype
# 35  HN057            2     Brown
srr["HN057"]
#                                 HN057 
# "SRR12190489;SRR12190490;SRR12190491" 


