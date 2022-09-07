# This script prepares the lookup tables used to prepare phenotypic data
# for GWAS analyses
stem_termination_all_lookup <- c("D – Determinate (stem abruptly terminating)" = 1,
				 "S – Semi-determinate (intermediate between determinate and indeterminate)" = 2,
				 "N – Indeterminate (stem tapering gradually toward tip)" = 3)

stem_termination_sn_lookup <- c("S – Semi-determinate (intermediate between determinate and indeterminate)" = 1,
				"N – Indeterminate (stem tapering gradually toward tip)" = 2)

pubescence_color_all_lookup <- c("G – Gray" = 1,
				"T – Tawny" = 2,
				"Lt – Light tawny" = 3)

pubescence_color_nogray_lookup <- c("T – Tawny" = 1,
				   "Lt – Light tawny" = 2)

pubescence_form_all_lookup <- c("A – Appressed (most hairs flat on leaf surface)" = 1,
				"Sa – Semi-appressed (between erect and appressed)" = 2,
				"E – Erect on leaf surface" = 3)

pubescence_form_noerect_lookup <- c("A – Appressed (most hairs flat on leaf surface)" = 1,
				    "Sa – Semi-appressed (between erect and appressed)" = 2)

pubescence_density_lookup <- c("N – Normal density" = 1,
			       "Ssp – Semi-sparse (slightly reduced density, especially on pulvinus)" = 2)

seed_coat_luster_all_lookup <- c("B – Bloom (heavy coating of powdery substance adhering to the seed coat)" = 1,
				 "D – Dull" = 2,
				 "I – Intermediate (between dull and shiny)" = 3,
				 "Lb – Light bloom (slight bloom)" = 4,
				 "S – Shiny" = 5)

seed_coat_luster_nointermediate_lookup <- c("B – Bloom (heavy coating of powdery substance adhering to the seed coat)" = 1,
					    "D – Dull" = 2,
					    "Lb – Light bloom (slight bloom)" = 3,
					    "S – Shiny" = 4)

seed_coat_luster_dullshiny_lookup <- c("D – Dull" = 1,
				       "S – Shiny" = 2)

pod_color_all_lookup <- c("Tn – Tan" = 1,
			  "Br – Brown" = 2,
			  "Dbr – Dark brown" = 3,
			  "Bl – Black" = 4)

pod_color_blbr_lookup <- c("Br – Brown" = 1,
			   "Bl – Black" = 2)

flower_color_lookup <- c("P – Purple" = 1,
			 "W – White" = 2)

seed_coat_color_all_lookup <- c("Y – Yellow" = 1,
				"Gn – Green" = 2,
				"Br – Brown" = 3,
				"Rbr – Reddish brown" = 4,
				"Bl – Black" = 5)

seed_coat_color_greenyellow_lookup <- c("Y – Yellow" = 1,
					"Gn – Green" = 2)

hilum_color_all_lookup <- c("Y – Yellow" = 1,
			    "G – Gray" = 2,
			    "Bf – Buff" = 3,
			    "Br – Brown" = 4,
			    "Rbr – Reddish brown" = 5,
			    "Ib – Imperfect black" = 6,
			    "Bl – Black" = 7)

hilum_color_blackbrown_lookup <- c("Br – Brown" = 1,
				   "Bl – Black" = 2)

hilum_color_rbr_lookup <- c("Br – Brown" = 1,
			    "Rbr – Reddish brown" = 2)

maturity_group_lookup <- c("II – MATURITY GROUP II" = 1,
			   "III – MATURITY GROUP III" = 2,
			   "IV – MATURITY GROUP IV" = 3,
			   "V – MATURITY GROUP V" = 4)

lookup_tables <- list(stem_termination_all_lookup = stem_termination_all_lookup,
		      stem_termination_sn_lookup = stem_termination_sn_lookup,
		      pubescence_color_all_lookup = pubescence_color_all_lookup,
		      pubescence_color_nogray_lookup = pubescence_color_nogray_lookup, 
		      pubescence_form_all_lookup = pubescence_form_all_lookup, 
		      pubescence_form_noerect_lookup = pubescence_form_noerect_lookup,
		      pubescence_density_lookup = pubescence_density_lookup,
		      seed_coat_luster_all_lookup = seed_coat_luster_all_lookup,
		      seed_coat_luster_nointermediate_lookup = seed_coat_luster_nointermediate_lookup,
		      seed_coat_luster_dullshiny_lookup = seed_coat_luster_dullshiny_lookup,
		      pod_color_all_lookup = pod_color_all_lookup,
		      pod_color_blbr_lookup = pod_color_blbr_lookup,
		      flower_color_lookup = flower_color_lookup,
		      seed_coat_color_all_lookup = seed_coat_color_all_lookup,
		      seed_coat_color_greenyellow_lookup = seed_coat_color_greenyellow_lookup,
		      hilum_color_all_lookup = hilum_color_all_lookup,
		      hilum_color_blackbrown_lookup = hilum_color_blackbrown_lookup,
		      hilum_color_rbr_lookup = hilum_color_rbr_lookup,
		      maturity_group_lookup = maturity_group_lookup)

# Saving the lookup tables in an .rds file for retrieval later on
saveRDS(lookup_tables, file = "phenotypic_data/lookup_tables.rds")

