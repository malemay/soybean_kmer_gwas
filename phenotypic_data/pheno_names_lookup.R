# Creating a set of lookup tables that will allow simplifying the
# names of the phenotypes for used in contigency tables in the
# manuscript
stem_termination <- c("N – Indeterminate (stem tapering gradually toward tip)" = "Indeterminate",
		      "D – Determinate (stem abruptly terminating)" = "Determinate",
		      "S – Semi-determinate (intermediate between determinate and indeterminate)" = "Semi-determinate")

pubescence_color <- c("G – Gray" = "Gray",
		      "T – Tawny" = "Tawny",
		      "Ng – Near gray" = "Near gray",
		      "- – Indicates curly PUB_FORM or glabrous PUB_DENSITY" = "Curly or glabrous",
		      "Lt – Light tawny" = "Light tawny")

pubescence_form <- c("E – Erect on leaf surface" = "Erect",
		     "Sa – Semi-appressed (between erect and appressed)" = "Semi-appressed",
		     "C – Curly (twisted and appressed)" = "Curly",
		     "A – Appressed (most hairs flat on leaf surface)" = "Appressed")

pubescence_density <- c("N – Normal density" = "Normal",
			"Ssp – Semi-sparse (slightly reduced density, especially on pulvinus)" = "Semi-sparse",
			"Dn – Dense (increased density most noticable on stem)" = "Dense",
			"Sp – Sparse (greatly reduced density, most noticable on stem)" = "Sparse",
			"G – Glabrous (no pubescense)" = "Glabrous",
			"Ssp – Semi-sparse (slightly reduced density, especially on pulvinus), N – Normal density" = "Semi-sparse/normal")

seed_coat_luster <- c("S – Shiny" = "Shiny",
		      "D – Dull" = "Dull",
		      "I – Intermediate (between dull and shiny)" = "Intermediate",
		      "B – Bloom (heavy coating of powdery substance adhering to the seed coat)" = "Bloom",
		      "I – Intermediate (between dull and shiny), B – Bloom (heavy coating of powdery substance adhering to the seed coat)" = "Intermediate/Bloom",
		      "S – Shiny, I – Intermediate (between dull and shiny)" = "Shiny/intermediate",
		      "Lb – Light bloom (slight bloom)" = "Light bloom")

pod_color <- c("Br – Brown" = "Brown",
	       "Tn – Tan" = "Tan",
	       "Dbr – Dark brown" = "Dark brown",
	       "Lbr – Light brown" = "Light brown",
	       "Bl – Black" = "Black",
	       "H – Heterogeneous mix of 2 or more colors" = "Heterogeneous")

flower_color <- c("W – White" = "White",
		  "P – Purple" = "Purple",
		  "Lp – Light purple" = "Light purple",
		  "Dp – Dark purple" = "Dark purple")

seed_coat_color <- c("Y – Yellow" = "Yellow",
		     "Br – Brown" = "Brown",
		     "Bl – Black" = "Black",
		     "Gn – Green" = "Green",
		     "Rbr – Reddish brown" = "Reddish brown",
		     "Gnbr – Greenish brown" = "Greenish brown",
		     "Lgn – Light green" = "Light green",
		     "G – Gray" = "Gray",
		     "Ggn – Grayish green" = "Grayish green")

hilum_color <- c("Bf – Buff" = "Buff",
		 "Ib – Imperfect black" = "Imperfect black",
		 "Y – Yellow" = "Yellow",
		 "Tn – Tan" = "Tan",
		 "Ig – Imperfect gray" = "Imperfect gray",
		 "Dg – Dark gray" = "Dark gray",
		 "Br – Brown" = "Brown",
		 "Bl – Black" = "Black",
		 "Brbl – Brown w/black" = "Brown with black",
		 "Gn – Green" = "Green",
		 "G – Gray" = "Gray",
		 "Lbf – Light buff" = "Light buff",
		 "Rbr – Reddish brown" = "Reddish brown",
		 "Gnbr – Greenish-brown" = "Greenish brown",
		 "Brbl – Brown w/black, Bl – Black" = "Brown with black/black",
		 "Y – Yellow, Bf – Buff" = "Yellow/Buff",
		 "Dbf – Dark buff" = "Dark buff",
		 "Blbr – Black hilum with brown outer ring" = "Black hilum with brown outer ring",
		 "Rbf – Reddish buff" = "Reddish buff",
		 "Dib – Dark imperfect black" = "Dark imperfect black")

maturity_group <- c("III – MATURITY GROUP III" = "III",
		    "IV – MATURITY GROUP IV" = "IV",
		    "V – MATURITY GROUP V" = "V",
		    "II – MATURITY GROUP II" = "II",
		    "III – MATURITY GROUP III, IV – MATURITY GROUP IV" = "III/IV")

pheno_names_lookup <- list("STEM.TERMINATION.TYPE" = stem_termination,
			   "PUBESCENCE.COLOR" = pubescence_color,
			   "PUBESCENCE.FORM" = pubescence_form,
			   "PUBESCENCE.DENSITY" = pubescence_density,
			   "SEED.COAT.LUSTER" = seed_coat_luster,
			   "POD.COLOR" = pod_color,
			   "FLOWER.COLOR" = flower_color,
			   "SEED.COAT.COLOR" = seed_coat_color,
			   "HILUM.COLOR" = hilum_color,
			   "MATURITY.GROUP" = maturity_group)

# Saving to file for use in scripts
saveRDS(pheno_names_lookup, file = "phenotypic_data/pheno_names_lookup.rds")

