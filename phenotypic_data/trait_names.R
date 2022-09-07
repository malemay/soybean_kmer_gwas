# A script that creates a table with the correspondence between
# trait names in the USDA GRIN database and the ones I used
trait_names <- c("STEM.TERMINATION.TYPE" = "stem_termination",
		 "PUBESCENCE.COLOR" = "pubescence_color",
		 "PROTEIN" = "protein",
		 "PUBESCENCE.FORM" = "pubescence_form",
		 "OIL" = "oil",
		 "PUBESCENCE.DENSITY" = "pubescence_density",
		 "SEED.COAT.LUSTER" = "seed_coat_luster",
		 "POD.COLOR" = "pod_color",
		 "FLOWER.COLOR" = "flower_color",
		 "SEED.COAT.COLOR" = "seed_coat_color",
		 "HILUM.COLOR" = "hilum_color",
		 "MATURITY.GROUP" = "maturity_group")

# Saving the table as an .rds file for retrieval by R scripts
saveRDS(trait_names, file = "phenotypic_data/trait_names.rds")

