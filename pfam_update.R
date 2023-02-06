library(tidyverse) 

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

pfam <- readRDS("inst/extdata/pfamDataBioMart.rds")
names(pfam)

kinase <- pfam %>%
  filter(grepl("kinase", NAME)) %>%
  select(pfam_id, NAME) %>%
  unique()

pfam2 <- pfam %>%
  group_by(hgnc_symbol, pfam_id, chromosome_name, gene_start, gene_end, strand, NAME) %>%
  summarise(domains = str_c(domain_start, domain_end, sep = "-", collapse = "; ")) %>%
  filter(grepl(";", domains)) %>%
  filter(pfam_id %in% kinase$pfam_id) %>%
  mutate(domain_start_to_remove = case_when(hgnc_symbol == "NTRK1" ~ "156841034",
                                            TRUE ~ NA_character_))
