library(tidyverse) 

# Set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

pfam <- readRDS("inst/extdata/pfamDataBioMart.rds")
names(pfam)

kinase <- pfam %>%
  filter("kinase", NAME)) %>%
  select(pfam_id, NAME) %>%
  unique()

updated <- read_tsv("~/Downloads/new_entries.log", col_names = F) %>%
  separate(X1, sep = ":", into = c(NA, "gene"))
intersect(updated$gene, kinase$NAME)

pfam2 <- pfam %>%
  group_by(hgnc_symbol, pfam_id, chromosome_name, gene_start, gene_end, strand, NAME) %>%
  summarise(domains = str_c(domain_start, domain_end, sep = "-", collapse = "; ")) %>%
  filter(";", domains)) %>%
  filter(pfam_id %in% kinase$pfam_id) %>%
  mutate(domain_start_to_remove = case_when(hgnc_symbol == "NTRK1" ~ "156841034",
                                            TRUE ~ NA_character_))

pfam_excel <- as.data.frame(readxl::read_excel("~/Downloads/kinase (2).xlsx")) %>%
  filter(Domain_to_be_updated != "correct" | !is.na(domain_start_to_remove))

dom_rm <- pfam_excel %>%
  filter(!is.na(domain_start_to_remove)) %>%
  dplyr::select(hgnc_symbol, domain_start_to_remove)

  
names(pfam_excel)

new_symbol_pfam <- read_tsv("~/Downloads/pfamDataBiomart_2023_01_01.tsv") %>%
  mutate(remove_row = case_when(hgnc_symbol == "NTRK1" & domain_start == "156841034" ~ "yes",
                                hgnc_symbol == "CAMKV" & domain_start == "49861189" ~ "yes",
                                hgnc_symbol == "CDK1" & domain_start == "60785784" ~ "yes",
                                hgnc_symbol == "CDK5" & domain_start == "151056615" ~ "yes",
                                hgnc_symbol == "EPHA6" & domain_start == "97592736" ~ "yes",
                                hgnc_symbol == "MINK1" & domain_start == "4884468" ~ "yes",
                                hgnc_symbol == "MKNK1" & domain_start == "46576582" ~ "yes",
                                hgnc_symbol == "NEK10" & domain_start == "27192196" ~ "yes",
                                hgnc_symbol == "NTRK1" & domain_start == "156841034" ~ "yes",
                                hgnc_symbol == "PRKCB" & domain_start == "24214672" ~ "yes",
                                hgnc_symbol == "STRADA" & domain_start == "63690289" ~ "yes",
                                hgnc_symbol == "TYK2" & domain_start == "10356586" ~ "yes",
                                TRUE ~ "false"))

