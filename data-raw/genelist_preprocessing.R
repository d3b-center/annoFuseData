# Read in the necessary data
library(tidyverse)

ann <- read_tsv("../inst/extdata/genelistreference.txt") %>%
  select(-classification)
onco_tsg <- read_tsv("../inst/extdata/cancerGeneList.tsv") %>%
  dplyr::rename(Gene_Symbol = `Hugo Symbol`) %>%
  mutate(classification = case_when(
    `Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "No" ~ "Oncogene",
    `Is Oncogene` == "No" & `Is Tumor Suppressor Gene` == "Yes" ~ "TumorSuppressorGene",
    `Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "Yes" ~ "Oncogene, TumorSuppressorGene",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(classification))
  
#' Internal preprocessing function
#'
#' This function does some preprocessing and is not exported.
#'
#' @keywords internal

# Function to add classification, update file column, and add missing genes
update_ann_rm <- function(gene_file, onco_tsg) {
  # Filter out the specified file and clean the file and type columns
  gene_file <- gene_file %>%
    dplyr::filter(file != "allOnco_Feb2017.tsv") %>%
    mutate(file = gsub(", allOnco_Feb2017.tsv|allOnco_Feb2017.tsv, ", "", file),
           type = gsub(", Oncogene|Oncogene, ", "", type))
  
  # Create a lookup table from onco_tsg
  classification_lookup <- onco_tsg %>% 
    select(Gene_Symbol, classification)
  
  # Perform a full join to ensure all genes from onco_tsg are included
  combined_data <- full_join(gene_file, classification_lookup, by = c("Gene_Symbol"))
  
  # Update the type and file columns
  combined_data <- combined_data %>%
    rowwise() %>%
    mutate(
      # Add classification to type if not already present
      type = ifelse(
        is.na(type), classification, 
        ifelse(!is.na(classification), 
               paste0(type, ifelse(type == "", "", ", "), classification), 
               type)
      ),
      # Ensure type contains unique classifications
      type = paste(unique(trimws(unlist(strsplit(type, ", ")))), collapse = ", "),
      # Add "OncoKB" to file if not already present
      file = ifelse(is.na(file) & Gene_Symbol %in% classification_lookup$Gene_Symbol, "OncoKBv20240704", 
                    ifelse(!grepl("OncoKBv20240704", file) & Gene_Symbol %in% classification_lookup$Gene_Symbol, 
                           paste0(file, ifelse(file == "", "", ", "), "OncoKBv20240704"), 
                           file))
    ) %>%
    ungroup() %>%
    distinct() # Ensure the rows are distinct
  
  return(combined_data)
}

# Apply the function to ann_rm
ann_rm <- update_ann_rm(ann, onco_tsg) %>%
  group_by(Gene_Symbol, classification) %>%
  mutate(file = case_when(is.na(classification) ~ gsub(", OncoKBv20240704", "", file), 
                          TRUE ~ file)) %>%
  reframe(type = paste(unique(trimws(unlist(strsplit(type, ", ")))), collapse = ", "),
          file = paste(unique(trimws(unlist(strsplit(file, ", ")))), collapse = ", ")) %>%
  select(-classification) %>%
  arrange(Gene_Symbol) %>%
  write_tsv("../inst/extdata/genelistreference.txt")
