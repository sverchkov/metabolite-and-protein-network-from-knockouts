# Make a t-sne plot of the data

library(tidyr)
library(dplyr)
library(tsne)
library(ggplot2)
library(glue)

# Metadata
molecule_ids <- readRDS('processed-data/molecule-ids.rds')

# Input is the t-statistic table
feature_df <- readRDS('processed-data/feature_t_df.rds')

# Data prep:

# Subsample for testing
active_ids <- (
  if( FALSE ){ # Subsampling switch
    molecule_ids %>%
    group_by(`Molecule Type`) %>%
    sample_frac(0.02) %>%
    ungroup()
  } else {
    molecule_ids
  } ) %>% select(`Molecule ID`, `Molecule Type`)

# High-NA content found in KO 23 (PPTC7)
nona_data <- inner_join(active_ids, feature_df, by="Molecule ID") %>%
  select(-`PPTC7`) %>% drop_na()

molecule_info <- nona_data %>% select(`Molecule ID`, `Molecule Type`)

molecules_matrix <- nona_data %>% select(-`Molecule ID`, -`Molecule Type`) %>% as.matrix()
kos_matrix <- t(molecules_matrix)

ko_info <- rownames(kos_matrix)

for (molecule_perplexity in c(16,32,64,128,256,512,1024)) {
  
  # Run molecule t-sne
  #molecule_perplexity <- 8
  tsne_molecules <- tsne(X=molecules_matrix, k=2, perplexity=molecule_perplexity)
  tsne_molecules_df <- molecule_info %>% mutate(x=tsne_molecules[,1], y=tsne_molecules[,2])
  
  # Plot molecule t-sne
  ggplot(tsne_molecules_df, aes(x=x, y=y, color=`Molecule Type`)) + geom_point()
  
  # Save pdf
  ggsave(glue('results/molecule-tsne-{molecule-perplexity}.pdf'))
}

for (knockout_perplexity in c(8,16,32,64)){
  
  # Run knockout t-sne
  tsne_knockouts <- tsne(X=kos_matrix, k=2, perplexity=knockout_perplexity)
  tsne_knockouts_df <- tibble( `KO`=ko_info, x=tsne_knockouts[,1], y=tsne_knockouts[,2])
  
  # Plot knockout t-sne
  ggplot(tsne_knockouts_df, aes(x=x, y=y, label=KO)) + geom_text()
  
  # Save pdf
  ggsave(glue('results/knockout-tsne-{knockout_perplexity}.pdf'))
}