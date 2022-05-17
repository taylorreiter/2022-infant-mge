library(dplyr)
library(readr)

cdbg_by_record <- read_csv(snakemake@input[["cdbg_by_record"]], show_col_types = F) %>%
  mutate(total_abundance = total_kmers * mean_abund,
         database_origin = ifelse(grepl("viridae|virales|virus", x = record_name), "chvd", "card")) 

total_kmers <- cdbg_by_record %>%
  group_by(database_origin) %>%
  slice_max(total_kmers, n = 10)

total_abundance <- cdbg_by_record %>%
  group_by(database_origin) %>%
  slice_max(total_abundance, n = 10)

mean_abund <- cdbg_by_record %>%
  group_by(database_origin) %>%
  slice_max(mean_abund, n = 10)

candidate_mges <- bind_rows(total_kmers, total_abundance, mean_abund) %>%
  ungroup() %>%
  distinct() %>%
  select(record_name) 

write_tsv(candidate_mges, snakemake@output[['candidate_mges']], col_names = F)
