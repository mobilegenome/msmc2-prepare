library(tidyverse)

files <- fs::dir_ls("logs/merge_genotypes/")

df <- read_tsv(files, col_names=c("variable", "value"), id="fpath")

df2 <- df |> mutate(sample = str_remove(fpath, "logs/merge_genotypes/"), chr=str_extract(fpath, "chr[0-9ALGE]+")) |> mutate(sample = str_remove(sample, "_chr.*")) |> select(-fpath)

df2 |> pivot_wider(id_cols=c(sample, chr), names_from=variable, values_from=value) |> group_by(sample) |> 
    summarise(frac_replaced = round(sum(replaced_genotypes)/sum(read_lines),3),
              n_inconsistent_ref_alt_gt = sum(inconsistent_ref_alt_genotypes),
              frac_inconsistent_ref_alt_gt = round(sum(inconsistent_ref_alt_genotypes)/sum(read_lines),3),
              n_inconsistent_ref_alt_allels = sum(inconsistent_ref_alt_alleles),
              frac_inconsistent_ref_alt_allels = round(sum(inconsistent_ref_alt_alleles)/sum(read_lines),3),
              ) -> df3

write_tsv(df3, "results/overview.tsv")
