## ----my setup------------------------------------------------------------
library(tidyverse)
library(broom)



slop <- commandArgs(trailingOnly = TRUE)[1] %>% as.numeric()
window <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric()
rep <- commandArgs(trailingOnly = TRUE)[3] %>% as.numeric()
bin <- commandArgs(trailingOnly = TRUE)[4] %>% as.numeric()
cell <- commandArgs(trailingOnly = TRUE)[5] %>% as.character()
qval <- commandArgs(trailingOnly = TRUE)[6] %>% as.numeric()
motif <- commandArgs(trailingOnly = TRUE)[7] %>% as.numeric()
tst <- commandArgs(trailingOnly = TRUE)[8]
n_sample <- commandArgs(trailingOnly = TRUE)[9] %>% as.numeric()



dir_output_table <- str_glue("./result/table/{tst}/slop{slop}/wlen{window}/cell_{cell}_q{qval}_motif{motif}")
dir_output_pdf <- str_glue("./result/pdf/{tst}/slop{slop}/wlen{window}/cell_{cell}_q{qval}_motif{motif}")

dir.create(dir_output_table,
           recursive = T,
           showWarnings = F)
dir.create(dir_output_pdf,
           recursive = T,
           showWarnings = F)


## ----load_{window * 2 + 1}bp----------------------------------------------------------

df <-
  read_tsv(
    str_glue(
      "/home/kogu/analysis/ner_analysis/result_count/{tst}/wlen{window}/motif_under_peak_{cell}_bl_q{qval}_motif{motif}__rep{rep}b_{bin}.slop{slop}__{window}bp_noCDS_callable__count.txt"
    ),
    col_names = c("chrom", "start", "end", "name", "count")
  ) %>%
  mutate(bin_length = end - start)
df

dfv <-
  read_tsv(
    str_glue(
      "/home/kogu/analysis/ner_analysis/result_count/{tst}/wlen{window}/motif_under_peak_{cell}_bl_q{qval}_motif{motif}__rep{rep}b_{bin}.slop{slop}_v__{window}bp_noCDS_callable__count.txt"
    ),
    col_names = c("chrom", "start", "end", "name", "count")
  ) %>%
  mutate(bin_length = end - start)
dfv

dfp <-
  read_tsv(
    str_glue(
      "/home/kogu/analysis/ner_analysis/result_count/{tst}/wlen{window}/motif_under_peak_{cell}_bl_q{qval}_motif{motif}__rep{rep}b_{bin}.slop{slop}_peri__{window}bp_noCDS_callable__count.txt"
    ),
    col_names = c("chrom", "start", "end", "name", "count")
  ) %>%
  mutate(bin_length = end - start)

dfp

# list of all TF names
df %>% .$name %>% sort() %>% unique()


## ----raw_rate_{window * 2 + 1}bp------------------------------------------------------
df_01 <-
  df %>%
  mutate(rate_by_bin = count / bin_length / n_sample) %>%
  arrange(-rate_by_bin)
df_01

dfv_01 <-
  dfv %>%
  mutate(rate_by_bin = count / bin_length / n_sample) %>%
  arrange(-rate_by_bin)
dfv_01


dfp_01 <-
  dfp %>%
  mutate(rate_by_bin = count / bin_length / n_sample) %>%
  arrange(-rate_by_bin)
dfp_01




df_01 %>% write_excel_csv(
  str_glue(
    "{dir_output_table}/01_rawtable_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
  )
)


dfv_01 %>% write_excel_csv(
  str_glue(
    "{dir_output_table}/01_v_rawtable_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
  )
)


dfp_01 %>% write_excel_csv(
  str_glue(
    "{dir_output_table}/01_p_rawtable_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
  )
)




## ----tfbs_rate_{window * 2 + 1}bp-----------------------------------------------------
df_02 <-
  df %>%
  group_by(name) %>%
  summarise(tfbs_count = sum(count),
            tfbs_length = sum(bin_length)) %>%
  mutate(
    rate_by_tfbs = tfbs_count / tfbs_length / n_sample,
    rate_by_tfbs_per_megabase = rate_by_tfbs * 1e6
  ) %>%
  arrange(-rate_by_tfbs)

df_02
df_02 %>% summary

df_02 %>%
  write_excel_csv(
    str_glue(
      "{dir_output_table}/02_tfbstable_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
    )
  )

df_02 %>%
  ggplot(aes(tfbs_length, tfbs_count)) + geom_point() -> gg
ggsave(
  str_glue(
    "{dir_output_pdf}/02_tfbs_length_vs_count_{window * 2 + 1}bp_rep{rep}b_{bin}.pdf"
  ),
  gg
)



## ----dfv_02_{window * 2 + 1}bp--------------------------------------------------------
dfv_02 <-
  dfv %>%
  group_by(name) %>%
  summarise(tfbs_count = sum(count),
            tfbs_length = sum(bin_length)) %>%
  mutate(
    rate_by_tfbs = tfbs_count / tfbs_length / n_sample,
    rate_by_tfbs_per_megabase = rate_by_tfbs * 1e6
  ) %>%
  arrange(-rate_by_tfbs)

dfv_02
dfv_02 %>% summary

dfv_02 %>%
  write_excel_csv(
    str_glue(
      "{dir_output_table}/02_v_tfbstable_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
    )
  )

dfv_02 %>%
  ggplot(aes(tfbs_length, tfbs_count)) +
  geom_point() -> gg
ggsave(
  str_glue(
    "{dir_output_pdf}/02_v_tfbs_length_vs_count_{window * 2 + 1}bp_rep{rep}b_{bin}.pdf"
  ),
  gg
)



## ----dfp_02_{window * 2 + 1}bp--------------------------------------------------------
dfp_02 <-
  dfp %>%
  group_by(name) %>%
  summarise(tfbs_count = sum(count),
            tfbs_length = sum(bin_length)) %>%
  mutate(
    rate_by_tfbs = tfbs_count / tfbs_length / n_sample,
    rate_by_tfbs_per_megabase = rate_by_tfbs * 1e6
  ) %>%
  arrange(-rate_by_tfbs)

dfp_02
dfp_02 %>% summary

dfp_02 %>%
  write_excel_csv(
    str_glue(
      "{dir_output_table}/02_p_tfbstable_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
    )
  )

dfp_02 %>%
  ggplot(aes(tfbs_length, tfbs_count)) +
  geom_point() -> gg
ggsave(
  str_glue(
    "{dir_output_pdf}/02_p_tfbs_length_vs_count_{window * 2 + 1}bp_rep{rep}b_{bin}.pdf"
  ),
  gg
)



## ----join_{window * 2 + 1}bp----------------------------------------------------------
df_02 %>% nrow
dfv_02 %>% nrow
dfp_02 %>% nrow

df_02
dfv_02
dfp_02


dfj <-
  full_join(df_02,
            dfv_02,
            by = "name",
            suffix = c("", "_v")) %>%
  full_join(dfp_02, by = "name", suffix = c("", "_p"))

# replace na of dataframe
dfj[is.na(dfj)] <- 0

dfj <- dfj %>%
  mutate(
    ratio_of_mut_rate_NEAR_per_DISTANT = rate_by_tfbs_per_megabase / rate_by_tfbs_per_megabase_v,
    ratio_of_mut_rate_NEAR_per_PERIPHERY = rate_by_tfbs_per_megabase / rate_by_tfbs_per_megabase_p,
    n_sample = n_sample
  )



# statistical test
dfj
dfj <- dfj %>%
  mutate(
    fisher_p_near_distant =
      pmap_dbl(list(
        a = tfbs_count,
        b = n_sample * tfbs_length - tfbs_count,
        c = tfbs_count_v,
        d = n_sample * tfbs_length_v - tfbs_count_v
      ),
      function(a, b, c, d) {
        fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
      }),
    fisher_q_bh_near_distant = p.adjust(fisher_p_near_distant, method = "BH")
  ) %>%
  mutate(
    fisher_p_near_periphery =
      pmap_dbl(list(
        a = tfbs_count,
        b = n_sample * tfbs_length - tfbs_count,
        c = tfbs_count_p,
        d = n_sample * tfbs_length_p - tfbs_count_p
      ),
      function(a, b, c, d) {
        fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
      }),
    fisher_q_bh_near_periphery = p.adjust(fisher_p_near_periphery, method = "BH")
  ) %>%
  return()



dfj %>%
  write_excel_csv(
    str_glue(
      "{dir_output_table}/03_fulljoin_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
    )
  )





gg <- ggplot(dfj) +
  geom_col(aes(".Near", rate_by_tfbs_per_megabase)) +
  geom_col(aes("Distant", rate_by_tfbs_per_megabase_v),
              fill = "red") +
  geom_col(aes("Periphery", rate_by_tfbs_per_megabase_p),
              fill = "blue") +
  ylab("Mutation rate in TFBS(re-defined)") +
  ## admiting some dots will go over upper limit -----------
ylim(0, 40) +
  ##--------------------------------------------------------
ggtitle(str_glue("Window:{window * 2 + 1}sbp, rep{rep}b_{bin}")) +
  NULL

ggsave(
  str_glue(
    "{dir_output_pdf}/03_mut_rate_by_TF_{window * 2 + 1}bp_rep{rep}b_{bin}.pdf"
  ),
  gg
)


dfj %>%
  select(starts_with("rate_by_tfbs_per_megabase")) %>%
  psych::describe() %>%
  rownames_to_column() %>%
  write_excel_csv(
    str_glue(
      "{dir_output_table}/03b_fulljoin_summary_of_rate_{tst}_{window * 2 + 1}bp_rep{rep}b_{bin}.csv"
    )
  )


## ------------------------------------------------------------------------
sessioninfo::session_info()
