
pspa_gc <- read_tsv("data/shiny/pspa-gc_lin_counts.txt")

pspa_gc.count <- pspa_gc %>%
  group_by(GenContext.norep) %>%
  summarise(totalcount=sum(count))

pspa_cum <- left_join(pspa_gc, pspa_gc.count, by="GenContext.norep")

pspa_cum %>%
  dplyr::filter(totalcount>=15)
