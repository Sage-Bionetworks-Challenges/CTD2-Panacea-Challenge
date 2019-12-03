library(tidyverse)

hgnc <- read_tsv("hgnc_ids.txt") %>% set_names("hgnc_id", "symbol", "alt_symbol")

universe <- read_csv("possible_targets.csv") %>% set_names(c("x", "target")) %>% select(-x)

klaeger_targets <- synGet("syn21036465")$path %>% 
  readr::read_csv() %>% 
  rename(cmpd = X1) %>% 
  gather(target, pkd, -cmpd) %>% 
  select(target) %>% 
  separate_rows(target, sep = '\\.') %>% 
  distinct() %>% 
  write_csv("targets_for_scoring.csv")

universe <- bind_rows(universe, klaeger_targets) %>% 
  distinct()

current <- universe %>% filter(target %in% hgnc$symbol) %>% 
  left_join(hgnc %>% 
              select(hgnc_id, symbol) %>% 
              distinct(), by = c('target'='symbol')) %>% 
  filter(hgnc_id != 'HGNC:10811') ##this is a duplicate of SGK2

alt <- universe %>% filter(!target %in% current$target) %>% 
  left_join(hgnc %>% 
              select(hgnc_id, alt_symbol) %>% 
              separate_rows(alt_symbol, sep = ",") %>% 
              mutate(symbol = trimws(alt_symbol)) %>% 
              select(-alt_symbol) %>% 
              distinct(), by = c('target'='symbol'))

universe_map <- bind_rows(current, alt) %>%
  distinct() %>% 
  arrange(target) %>% 
  write_csv("target_list.csv")

id_map <- synGet("syn21197825")$path %>% 
  readr::read_csv()

cols <- c(rep(NA, length(cmpd_ids$cmpd_id))) %>% 
  set_names(cmpd_ids$cmpd_id)

template <- universe_map %>% 
  select(-hgnc_id)

for(i in names(cols)){
  template <- template %>% add_column(!!i := NA)
}

write_csv(template, 'template.csv')






