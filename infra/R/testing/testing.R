source('validation_scoring.R')

#leaderboard

template <- read_csv('template.csv')
gold <- read_csv("panacea_gold_standard.csv")

random_pred <- template %>% 
  mutate_at(vars(-target), ~ runif(length(.))) %>% 
  write_csv('random_pred.csv')

i <- 'cmpd_YM'

perfect_pred <- template

for(i in colnames(perfect_pred)[-1]){
  
  targets <- template$target
  foo_g <- filter(gold, cmpd_id == i)
  foo_p <- select(perfect_pred, target, i) %>% 
    rename("pred" = i)
  
  foo_p <- foo_p %>% 
    rowwise() %>% 
    mutate(pred = case_when(target %in% foo_g$target ~ runif(1, min = 0.9, max = 1),
                              !target %in% foo_g$target ~ runif(1, min = 0, max = 0.9))) %>% 
    set_names(c("target", eval(i)))
                  
  perfect_pred <- perfect_pred %>% 
    select(-i) %>% 
    left_join(foo_p)
  
}

write_csv(perfect_pred, 'perfect_pred.csv')

good_pred <- template

for(i in colnames(good_pred)[-1]){
  
  targets <- template$target
  foo_g <- filter(gold, cmpd_id == i)
  foo_p <- select(good_pred, target, i) %>% 
    rename("pred" = i)
  
  foo_p <- foo_p %>% 
    rowwise() %>% 
    mutate(pred = case_when(target %in% foo_g$target ~ runif(1, min = 0.85, max = 1),
                            !target %in% foo_g$target ~ runif(1, min = 0, max = 0.95))) %>% 
    set_names(c("target", eval(i)))
  
  good_pred <- good_pred %>% 
    select(-i) %>% 
    left_join(foo_p)
  
}

write_csv(good_pred, 'good_pred.csv')

ok_pred <- template

for(i in colnames(ok_pred)[-1]){
  
  targets <- template$target
  foo_g <- filter(gold, cmpd_id == i)
  foo_p <- select(ok_pred, target, i) %>% 
    rename("pred" = i)
  
  foo_p <- foo_p %>% 
    rowwise() %>% 
    mutate(pred = case_when(target %in% foo_g$target ~ runif(1, min = 0.7, max = 1),
                            !target %in% foo_g$target ~ runif(1, min = 0, max = 0.95))) %>% 
    set_names(c("target", eval(i)))
  
  ok_pred <- ok_pred %>% 
    select(-i) %>% 
    left_join(foo_p)
  
}

write_csv(ok_pred, 'ok_pred.csv')


sc1_null_final <- synGet("syn21341654")$path
sc2_null_final <- synGet("syn21341657")$path
sc1_null_lead <- synGet('syn21341653')$path
sc2_null_lead <- synGet('syn21341655')$path

foo <- tribble(~sc1, ~sc2, ~null)

for(i in 1:10){
  bar <- score(prediction_path = "random_pred.csv",
               gold_path = "panacea_gold_standard.csv",
               null_model_path_sc1 = sc1_null_lead,
               null_model_path_sc2 = sc2_null_lead,
               round = "leaderboard")
  
  bar <- tribble(~sc1, ~sc2, ~null,  bar[1], bar[2], "leaderboard")
  
  print(bar)
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "random_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null,
              null_model_path_sc2 = sc2_null,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null,  bar2[1], bar2[2], "final")
foo_plot_1 <- foo %>% bind_rows(bar) %>% gather(sc, score, -null) %>% mutate(pred = "random")

foo <- tribble(~sc1, ~sc2, ~null)

for(i in 1:10){
  bar <- score(prediction_path = "ok_pred.csv",
               gold_path = "panacea_gold_standard.csv",
               null_model_path_sc1 = sc1_null_lead,
               null_model_path_sc2 = sc2_null_lead,
               round = "leaderboard")
  
  bar <- tribble(~sc1, ~sc2, ~null,  bar[1], bar[2], "leaderboard")
  
  print(bar)
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "ok_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null,
              null_model_path_sc2 = sc2_null,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null,  bar2[1], bar2[2], "final")

foo_plot_2 <- foo %>% bind_rows(bar) %>% gather(sc, score, -null) %>% mutate(pred = "ok")



foo <- tribble(~sc1, ~sc2, ~null)

for(i in 1:10){
  bar <- score(prediction_path = "good_pred.csv",
               gold_path = "panacea_gold_standard.csv",
               null_model_path_sc1 = sc1_null_lead,
               null_model_path_sc2 = sc2_null_lead,
               round = "leaderboard")
  
  bar <- tribble(~sc1, ~sc2, ~null,  bar[1], bar[2], "leaderboard")
  
  print(bar)
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "good_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null,
              null_model_path_sc2 = sc2_null,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null,  bar2[1], bar2[2], "final")

foo_plot_3 <- foo %>% bind_rows(bar) %>% gather(sc, score, -null) %>% mutate(pred = "good")


foo <- tribble(~sc1, ~sc2, ~null)

for(i in 1:10){
  bar <- score(prediction_path = "perfect_pred.csv",
               gold_path = "panacea_gold_standard.csv",
               null_model_path_sc1 = sc1_null_lead,
               null_model_path_sc2 = sc2_null_lead,
               round = "leaderboard")
  
  bar <- tribble(~sc1, ~sc2, ~null,  bar[1], bar[2], "leaderboard")
  
  print(bar)
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "perfect_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null,
              null_model_path_sc2 = sc2_null,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null,  bar2[1], bar2[2], "final")

foo_plot_4 <- foo %>% bind_rows(bar) %>% gather(sc, score, -null) %>% mutate(pred = "perfect")

foo_all <- bind_rows(foo_plot_1, foo_plot_2, foo_plot_3, foo_plot_4)

ggplot(foo_all %>% filter(sc == "sc1")) +
  geom_point(aes(x = pred, y = score, color = null), stat = "identity")


ggplot(foo_all %>% filter(sc == "sc2")) +
  geom_point(aes(x = pred, y = score, color = null), stat = "identity")
