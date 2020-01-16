source('../validation_scoring_subset.R')
library(synapser)
library(ggplot2)
synLogin()
#leaderboard

template <- read_csv('template.csv')
gold <- read_csv("panacea_gold_standard.csv")

# random_pred <- template %>% 
#   mutate_at(vars(-target), ~ runif(length(.))) %>% 
#   write_csv('random_pred.csv')
# 
# i <- 'cmpd_YM'
# 
# perfect_pred <- template
# 
# for(i in colnames(perfect_pred)[-1]){
#   
#   targets <- template$target
#   foo_g <- filter(gold, cmpd_id == i)
#   foo_p <- select(perfect_pred, target, i) %>% 
#     rename("pred" = i)
#   
#   foo_p <- foo_p %>% 
#     rowwise() %>% 
#     mutate(pred = case_when(target %in% foo_g$target ~ runif(1, min = 0.9, max = 1),
#                               !target %in% foo_g$target ~ runif(1, min = 0, max = 0.9))) %>% 
#     set_names(c("target", eval(i)))
#                   
#   perfect_pred <- perfect_pred %>% 
#     select(-i) %>% 
#     left_join(foo_p)
#   
# }
# 
# write_csv(perfect_pred, 'perfect_pred.csv')
# 
# good_pred <- template
# 
# for(i in colnames(good_pred)[-1]){
#   
#   targets <- template$target
#   foo_g <- filter(gold, cmpd_id == i)
#   foo_p <- select(good_pred, target, i) %>% 
#     rename("pred" = i)
#   
#   foo_p <- foo_p %>% 
#     rowwise() %>% 
#     mutate(pred = case_when(target %in% foo_g$target ~ runif(1, min = 0.85, max = 1),
#                             !target %in% foo_g$target ~ runif(1, min = 0, max = 0.95))) %>% 
#     set_names(c("target", eval(i)))
#   
#   good_pred <- good_pred %>% 
#     select(-i) %>% 
#     left_join(foo_p)
#   
# }
# 
# write_csv(good_pred, 'good_pred.csv')
# 
# ok_pred <- template
# 
# for(i in colnames(ok_pred)[-1]){
#   
#   targets <- template$target
#   foo_g <- filter(gold, cmpd_id == i)
#   foo_p <- select(ok_pred, target, i) %>% 
#     rename("pred" = i)
#   
#   foo_p <- foo_p %>% 
#     rowwise() %>% 
#     mutate(pred = case_when(target %in% foo_g$target ~ runif(1, min = 0.7, max = 1),
#                             !target %in% foo_g$target ~ runif(1, min = 0, max = 0.95))) %>% 
#     set_names(c("target", eval(i)))
#   
#   ok_pred <- ok_pred %>% 
#     select(-i) %>% 
#     left_join(foo_p)
#   
# }
# 
# write_csv(ok_pred, 'ok_pred.csv')




sc1_null_final <- synGet("syn21362374")$path
sc2_null_final <- synGet("syn21362377")$path
sc1_null_lead <- synGet('syn21362373')$path
sc2_null_lead <- synGet('syn21362375')$path

foo <- tribble(~sc1, ~sc2, ~null)

for(i in 1:10){
  bar <- score(prediction_path = "random_pred.csv",
               gold_path = "panacea_gold_standard.csv",
               null_model_path_sc1 = sc1_null_lead,
               null_model_path_sc2 = sc2_null_lead,
               round = "leaderboard")
  
  bar <- tribble(~sc1, ~sc2, ~null,  bar[1], bar[2], "leaderboard")
  
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "random_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
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
  
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "ok_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
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
  
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "good_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc1_null_final,
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
  
  foo <- foo %>% bind_rows(bar)
}

bar2 <- score(prediction_path = "perfect_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null,  bar2[1], bar2[2], "final")

foo_plot_4 <- foo %>% bind_rows(bar) %>% gather(sc, score, -null) %>% mutate(pred = "perfect")

foo_all <- bind_rows(foo_plot_1, foo_plot_2, foo_plot_3, foo_plot_4)

ggplot(foo_all %>% filter(sc == "sc1")) +
  geom_point(aes(x = pred, y = score, color = null), stat = "identity")

ggplot(foo_all %>% filter(sc == "sc2")) +
  geom_point(aes(x = pred, y = score, color = null), stat = "identity")



##### Test number of null models. 

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
for(i in 1:10){
  bar <- score(prediction_path = "ok_pred.csv",
               gold_path = "panacea_gold_standard.csv",
               null_model_path_sc1 = sc1_null_lead,
               null_model_path_sc2 = sc2_null_lead,
               round = "leaderboard",
               no_rand = j)
  
  bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
  
  foo <- foo %>% bind_rows(bar)
}
}

bar2 <- score(prediction_path = "ok_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 100)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "ok")


ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 100)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 100),
             aes(yintercept = score), color = 'red')

ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 100)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 100),
             aes(yintercept = score), color = 'red')


##with perfect prediction

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "perfect_pred.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}

bar2 <- score(prediction_path = "perfect_pred.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 500)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "perfect")

ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')

ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')


##Redo with Bence's best submission
## at this point I am really regretting not writing a cleaner parallelized function so that this was 30 lines instead of 500
## hashtag technicaldebt

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "sim_rep.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}


bar2 <- score(prediction_path = "sim_rep.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 100)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "dry_run")


p1 <- ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 100),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')


p2 <- ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 100),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')

##shRNA prediction

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "sim_shrna.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}

bar2 <- score(prediction_path = "sim_shrna.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 500)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "dry_run")

p3 <- ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')

p4 <- ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "sim_rep_0.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}

bar2 <- score(prediction_path = "sim_rep_0.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 100)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "dry_run")

p5<-ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 100),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')


p6<-ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "sim_rep_1.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}

bar2 <- score(prediction_path = "sim_rep_1.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 500)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "dry_run")

p7<-ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')

p8<-ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')

foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "sim_rep_2.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}

bar2 <- score(prediction_path = "sim_rep_2.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 500)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "dry_run")

p9<-ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')

p10<-ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')


foo <- tribble(~sc1, ~sc2, ~null)

for(j in c(1,5,10,15)){
  for(i in 1:25){
    bar <- score(prediction_path = "sim_ctrp.csv",
                 gold_path = "panacea_gold_standard.csv",
                 null_model_path_sc1 = sc1_null_lead,
                 null_model_path_sc2 = sc2_null_lead,
                 round = "leaderboard",
                 no_rand = j)
    
    bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar[1], bar[2], "leaderboard", j)
    
    foo <- foo %>% bind_rows(bar)
  }
}

bar2 <- score(prediction_path = "sim_ctrp.csv",
              gold_path = "panacea_gold_standard.csv",
              null_model_path_sc1 = sc1_null_final,
              null_model_path_sc2 = sc2_null_final,
              round = "final") 

bar <- tribble(~sc1, ~sc2, ~null, ~no_rand,  bar2[1], bar2[2], "final", 500)

foo_plot_testing_null <- foo %>% bind_rows(bar) %>% gather(sc, score, -null, -no_rand) %>% mutate(pred = "dry_run")

p11<-ggplot(foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc1") %>% filter(no_rand == 500),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc1')

p12<-ggplot(foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand != 500)) +
  ggbeeswarm::geom_beeswarm(aes(x = factor(no_rand), y = score, color = no_rand, group= no_rand)) +
  geom_hline(data = foo_plot_testing_null %>% filter(sc == "sc2") %>% filter(no_rand == 100),
             aes(yintercept = score), color = 'red') +
  labs(x = "number of null models", title = 'sc2')

cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
                   labels = c("sim_rep", "sim_rep",
                              "sim_shrna", "sim_shrna",
                              "sim_rep_0", "sim_rep_0",
                              "sim_rep_1", "sim_rep_1",
                              "sim_rep_2", "sim_rep_2",
                              "sim_ctrp", "sim_ctrp"),
                 scale = 0.9)
