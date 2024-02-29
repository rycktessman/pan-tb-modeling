library(dplyr)
library(tidyr)
library(data.table)

path_out <- "output/"
analyses <- c("main", "resistN", "dstN", "worse_adhereLTFU", "better_adhereLTFU", "worse_Pan", "short_RS")  
analysis <- analyses[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

#load main output
out_comb <- fread(paste0(path_out, "final_outputs", "_", analysis, ".csv"))

#add total notifications per country so we can estimate absolute outcomes
out_comb <- out_comb %>% mutate(notif_abs=case_when(country=="India"~1862969, #new_ep + new_labconf + new_clindx for 2021
                                                    country=="South Africa"~156699,
                                                    country=="Philippines"~287104)) 

#calculate additional outcomes
out_comb <- out_comb %>% mutate(pop = pop_none + pop_R + pop_N + pop_RN,
                                prop_R_all=(pop_R + pop_RN)/pop,
                                prop_N_RS=pop_N/(pop_N + pop_none),
                                prop_N_RR=pop_RN/(pop_R + pop_RN),
                                deaths=deaths_primary+
                                  deaths_secondary+
                                  deaths_retx,
                                cures=cures_none + cures_R + cures_N + cures_RN,
                                cures_ratio=cures/pop,
                                cures_ratio_none=cures_none/pop_none,
                                cures_ratio_R=cures_R/pop_R,
                                cures_ratio_N=cures_N/pop_N,
                                cures_ratio_RN=cures_RN/pop_RN,
                                cures_init = cures_init_none + cures_init_R + cures_init_N + cures_init_RN,
                                cures_init_ratio=cures_init/pop,
                                cures_init_ratio_none=cures_init_none/pop_none,
                                cures_init_ratio_R=cures_init_R/pop_R,
                                cures_init_ratio_N=cures_init_N/pop_N,
                                cures_init_ratio_RN=cures_init_RN/pop_RN,
                                deaths_ratio = deaths/pop,
                                secondary_ratio = secondary/pop,
                                cures_abs=cures*notif_abs,
                                cures_init_abs=cures_init*notif_abs,
                                noncures_abs=(1-cures)*notif_abs,
                                noncures_init_abs=(1-cures_init)*notif_abs,
                                deaths_abs=deaths*notif_abs,
                                deaths_primary_abs=deaths_primary*notif_abs,
                                secondary_abs=secondary*notif_abs)
out_comb <- out_comb %>% ungroup() %>% group_by(year, country, sim) %>%
  mutate(inc_cures=cures-cures[scenario=="SQ"],
         inc_cures_ratio=cures_ratio-cures_ratio[scenario=="SQ"],
         inc_cures_init=cures_init-cures_init[scenario=="SQ"],
         inc_cures_init_ratio=cures_init_ratio-cures_init_ratio[scenario=="SQ"],
         inc_deaths=deaths-deaths[scenario=="SQ"],
         inc_deaths_ratio=deaths_ratio-deaths_ratio[scenario=="SQ"],
         inc_secondary=secondary-secondary[scenario=="SQ"],
         inc_secondary_ratio=secondary_ratio-secondary_ratio[scenario=="SQ"])

out_sum_yr <- out_comb %>% ungroup()  %>% group_by(year, country, scenario, scenario_lab) %>%
  select(-sim) %>% summarise_all(list("mean"=~mean(.), "lb"=~quantile(., p=0.025),
                                      "ub"=~quantile(., p=0.975)))

#version that has totals across years
out_comb_sum <- out_comb %>% ungroup() %>% select(-c(year, cdr, cfr, cfr_notx, disc_life_exp)) %>% 
  group_by(sim, country, scenario, scenario_lab) %>%
  summarise_all(sum)  %>% 
  #need to recalculate ratios because they can't be summed
  mutate(cures_ratio=cures/pop,
         noncures_ratio=1-cures_ratio,
         cures_ratio_none=cures_none/pop_none,
         cures_ratio_R=cures_R/pop_R,
         cures_ratio_N=cures_N/pop_N,
         cures_ratio_RN=cures_RN/pop_RN,
         cures_init_ratio=cures_init/pop,
         noncures_init_ratio=1-cures_init_ratio,
         cures_init_ratio_none=cures_init_none/pop_none,
         cures_init_ratio_R=cures_init_R/pop_R,
         cures_init_ratio_N=cures_init_N/pop_N,
         cures_init_ratio_RN=cures_init_RN/pop_RN,
         prop_R_all=(pop_R + pop_RN)/pop,
         prop_N_RS=pop_N/(pop_N + pop_none),
         prop_N_RR=pop_RN/(pop_R + pop_RN),
         rel_cohort=rel_cohort/10,
         deaths_ratio=deaths/pop,
         secondary_ratio=secondary/pop)
out_comb_sum <- out_comb_sum %>% ungroup() %>% group_by(country, sim) %>%
  mutate(inc_cures=cures-cures[scenario=="SQ"],
         inc_cures_ratio=cures_ratio-cures_ratio[scenario=="SQ"],
         inc_cures_abs=cures_abs-cures_abs[scenario=="SQ"],
         inc_cures_pct_change=inc_cures_abs/cures_abs[scenario=="SQ"],
         inc_noncures_abs=noncures_abs-noncures_abs[scenario=="SQ"],
         inc_noncures_pct_change=inc_noncures_abs/noncures_abs[scenario=="SQ"],
         inc_cures_init=cures_init-cures_init[scenario=="SQ"],
         inc_cures_init_ratio=cures_init_ratio-cures_init_ratio[scenario=="SQ"],
         inc_cures_init_abs=cures_init_abs-cures_init_abs[scenario=="SQ"],
         inc_cures_init_pct_change=inc_cures_init_abs/cures_init_abs[scenario=="SQ"],
         inc_noncures_init_abs=noncures_init_abs-noncures_init_abs[scenario=="SQ"],
         inc_noncures_init_pct_change=inc_noncures_init_abs/noncures_init_abs[scenario=="SQ"],
         inc_deaths=deaths-deaths[scenario=="SQ"],
         inc_deaths_ratio=deaths_ratio-deaths_ratio[scenario=="SQ"],
         inc_deaths_abs=deaths_abs-deaths_abs[scenario=="SQ"],
         inc_deaths_pct_change=inc_deaths_abs/deaths_abs[scenario=="SQ"],
         inc_secondary=secondary-secondary[scenario=="SQ"],
         inc_secondary_ratio=secondary_ratio-secondary_ratio[scenario=="SQ"],
         inc_secondary_abs=secondary_abs-secondary_abs[scenario=="SQ"],
         inc_secondary_pct_change=inc_secondary_abs/secondary_abs[scenario=="SQ"])
out_comb_sum <- out_comb_sum %>% 
  mutate(cost_nondrug_tot = cost_outpatient + cost_inpatient + cost_labs +
           cost_aes + cost_support + cost_oop + cost_indirect)

out_sum <-  out_comb_sum %>% ungroup() %>% group_by(country, scenario, scenario_lab) %>%
  select(-sim) %>% summarise_all(list("mean"=~mean(.), "lb"=~quantile(., p=0.025),
                                      "ub"=~quantile(., p=0.975)))

print("Out sum done")

#COSTS VERSION
costs <- out_comb %>% select(country, scenario, scenario_lab, sim, year, notif_abs, starts_with("cost")) %>%
  select(-c(cost_retx_hs, cost_secondary_hs))
costs <- pivot_longer(costs, cols=starts_with("cost"), names_sep="_", values_to="cost",
                      names_to=c("drop", "cost_type")) %>% select(-drop)

costs_hs <- out_comb %>% select(country, scenario, scenario_lab, sim, year, starts_with("cost")) %>%
  select(-c(cost_oop, cost_indirect, cost_retx, cost_secondary)) %>%
  rename(cost_retx=cost_retx_hs, cost_secondary=cost_secondary_hs)
costs_hs <- pivot_longer(costs_hs, cols=starts_with("cost"), names_sep="_", values_to="cost_hs",
                         names_to=c("drop", "cost_type")) %>% select(-drop)

print("Cost pivoting done")

#combine HS and societal perspective costs
costs <- left_join(costs, costs_hs, by=c("country", "scenario", "scenario_lab",
                                         "sim", "year", "cost_type"))
costs <- costs %>% mutate(cost_hs=if_else(is.na(cost_hs) & cost_type %in% c("oop", "indirect"), 0, cost_hs))


#costs of drugs are currently meaningless in the pan-TB scenario - set to 0 or used to quantify person-time on tx
costs <- costs %>% mutate(cost=if_else(cost_type=="drugs" & scenario %in% c("Pan", "Pan_inject"),
                                       0, cost),
                          cost_hs=if_else(cost_type=="drugs" & scenario %in% c("Pan", "Pan_inject"),
                                          0, cost_hs))

costs <- costs %>% group_by(country, year, scenario, sim) %>% 
  mutate(cost_tot=sum(cost),
         cost_tot_nodrugs=sum(cost[cost_type!="drugs"]),
         cost_hs_tot=sum(cost_hs),
         cost_short_hs_tot=sum(cost_hs[!(cost_type %in% c("secondary", "retx"))]),
         cost_short_tot=sum(cost[!cost_type %in% c("secondary", "retx")]))
print("Cost totals done")

costs <- costs %>% mutate(cost_tot_nodrugs_abs=cost_tot_nodrugs*notif_abs*10)

#incremental costs
costs <- costs %>% ungroup() %>% group_by(country, year, sim) %>%
  mutate(inc_cost_tot_nodrugs_abs=cost_tot_nodrugs_abs[scenario=="SQ"]-cost_tot_nodrugs_abs,
         cost_tot_pct_change=(cost_tot_nodrugs-cost_tot_nodrugs[scenario=="SQ"])/cost_tot_nodrugs[scenario=="SQ"])
print("Incremental costs done")

costs_sum <- costs %>% group_by(country, scenario, scenario_lab, cost_type) %>% 
  summarise(cost_mean=mean(cost), cost_lb=quantile(cost, p=0.025), 
            cost_ub=quantile(cost, p=0.975), cost_tot_mean=mean(cost_tot), 
            cost_tot_lb=quantile(cost_tot, p=0.025), cost_tot_ub=quantile(cost_tot, p=0.975),
            cost_hs_mean=mean(cost_hs), cost_hs_lb=quantile(cost_hs, p=0.025), 
            cost_hs_ub=quantile(cost_hs, p=0.975), cost_hs_tot_mean=mean(cost_hs_tot), 
            cost_hs_tot_lb=quantile(cost_hs_tot, p=0.025), cost_hs_tot_ub=quantile(cost_hs_tot, p=0.975),
            cost_short_hs_tot_mean=mean(cost_short_hs_tot), cost_short_hs_tot_lb=quantile(cost_short_hs_tot, p=0.025), 
            cost_short_hs_tot_ub=quantile(cost_short_hs_tot, p=0.975), 
            cost_short_tot_mean=mean(cost_short_tot), cost_short_tot_lb=quantile(cost_short_tot, p=0.025), 
            cost_short_tot_ub=quantile(cost_short_tot, p=0.975),
            cost_tot_pct_change_mean=mean(cost_tot_pct_change), 
            cost_tot_pct_change_lb=quantile(cost_tot_pct_change, p=0.025), 
            cost_tot_pct_change_ub=quantile(cost_tot_pct_change, p=0.975),
            cost_tot_nodrugs_abs_mean=mean(cost_tot_nodrugs_abs),
            cost_tot_nodrugs_abs_lb=quantile(cost_tot_nodrugs_abs, p=0.025),
            cost_tot_nodrugs_abs_ub=quantile(cost_tot_nodrugs_abs, p=0.975),
            inc_cost_tot_nodrugs_abs_mean=mean(inc_cost_tot_nodrugs_abs),
            inc_cost_tot_nodrugs_abs_lb=quantile(inc_cost_tot_nodrugs_abs, p=0.025),
            inc_cost_tot_nodrugs_abs_ub=quantile(inc_cost_tot_nodrugs_abs, p=0.975))
print("Costs sum done")

#DALYS VERSION
dalys <- out_comb %>% select(country, scenario, scenario_lab, sim, year, starts_with("dalys"))
dalys <- pivot_longer(dalys, cols=starts_with("dalys"), names_sep="_", values_to="dalys",
                      names_to=c("drop", "daly_type")) %>% select(-drop)
dalys <- dalys %>% group_by(country, year, scenario, sim) %>% 
  mutate(dalys_tot=sum(dalys))
dalys_sum <- dalys %>% group_by(country, scenario, scenario_lab, daly_type) %>% 
  summarise(dalys_mean=mean(dalys), dalys_lb=quantile(dalys, p=0.025), 
            dalys_ub=quantile(dalys, p=0.975), dalys_tot_mean=mean(dalys_tot), 
            dalys_tot_lb=quantile(dalys_tot, p=0.025), dalys_tot_ub=quantile(dalys_tot, p=0.975))

print("DALYs sum done")

# SAVE ALL TO FILE
write.csv(out_sum_yr, file=paste0(path_out, "out_sum_yr_", analysis, ".csv"), row.names=F)
write.csv(out_sum, file=paste0(path_out, "out_sum_", analysis, ".csv"), row.names=F)
write.csv(costs_sum, file=paste0(path_out, "costs_sum_", analysis, ".csv"), row.names=F)
write.csv(dalys_sum, file=paste0(path_out, "dalys_sum_", analysis, ".csv"), row.names=F)