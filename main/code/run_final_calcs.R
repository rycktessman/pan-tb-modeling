#adjust output from run_model for cohort sizes and secondary cases in estimate_incidence
library(dplyr)

#basic info
analyses <- c("main", "resistN", "dstN", "worse_adhereLTFU", "better_adhereLTFU", "worse_Pan", "short_RS",
              "bad_efficacy", "good_efficacy", "bad_duration", "good_duration", "bad_adherence", "good_adherence",
              "bad_forgiveness", "good_forgiveness", "bad_safety", "good_safety",
              "bad_baseres", "good_baseres", "bad_restrend", "good_restrend")  
analysis <- analyses[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
print(analysis)

countries <- c("india", "southafrica", "philippines")
names(countries) <- c("India", "South Africa", "Philippines")

if(analysis %in% analyses[1:7]) {
  scen_names <- c("Status Quo",
                  "Pan-TB TRP: Oral",
                  "Improved RS Regimen",
                  "Improved RR Regimen",
                  "Improved RR Retention",
                  "Improved Rif DST",
                  "Improved RS & RR Regimens",
                  "Pan-TB TRP: LAI",
                  "Pan-TB Oral, with DST")
} else {
  scen_names <- c("Status Quo",
                  "Pan-TB TRP: Oral",
                  "Pan-TB TRP: LAI")
}

#load parameters
load("data/params.Rda")
load("data/cost_params.Rda") #only for life expectancy, so don't need to use LAI version since it's the same

print(analysis)

#load output from run_model and estimate_incidence
path_out <- "output/"
out_comb <- data.frame()
for(i in names(countries)) {
  load(paste0(path_out, "out_", countries[i], "_", analysis,  ".Rda"))
  scenarios <- names(out_all)
  names(scenarios) <- scen_names
  if(analysis %in% analyses[1:7]) {
    scenarios <- scenarios[c(1, 2, 8, 9, 3:7)] #reorder for factor
  }
  out_all <- bind_rows(out_all, .id="scenario")
  out_all <- out_all %>% mutate(country=i)
  out_comb <- bind_rows(out_comb, out_all)
}

out_comb <- out_comb %>% mutate(scenario_lab=names(scenarios)[match(scenario, scenarios)])
out_comb <- out_comb %>% mutate(scenario_lab=factor(scenario_lab, levels=names(scenarios)))

out_secondary <- read.csv(paste0(path_out, "secondary_out_all", "_", analysis, ".csv")) %>% rename("year"="year_int")

out_comb <- left_join(out_comb, out_secondary)

#merge with relevant parameters: case fatality ratio, case detection ratio, life expectancy
cfr <- lapply(countries, function(x)
  data.frame(sim=1:length(params[[x]][["cfr"]]),
             cfr=params[[x]][["cfr"]]))
cfr <- bind_rows(cfr, .id="country")
out_comb <- left_join(out_comb, cfr)

cdr <- lapply(countries, function(x)
  data.frame(sim=1:length(params[[x]][["cdr"]]),
             cdr=params[[x]][["cdr"]]))
cdr <- bind_rows(cdr, .id="country")
out_comb <- left_join(out_comb, cdr)

disc_life_exp <- lapply(countries, function(x)
  data.frame(sim=1:length(cost_params[[x]][["disc_life_exp"]]),
             disc_life_exp=cost_params[[x]][["disc_life_exp"]]))
disc_life_exp <- bind_rows(disc_life_exp, .id="country")
out_comb <- left_join(out_comb, disc_life_exp)

#estimate deaths, DALYs, and costs among secondary cases, based on avgs among primary cases
#avgs among primary cases are per person, # secondary cases are per 100
out_comb <- out_comb %>% mutate(secondary=secondary/100)
#deaths during and post treatment (primary), during retreatment (retx) and before tx (apply CFR)
out_comb <- out_comb %>% 
  mutate(deaths_secondary=deaths_primary*cdr*secondary + 
           deaths_retx*cdr*secondary + secondary*cfr_notx)
#DALYs
out_comb <- out_comb %>% 
  mutate(dalys_secondary=dalys_tb*secondary + 
           dalys_ae*cdr*secondary +
           dalys_posttb*secondary + 
           dalys_retx*cdr*secondary +
           deaths_secondary*disc_life_exp)
#costs
out_comb <- out_comb %>% 
  mutate(cost_secondary=cost_outpatient*cdr*secondary +
           cost_inpatient*cdr*secondary +
           cost_labs*cdr*secondary +
           cost_aes*cdr*secondary +
           cost_support*cdr*secondary + 
           cost_oop*cdr*secondary +
           cost_indirect*cdr*secondary +
           cost_drugs*cdr*secondary +
           cost_retx*cdr*secondary,
         cost_secondary_hs=cost_outpatient*cdr*secondary +
           cost_inpatient*cdr*secondary +
           cost_labs*cdr*secondary +
           cost_aes*cdr*secondary +
           cost_support*cdr*secondary + 
           cost_drugs*cdr*secondary +
           cost_retx_hs*cdr*secondary)

#now adjust the primary outcomes by relative cohort size under non-status quo scenarios
out_comb <- out_comb %>% group_by(country, sim, year) %>%
  mutate(rel_cohort=cohort/cohort[scenario=="SQ"])
out_comb <- out_comb %>% 
  mutate_at(c("pop_none", "pop_R", "pop_N", "pop_RN", 
              "cures_init_none", "cures_init_R", "cures_init_N", "cures_init_RN", 
              "cures_none", "cures_R", "cures_N", "cures_RN", 
              "failures", "deaths_primary" ,"deaths_retx", "cost_outpatient", "cost_inpatient",
              "cost_labs", "cost_aes", "cost_support", "cost_drugs", "cost_oop", "cost_indirect",
              "cost_retx", "cost_retx_hs", "dalys_ae", "dalys_tb", "dalys_posttb", "dalys_deaths", "dalys_retx",
              "initiate_primary_RS", "initiate_primary_RR", "initiate_primary_P", "initiate_primary_I",
              "initiate_all_RS", "initiate_all_RR", "initiate_all_P", "initiate_all_I"), ~(.x)*rel_cohort)

#save to file
write.csv(out_comb, file=paste0(path_out, "final_outputs", "_", analysis, ".csv"), row.names=F)

