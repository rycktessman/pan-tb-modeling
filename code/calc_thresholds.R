library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

#scenario_use <- "Pan" #Pan or Pan_inject
#analysis <- "main" #main, resistN, dstN, worse_adhereLTFU, better_adhereLTFU, worse_Pan, short_RS  
scenario_use <- Sys.getenv('scenario_use') #Pan or Pan_inject
analyses <- c("main", "resistN", "dstN", "worse_adhereLTFU", "better_adhereLTFU", "worse_Pan", "short_RS")  
analysis <- analyses[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
print(analysis)

load("data/cost_params.Rda") #only used for wastage and drugs schedule, so no need to load LAI params (since same)
path_out <- "output/"
out_all <- fread(paste0(path_out, "final_outputs", "_", analysis, ".csv")) 

#add WTP thresholds
out_all <- out_all %>% mutate(wtp=case_when(country=="India"~430,
                                            country=="Philippines"~1060,
                                            country=="South Africa"~3400))

#adjust initiate_all_P to incorporate 10% inflation in drug costs under the Pan regimen
out_all <- out_all %>% mutate(initiate_all_RS = initiate_primary_RS + 1.1*(initiate_all_RS - initiate_primary_RS),
                              initiate_all_RR = initiate_primary_RR + 1.1*(initiate_all_RR - initiate_primary_RR),
                              initiate_all_P = initiate_primary_P + 1.1*(initiate_all_P - initiate_primary_P),
                              initiate_all_I = initiate_primary_I + 1.1*(initiate_all_I - initiate_primary_I))

#calculate average (with or without discounting - except DALYs are already partly discounted) across years
out_all <- out_all %>% mutate(pop = pop_none + pop_R + pop_N + pop_RN)
out_all <- out_all %>% select(country, scenario, scenario_lab, sim, year, wtp, pop,
                              starts_with("cost"), starts_with("dalys"), starts_with("initiate"))
out_all <- out_all %>% group_by(country, wtp, scenario, scenario_lab, sim) %>%
  mutate_all(list("disc"=~(./(1.03^(year-1)))))
out_all <- out_all %>% group_by(country, wtp, scenario, scenario_lab, sim) %>% 
  summarise_all(mean)
out_all <- out_all %>% select(-c(year_disc, year, pop_disc))

out_all <- out_all %>% mutate(cost_short=cost_outpatient + cost_inpatient + cost_labs + 
                                cost_aes + cost_support + cost_drugs + cost_oop + cost_indirect,
                              cost_short_hs=cost_outpatient + cost_inpatient + cost_labs + 
                                cost_aes + cost_support + cost_drugs,
                              cost_med=cost_short + cost_retx + cost_secondary,
                              cost_med_hs=cost_short + cost_retx_hs + cost_secondary_hs,
                              cost_cea=cost_outpatient_disc + cost_inpatient_disc + cost_labs_disc +
                                cost_aes_disc + cost_support_disc + cost_drugs_disc + cost_oop_disc +
                                cost_indirect_disc + cost_retx_disc + cost_secondary_disc,
                              cost_cea_hs=cost_outpatient_disc + cost_inpatient_disc + cost_labs_disc +
                                cost_aes_disc + cost_support_disc + cost_drugs_disc + 
                                cost_retx_hs_disc + cost_secondary_hs_disc,
                              dalys_cea=dalys_ae_disc + dalys_tb_disc + dalys_posttb_disc +
                                dalys_deaths_disc + dalys_secondary_disc + dalys_retx_disc)

#merge drug costs with parameters
drug_costs <- sapply(names(cost_params), function(x) 
  as.data.frame(t(sapply(1:max(out_all$sim), function(y) rowSums(cost_params[[x]][["drugs_sched"]][[y]])[1:4]))), 
  simplify=F, USE.NAMES=T)
drug_costs <- bind_rows(drug_costs, .id="country")
drug_costs <- drug_costs %>% mutate(country=case_when(country=="india"~"India",
                                                      country=="southafrica"~"South Africa",
                                                      country=="philippines"~"Philippines"))
colnames(drug_costs) <- c("country", "unitcost_RS", "unitcost_RR", "unitcost_P", "unitcost_I")
drug_costs <- drug_costs %>% group_by(country) %>% mutate(sim=1:n())
out_all <- left_join(out_all, drug_costs, by=c("country", "sim"))

wastage <- sapply(names(cost_params), function(x) as.data.frame(cost_params[[x]][["wastage"]]), simplify=F, USE.NAMES=T)
wastage <- bind_rows(wastage, .id="country")
wastage <- wastage %>% mutate(country=case_when(country=="india"~"India",
                                                country=="southafrica"~"South Africa",
                                                country=="philippines"~"Philippines"))
names(wastage) <- c("country", "wastage")
wastage <- wastage %>% group_by(country) %>% mutate(sim=1:n())
out_all <- left_join(out_all, wastage, by=c("country", "sim"))

# THRESHOLD 1: PAN vs. STATUS QUO

# estimate cost-neutral threshold for Pan regimen, vs. status quo with current prices
# first calculate cost of drugs under Pan implied by cost neutrality
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_1=cost_short[scenario=="SQ"] - cost_short[scenario==scenario_use],
         drugs_med_1=cost_med[scenario=="SQ"]-cost_med[scenario==scenario_use],
         drugs_cea_1=wtp*(dalys_cea[scenario=="SQ"] - dalys_cea[scenario==scenario_use]) +
           cost_cea[scenario=="SQ"] - cost_cea[scenario==scenario_use],
         drugs_short_hs_1=cost_short_hs[scenario=="SQ"] - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_1=cost_med_hs[scenario=="SQ"]-cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_1=wtp*(dalys_cea[scenario=="SQ"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_hs[scenario=="SQ"] - cost_cea_hs[scenario==scenario_use])
# next use initiation output to translate the drugs cost into a drug price
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_1=drugs_short_1/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_1=drugs_med_1/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_1=drugs_cea_1/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_1=drugs_short_hs_1/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_1=drugs_med_hs_1/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_1=drugs_cea_hs_1/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))
# note that some retx in Pan scenario take HRZE or Ind. (this is included under Retreatments, not Drugs)

# THRESHOLD 2A: PAN VS. IMPROVED RS (RS cost = HRZE)
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_2a=cost_short[scenario=="RS"] - cost_short[scenario==scenario_use],
         drugs_med_2a=cost_med[scenario=="RS"]-cost_med[scenario==scenario_use],
         drugs_cea_2a=wtp*(dalys_cea[scenario=="RS"] - dalys_cea[scenario==scenario_use]) +
           cost_cea[scenario=="RS"] - cost_cea[scenario==scenario_use],
         drugs_short_hs_2a=cost_short_hs[scenario=="RS"] - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_2a=cost_med_hs[scenario=="RS"]-cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_2a=wtp*(dalys_cea[scenario=="RS"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_hs[scenario=="RS"] - cost_cea_hs[scenario==scenario_use])
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_2a=drugs_short_2a/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_2a=drugs_med_2a/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_2a=drugs_cea_2a/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_2a=drugs_short_hs_2a/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_2a=drugs_med_hs_2a/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_2a=drugs_cea_hs_2a/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))

# THRESHOLD 2B: PAN VS. IMPROVED RS (RS cost = cost-neutral vs. status quo)
#part 1: estimate cost-neutral costs for improved RS regimen vs. status quo
#first estimate costs under the improved RS scenario, removing RS regimen costs
out_all <- out_all %>% 
  mutate(cost_short_noRSdrugs=cost_short - unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_med_noRSdrugs=cost_med - unitcost_RS*initiate_all_RS*(1+wastage),
         cost_cea_noRSdrugs=cost_cea - unitcost_RS*initiate_all_RS_disc*(1+wastage),
         cost_short_hs_noRSdrugs=cost_short_hs - unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_med_hs_noRSdrugs=cost_med_hs - unitcost_RS*initiate_all_RS*(1+wastage),
         cost_cea_hs_noRSdrugs=cost_cea_hs - unitcost_RS*initiate_all_RS_disc*(1+wastage))
#next estimate improved RS drug costs when priced cost-neutrally, as above
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_RS2b=cost_short[scenario=="SQ"]-cost_short_noRSdrugs[scenario=="RS"],
         drugs_med_RS2b=cost_med[scenario=="SQ"]-cost_med_noRSdrugs[scenario=="RS"],
         drugs_cea_RS2b=cost_cea[scenario=="SQ"]-cost_cea_noRSdrugs[scenario=="RS"],
         drugs_short_hs_RS2b=cost_short_hs[scenario=="SQ"]-cost_short_hs_noRSdrugs[scenario=="RS"],
         drugs_med_hs_RS2b=cost_med_hs[scenario=="SQ"]-cost_med_hs_noRSdrugs[scenario=="RS"],
         drugs_cea_hs_RS2b=cost_cea_hs[scenario=="SQ"]-cost_cea_hs_noRSdrugs[scenario=="RS"])
#now estimate Pan drug costs when priced cost-neutrally vs. improved RS scenario w/ cost-neutral RS price
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_2b=cost_short_noRSdrugs[scenario=="RS"] + drugs_short_RS2b - cost_short[scenario==scenario_use],
         drugs_med_2b=cost_med_noRSdrugs[scenario=="RS"] + drugs_med_RS2b - cost_med[scenario==scenario_use],
         drugs_cea_2b=wtp*(dalys_cea[scenario=="RS"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_noRSdrugs[scenario=="RS"] + drugs_cea_RS2b - cost_cea[scenario==scenario_use],
         drugs_short_hs_2b=cost_short_hs_noRSdrugs[scenario=="RS"] + drugs_short_hs_RS2b - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_2b=cost_med_hs_noRSdrugs[scenario=="RS"] + drugs_med_hs_RS2b - cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_2b=wtp*(dalys_cea[scenario=="RS"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_hs_noRSdrugs[scenario=="RS"] + drugs_cea_hs_RS2b - cost_cea_hs[scenario==scenario_use])
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_2b=drugs_short_2b/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_2b=drugs_med_2b/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_2b=drugs_cea_2b/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_2b=drugs_short_hs_2b/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_2b=drugs_med_hs_2b/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_2b=drugs_cea_hs_2b/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))
#THIS IS THE SAME AS THE PAN THRESHOLD #1 - IF IMPROVED RS SCENARIO IS COST-NEUTRAL TO STATUS QUO, THEN IT SHOULD BE THE SAME.

# THRESHOLD 3A: PAN VS. IMPROVED RR (RR cost = BPaL[M] cost)
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_3a=cost_short[scenario=="RR"] - cost_short[scenario==scenario_use],
         drugs_med_3a=cost_med[scenario=="RR"]-cost_med[scenario==scenario_use],
         drugs_cea_3a=wtp*(dalys_cea[scenario=="RR"] - dalys_cea[scenario==scenario_use]) +
           cost_cea[scenario=="RR"] - cost_cea[scenario==scenario_use],
         drugs_short_hs_3a=cost_short_hs[scenario=="RR"] - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_3a=cost_med_hs[scenario=="RR"]-cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_3a=wtp*(dalys_cea[scenario=="RR"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_hs[scenario=="RR"] - cost_cea_hs[scenario==scenario_use])
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_3a=drugs_short_3a/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_3a=drugs_med_3a/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_3a=drugs_cea_3a/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_3a=drugs_short_hs_3a/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_3a=drugs_med_hs_3a/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_3a=drugs_cea_hs_3a/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))

# THRESHOLD 3B: PAN VS. IMPROVED RR (RR cost = cost-neutral vs. status quo)
#part 1: estimate cost-neutral costs for improved RR regimen vs. status quo
#first estimate costs under the improved RR scenario, removing RR regimen costs
out_all <- out_all %>% 
  mutate(cost_short_noRRdrugs=cost_short - unitcost_RR*initiate_primary_RR*(1+wastage),
         cost_med_noRRdrugs=cost_med - unitcost_RR*initiate_all_RR*(1+wastage),
         cost_cea_noRRdrugs=cost_cea - unitcost_RR*initiate_all_RR_disc*(1+wastage),
         cost_short_hs_noRRdrugs=cost_short_hs - unitcost_RR*initiate_primary_RR*(1+wastage),
         cost_med_hs_noRRdrugs=cost_med_hs - unitcost_RR*initiate_all_RR*(1+wastage),
         cost_cea_hs_noRRdrugs=cost_cea_hs - unitcost_RR*initiate_all_RR_disc*(1+wastage))
#next estimate improved RR drug costs when priced cost-neutrally, as above
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_RR3b=cost_short[scenario=="SQ"]-cost_short_noRRdrugs[scenario=="RR"],
         drugs_med_RR3b=cost_med[scenario=="SQ"]-cost_med_noRRdrugs[scenario=="RR"],
         drugs_cea_RR3b=cost_cea[scenario=="SQ"]-cost_cea_noRRdrugs[scenario=="RR"],
         drugs_short_hs_RR3b=cost_short_hs[scenario=="SQ"]-cost_short_hs_noRRdrugs[scenario=="RR"],
         drugs_med_hs_RR3b=cost_med_hs[scenario=="SQ"]-cost_med_hs_noRRdrugs[scenario=="RR"],
         drugs_cea_hs_RR3b=cost_cea_hs[scenario=="SQ"]-cost_cea_hs_noRRdrugs[scenario=="RR"])
#now estimate Pan drug costs when priced cost-neutrally vs. improved RR scenario w/ cost-neutral RR price
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_3b=cost_short_noRRdrugs[scenario=="RR"] + drugs_short_RR3b - cost_short[scenario==scenario_use],
         drugs_med_3b=cost_med_noRRdrugs[scenario=="RR"] + drugs_med_RR3b - cost_med[scenario==scenario_use],
         drugs_cea_3b=wtp*(dalys_cea[scenario=="RR"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_noRRdrugs[scenario=="RR"] + drugs_cea_RR3b - cost_cea[scenario==scenario_use],
         drugs_short_hs_3b=cost_short_hs_noRRdrugs[scenario=="RR"] + drugs_short_hs_RR3b - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_3b=cost_med_hs_noRRdrugs[scenario=="RR"] + drugs_med_hs_RR3b - cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_3b=wtp*(dalys_cea[scenario=="RR"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_hs_noRRdrugs[scenario=="RR"] + drugs_cea_hs_RR3b - cost_cea_hs[scenario==scenario_use])
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_3b=drugs_short_3b/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_3b=drugs_med_3b/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_3b=drugs_cea_3b/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_3b=drugs_short_hs_3b/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_3b=drugs_med_hs_3b/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_3b=drugs_cea_hs_3b/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))
#This is again the same as the pan threshold #1

# THRESHOLD 3C: PAN VS. IMPROVED RR (RR cost = Pan cost)
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_3c=(cost_short[scenario==scenario_use] - cost_short_noRRdrugs[scenario=="RR"])/
           ((1+wastage)*(initiate_primary_RR[scenario=="RR"] - initiate_primary_P[scenario==scenario_use])),
         threshold_med_soc_3c=(cost_med[scenario==scenario_use] - cost_med_noRRdrugs[scenario=="RR"])/
           ((1+wastage)*(initiate_all_RR[scenario=="RR"] - initiate_all_P[scenario==scenario_use])),
         threshold_cea_soc_3c=(cost_cea[scenario==scenario_use] - cost_cea_noRRdrugs[scenario=="RR"])/
           ((1+wastage)*(initiate_all_RR_disc[scenario=="RR"] - initiate_all_P_disc[scenario==scenario_use])),
         threshold_short_hs_3c=(cost_short_hs[scenario==scenario_use] - cost_short_hs_noRRdrugs[scenario=="RR"])/
           ((1+wastage)*(initiate_primary_RR[scenario=="RR"] - initiate_primary_P[scenario==scenario_use])),
         threshold_med_hs_3c=(cost_med_hs[scenario==scenario_use] - cost_med_hs_noRRdrugs[scenario=="RR"])/
           ((1+wastage)*(initiate_all_RR[scenario=="RR"] - initiate_all_P[scenario==scenario_use])),
         threshold_cea_hs_3c=(cost_cea_hs[scenario==scenario_use] - cost_cea_hs_noRRdrugs[scenario=="RR"])/
           ((1+wastage)*(initiate_all_RR_disc[scenario=="RR"] - initiate_all_P_disc[scenario==scenario_use])))
#if price is less than this threshold, improved RR scenario costs less
#if price is greater than this threshold, Pan scenario costs less

# THRESHOLD 4A: PAN VS. IMPROVED RS AND RR (RS price = 6HRZE, RR price = BPaLM)
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_4a=cost_short[scenario=="RS_RR"] - cost_short[scenario==scenario_use],
         drugs_med_4a=cost_med[scenario=="RS_RR"]-cost_med[scenario==scenario_use],
         drugs_cea_4a=wtp*(dalys_cea[scenario=="RS_RR"] - dalys_cea[scenario==scenario_use]) +
           cost_cea[scenario=="RS_RR"] - cost_cea[scenario==scenario_use],
         drugs_short_hs_4a=cost_short_hs[scenario=="RS_RR"] - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_4a=cost_med_hs[scenario=="RS_RR"]-cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_4a=wtp*(dalys_cea[scenario=="RS_RR"] - dalys_cea[scenario==scenario_use]) +
           cost_cea_hs[scenario=="RS_RR"] - cost_cea_hs[scenario==scenario_use])
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_4a=drugs_short_4a/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_4a=drugs_med_4a/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_4a=drugs_cea_4a/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_4a=drugs_short_hs_4a/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_4a=drugs_med_hs_4a/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_4a=drugs_cea_hs_4a/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))

# THRESHOLD 4B: PAN VS. IMPROVED RS AND RR (RS price = cost-neutral vs. SOC as in 2b; RR price = cost-neutral vs. SOC as in 3A)
#first estimate costs under the improved RS/RR scenario without either drugs
out_all <- out_all %>% 
  mutate(cost_short_noRSRRdrugs = cost_short - unitcost_RR*initiate_primary_RR*(1+wastage) -
           unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_med_noRSRRdrugs = cost_med - unitcost_RR*initiate_all_RR*(1+wastage) -
           unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_cea_noRSRRdrugs = cost_cea - unitcost_RR*initiate_all_RR_disc*(1+wastage) -
           unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_short_hs_noRSRRdrugs = cost_short_hs - unitcost_RR*initiate_primary_RR*(1+wastage) -
           unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_med_hs_noRSRRdrugs = cost_med_hs - unitcost_RR*initiate_all_RR*(1+wastage) -
           unitcost_RS*initiate_primary_RS*(1+wastage),
         cost_cea_hs_noRSRRdrugs = cost_cea_hs - unitcost_RR*initiate_all_RR_disc*(1+wastage) -
           unitcost_RS*initiate_primary_RS*(1+wastage))
#next estimate price at which improved RS is cost-neutral vs. SOC and improved RR is cost-neutral vs. SOC
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_RS=drugs_short_RS2b/((1+wastage)*initiate_primary_RS[scenario=="RS"]),
         threshold_med_soc_RS=drugs_med_RS2b/((1+wastage)*initiate_all_RS[scenario=="RS"]),
         threshold_cea_soc_RS=drugs_cea_RS2b/((1+wastage)*initiate_all_RS_disc[scenario=="RS"]),
         threshold_short_hs_RS=drugs_short_hs_RS2b/((1+wastage)*initiate_primary_RS[scenario=="RS"]),
         threshold_med_hs_RS=drugs_med_hs_RS2b/((1+wastage)*initiate_all_RS[scenario=="RS"]),
         threshold_cea_hs_RS=drugs_cea_hs_RS2b/((1+wastage)*initiate_all_RS_disc[scenario=="RS"]))
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_RR=drugs_short_RR3b/((1+wastage)*initiate_primary_RR[scenario=="RR"]),
         threshold_med_soc_RR=drugs_med_RR3b/((1+wastage)*initiate_all_RR[scenario=="RR"]),
         threshold_cea_soc_RR=drugs_cea_RR3b/((1+wastage)*initiate_all_RR_disc[scenario=="RR"]),
         threshold_short_hs_RR=drugs_short_hs_RR3b/((1+wastage)*initiate_primary_RR[scenario=="RR"]),
         threshold_med_hs_RR=drugs_med_hs_RR3b/((1+wastage)*initiate_all_RR[scenario=="RR"]),
         threshold_cea_hs_RR=drugs_cea_hs_RR3b/((1+wastage)*initiate_all_RR_disc[scenario=="RR"]))
#next estimate price at which Pan is cost-neutral vs. improved RR and RS
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_4b=cost_short_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_short_soc_RR*(1+wastage)*initiate_primary_RR[scenario=="RS_RR"] +
           threshold_short_soc_RS*(1+wastage)*initiate_primary_RS[scenario=="RS_RR"] - cost_short[scenario==scenario_use],
         drugs_med_4b=cost_med_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_med_soc_RR*(1+wastage)*initiate_all_RR[scenario=="RS_RR"] +
           threshold_med_soc_RS*(1+wastage)*initiate_all_RS[scenario=="RS_RR"] - cost_med[scenario==scenario_use],
         drugs_cea_4b=cost_cea_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_cea_soc_RR*(1+wastage)*initiate_all_RR_disc[scenario=="RS_RR"] +
           threshold_cea_soc_RS*(1+wastage)*initiate_all_RS_disc[scenario=="RS_RR"] - cost_cea[scenario==scenario_use],
         drugs_short_hs_4b=cost_short_hs_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_short_hs_RR*(1+wastage)*initiate_primary_RR[scenario=="RS_RR"] +
           threshold_short_hs_RS*(1+wastage)*initiate_primary_RS[scenario=="RS_RR"] - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_4b=cost_med_hs_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_med_hs_RR*(1+wastage)*initiate_all_RR[scenario=="RS_RR"] +
           threshold_med_hs_RS*(1+wastage)*initiate_all_RS[scenario=="RS_RR"] - cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_4b=cost_cea_hs_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_cea_hs_RR*(1+wastage)*initiate_all_RR_disc[scenario=="RS_RR"] +
           threshold_cea_hs_RS*(1+wastage)*initiate_all_RS_disc[scenario=="RS_RR"] - cost_cea_hs[scenario==scenario_use])
#finally estimate corresponding price
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_4b=drugs_short_4b/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_4b=drugs_med_4b/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_4b=drugs_cea_4b/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_4b=drugs_short_hs_4b/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_4b=drugs_med_hs_4b/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_4b=drugs_cea_hs_4b/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))

# THRESHOLD 4C: PAN VS. IMPROVED RS AND RR (RS price = 6HRZE; RR price = cost-neutral vs. SOC as in 3A)
#estimate price at which Pan is cost-neutral vs. improved RR and RS
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(drugs_short_4c=cost_short_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_short_soc_RR*(1+wastage)*initiate_primary_RR[scenario=="RS_RR"] +
           unitcost_RS*(1+wastage)*initiate_primary_RS[scenario=="RS_RR"] - cost_short[scenario==scenario_use],
         drugs_med_4c=cost_med_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_med_soc_RR*(1+wastage)*initiate_all_RR[scenario=="RS_RR"] +
           unitcost_RS*(1+wastage)*initiate_all_RS[scenario=="RS_RR"] - cost_med[scenario==scenario_use],
         drugs_cea_4c=cost_cea_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_cea_soc_RR*(1+wastage)*initiate_all_RR_disc[scenario=="RS_RR"] +
           unitcost_RS*(1+wastage)*initiate_all_RS_disc[scenario=="RS_RR"] - cost_cea[scenario==scenario_use],
         drugs_short_hs_4c=cost_short_hs_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_short_hs_RR*(1+wastage)*initiate_primary_RR[scenario=="RS_RR"] +
           unitcost_RS*(1+wastage)*initiate_primary_RS[scenario=="RS_RR"] - cost_short_hs[scenario==scenario_use],
         drugs_med_hs_4c=cost_med_hs_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_med_hs_RR*(1+wastage)*initiate_all_RR[scenario=="RS_RR"] +
           unitcost_RS*(1+wastage)*initiate_all_RS[scenario=="RS_RR"] - cost_med_hs[scenario==scenario_use],
         drugs_cea_hs_4c=cost_cea_hs_noRSRRdrugs[scenario=="RS_RR"] + 
           threshold_cea_hs_RR*(1+wastage)*initiate_all_RR_disc[scenario=="RS_RR"] +
           unitcost_RS*(1+wastage)*initiate_all_RS_disc[scenario=="RS_RR"] - cost_cea_hs[scenario==scenario_use])
#next estimate corresponding price
out_all <- out_all %>% group_by(country, sim) %>%
  mutate(threshold_short_soc_4c=drugs_short_4c/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_soc_4c=drugs_med_4c/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_soc_4c=drugs_cea_4c/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)),
         threshold_short_hs_4c=drugs_short_hs_4c/(initiate_primary_P[scenario==scenario_use]*(1+wastage)),
         threshold_med_hs_4c=drugs_med_hs_4c/(initiate_all_P[scenario==scenario_use]*(1+wastage)),
         threshold_cea_hs_4c=drugs_cea_hs_4c/(initiate_all_P_disc[scenario==scenario_use]*(1+wastage)))

# TAKES AVERAGES AND PERCENTILES ACROSS SIMS #
out_sum <- out_all %>% ungroup() %>% group_by(country, wtp) %>% select(starts_with("threshold")) %>%
  summarise_all(list("mean"=~mean(.), "lb"=~quantile(., p=0.025), "ub"=~quantile(., p=0.975)))

# PIVOT FOR TABLES #
out_sum_long <- out_sum %>% rename_with(~str_remove(., "threshold_"))
out_sum_long <- pivot_longer(out_sum_long, cols=names(out_sum_long)[3:ncol(out_sum_long)], names_sep="_",
                             names_to=c("type", "perspective", "comparison", "estimate"), values_to="threshold")
out_sum_long <- out_sum_long %>% group_by(country, wtp, type, perspective, comparison) %>%
  mutate(value_lab=paste0(round(threshold[estimate=="mean"], -1), " [", round(threshold[estimate=="lb"], -1), "-",
                          round(threshold[estimate=="ub"], -1), "]"))
out_sum_long <- out_sum_long %>% 
  mutate(perspective=if_else(perspective=="soc", "Societal", "Health Systems"))

# SAVE TO FILE #
write.csv(out_sum_long, file=paste0(path_out, "thresholds_", scenario_use, "_", analysis, ".csv"),
          row.names=F) 


