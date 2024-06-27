library(tidyverse)
library(flextable)
library(officer)

path_out <- "output/"

#LOAD OUTPUT FROM MAIN MODEL
countries <- c("india")
cures_reg_comb <- data.frame()
cures_reg_init_comb <- data.frame()
for(i in countries) {
  load(paste0(path_out, "cures_reg_", i, ".Rda"))
  load(paste0(path_out, "cures_reg_init_", i, ".Rda"))
  cures_reg_all <- bind_rows(cures_reg_all, .id="scenario") %>% mutate(country=i)
  cures_reg_all_init <- bind_rows(cures_reg_all_init, .id="scenario") %>% mutate(country=i)
  cures_reg_comb <- rbind(cures_reg_comb, cures_reg_all)
  cures_reg_init_comb <- rbind(cures_reg_init_comb, cures_reg_all_init)
}


scenarios <- names(cures_reg_all)
names(scenarios) <- c("Status Quo",
                      "Pan-TB TRP",
                      "Improved RS Regimen",
                      "Improved RR Regimen",
                      "Improved RR Retention",
                      "Improved Rif DST",
                      "Improved RS & RR Regimens",
                      "Pan-TB TRP: Long-Acting Injectable")

#cures by resistance and regimen assignment among all diagnoses
cures_reg_init_comb  <- cures_reg_init_comb %>% mutate(scenario_lab=names(scenarios)[match(scenario, scenarios)])
cures_reg_init_comb <- cures_reg_init_comb %>% mutate(scenario_lab=factor(scenario_lab, levels=names(scenarios)))
cures_reg_init_sum <- cures_reg_init_comb %>% group_by(scenario_lab, scenario) %>%
  select(-c(sim, year, country)) %>%
  summarise_all(list(mean=~mean(., na.rm=T),
                     lb=~quantile(., p=0.025, na.rm=T),
                     ub=~quantile(., p=0.975, na.rm=T)))
cures_reg_init_sub <- cures_reg_init_sum %>% filter(scenario %in% c("SQ", "Pan")) %>%
  ungroup() %>% select(-c(scenario, scenario_lab)) %>%
  mutate_all(~sum(., na.rm=T))
cures_reg_init_sub <- unique(cures_reg_init_sub)
means <- cures_reg_init_sub %>% select(ends_with("mean"))
lbs <- cures_reg_init_sub %>% select(ends_with("lb"))
ubs <- cures_reg_init_sub %>% select(ends_with("ub"))


table_init <- cbind(paste0(format(round(means %>% select(starts_with("none")), 3)*100, nsmall=1),
                      "% [",
                      format(round(lbs %>% select(starts_with("none")), 3)*100, nsmall=1),
                      "-",
                      format(round(ubs %>% select(starts_with("none")), 3)*100, nsmall=1),
                      "%]"),
               paste0(format(round(means %>% select(starts_with("R_")), 3)*100, nsmall=1),
                      "% [",
                      format(round(lbs %>% select(starts_with("R_")), 3)*100, nsmall=1),
                      "-",
                      format(round(ubs %>% select(starts_with("R_")), 3)*100, nsmall=1),
                      "%]"
               ),
               paste0(format(round(means %>% select(starts_with("N_")), 3)*100, nsmall=1),
                      "% [",
                      format(round(lbs %>% select(starts_with("N_")), 3)*100, nsmall=1),
                      "-",
                      format(round(ubs %>% select(starts_with("N_")), 3)*100, nsmall=1),
                      "%]"
               ),
               paste0(format(round(means %>% select(starts_with("RN")), 3)*100, nsmall=1),
                      "% [",
                      format(round(lbs %>% select(starts_with("RN")), 3)*100, nsmall=1),
                      "-",
                      format(round(ubs %>% select(starts_with("RN")), 3)*100, nsmall=1),
                      "%]"
               )
)
colnames(table_init) <- c("none", "R", "N", "RN")
rownames(table_init) <- c("RS SOC", "RR SOC", "Pan-TB TRP", "Individualized")
table_init[table_init=="0.0% [0.0-0.0%]"] <- "NA"
ft_init <- flextable(data = as.data.frame(table_init) %>% rownames_to_column("Regimen")) %>% 
  theme_box %>% 
  fontsize(size=8, part="all") %>%
  autofit

# Create a docx file with the table
read_docx() %>% 
  body_add_par("Initial treatment success, all diagnoses") %>%
  body_add_flextable(ft_init) %>% 
  body_add_par("Columns indicate resistance phenotypes (R = rif resistance; N = novel/BDQ resistance), rows indicate regimens. “NA” means nobody of a given resistance phenotype was assigned to a given regimen in our model.",
               style="Normal") %>%
  print(target = paste0(path_out, "cures_table_pooled", ".docx")) 
