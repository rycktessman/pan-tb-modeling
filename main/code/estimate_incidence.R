library(dplyr)
library(tidyr)
library(dampack)

path_out <- "output/"

analyses <- c("main", "resistN", "dstN", "worse_adhereLTFU", "better_adhereLTFU", "worse_Pan", "short_RS",
              "bad_efficacy", "good_efficacy", "bad_duration", "good_duration", "bad_adherence", "good_adherence",
              "bad_forgiveness", "good_forgiveness", "bad_safety", "good_safety",
              "bad_baseres", "good_baseres", "bad_restrend", "good_restrend")  
analysis <- analyses[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
print(analysis)

#basic info
countries <- c("india", "southafrica", "philippines")
names(countries) <- c("India", "South Africa", "Philippines")
#country <- countries[as.numeric(Sys.getenv('country_num'))]
#print(country)

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


#load output for all countries into 1 dataframe
out_comb <- data.frame()
for(i in names(countries)) {
  load(paste0("output/out_", countries[i], "_", analysis, ".Rda"))
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
out_comb <- out_comb %>% select(country, scenario, scenario_lab, sim, year, failures) %>%
  arrange(country, scenario_lab, sim, year)
n.sim <- max(out_comb$sim)

#load estimated incidence and case detection ratios for each country - use pre-covid (2019 for CDR, 2020 for inc)
inc <- read.csv("data/TB_burden_countries_2023-08-29.csv")
inc <- inc %>% filter(country %in% names(countries))
inc <- inc %>% select(country, year, starts_with("e_inc_100k"), starts_with("c_cdr"))
inc <- inc %>% group_by(country) %>% 
  summarise(inc=e_inc_100k[year==2020], 
            inc_lb=e_inc_100k_lo[year==2020],
            inc_ub=e_inc_100k_hi[year==2020],
            cdr=c_cdr[year==2019],
            cdr_lb=c_cdr_lo[year==2019],
            cdr_ub=c_cdr_hi[year==2019])
inc <- inc %>% mutate(cdr_ub=pmin(95, cdr_ub)) #> 100 reflects non-TB cases diagnosed - not relevant to our model, so cap at 95

#sample from incidence distributions
inc_samples <- sapply(names(countries), function(x)
  rnorm(n=n.sim, mean=inc %>% filter(country==x) %>% pull(inc),
        sd=((inc %>% filter(country==x) %>% pull(inc_ub))-
              (inc %>% filter(country==x) %>% pull(inc_lb)))/(2*1.96)),
  simplify=F, USE.NAMES=T)
inc_samples <- bind_rows(inc_samples) %>% mutate(sim=1:n())
inc_samples <- pivot_longer(inc_samples, cols=names(countries), names_to="country", values_to="inc")

#load CDR samples
load("data/params.Rda")
cdr_samples <- lapply(countries, function(x)
  data.frame(sim=1:length(params[[x]][["cdr"]]),
             cdr=params[[x]][["cdr"]]))
cdr_samples <- bind_rows(cdr_samples, .id="country")

#merge with modeled output
out_comb <- left_join(out_comb, inc_samples, by=c("country", "sim"))
out_comb <- left_join(out_comb, cdr_samples, by=c("country", "sim"))

#load serial interval distribution
serial <- read.csv("data/serial_interval.csv") %>% filter(year<11 & year>0) %>%
  pull(p_tb_inc)

#expand output into half year increments
out_comb2 <- out_comb %>% mutate(year=year+0.5) %>%
  group_by(country, scenario, sim) %>%
  mutate(failures=if_else(year<10,
                          (failures+lead(failures))/2,
                          failures)) 
out_comb <- rbind(out_comb, out_comb2) %>%
  arrange(country, scenario, sim, year)

out_comb <- out_comb %>% group_by(country, sim, year) %>%
  mutate(fail_SQ=failures[scenario=="SQ"])
out_comb <- out_comb %>% mutate(cohort=50, #biannual cohorts of 50
                                cases=0,
                                secondary=0)

#apply serial interval distribution to estimate % and size of cohort that is from pre-year 1 cases
#add extra 0.5 months - to account for delay to infectiousness after tx failure
out_comb <- left_join(out_comb, 
                      data.frame("year"=seq(from=1, to=12, by=0.5), 
                                 "serial"=c(0, 0, serial)), by="year")
out_comb <- out_comb %>% group_by(country, scenario, sim) %>%
  mutate(serial_sum=cumsum(serial))
out_comb <- out_comb %>% mutate(cohort0=cohort*(1-serial_sum))
out_comb <- out_comb %>% select(-c(serial, serial_sum))

serial <- c(0, serial)
names(serial) <- as.character(seq(from=0.5, to=11, by=0.5))

for(yr in seq(from=1, to=10, by=0.5)) {
  print(yr)
  for(yrs_out in seq(from=0.5, to=10.5-yr, by=0.5)) {
    out_comb <- out_comb %>% ungroup() %>% group_by(country, sim, scenario) %>%
      mutate(cases=if_else(year==yr+yrs_out, 
                           cases + cohort[year==yr]*serial[as.character(yrs_out)]*
                             (1+cdr*failures[year==yr])/(1+cdr*fail_SQ[year==yr]), 
                           cases),
             #for estimation of costs and deaths, also track just 2ndary cases resulting from failures each yr
             #only a % of these secondary cases would actually get treated (need to apply CDR again)
             secondary=if_else(year==yr+yrs_out,
                               secondary + cohort[year==yr]*serial[as.character(yrs_out)]*
                                 cdr*failures[year==yr]/(1+cdr*fail_SQ[year==yr]),
                               secondary))
  }
  out_comb <- out_comb %>% mutate(cohort=if_else(year==(yr+0.5), cohort0+cases, cohort))
}
out_comb <- out_comb %>% ungroup() %>% group_by(country, year, sim) %>% 
  mutate(cases_averted=cohort-cohort[scenario=="SQ"],
         secondary_averted=secondary-secondary[scenario=="SQ"])
#cases_averted implicitly already includes tertiary, 4th generation, etc. cases - bc they are included in subsequent cohorts

#now apply to incidence
out_comb <- out_comb %>% group_by(country, sim) %>%
  mutate(cases_change=if_else(year==1, 1, (50+cases_averted)/50))
out_comb <- out_comb %>% ungroup() %>% mutate(inc=inc*cases_change)
out_comb <- out_comb %>% group_by(country, sim) %>%
  mutate(inc_change=(inc-inc[scenario=="SQ"])/inc[scenario=="SQ"])

#plot and save incidence output
#take means and CIs
out_sum <- out_comb %>% group_by(country, year, scenario, scenario_lab) %>%
  summarise(inc_mean=mean(inc), inc_lb=quantile(inc, 0.025), inc_ub=quantile(inc, 0.975),
            inc_change_mean=mean(inc_change), inc_change_lb=quantile(inc_change, 0.025), 
            inc_change_ub=quantile(inc_change, 0.975),
            cases_change=mean(cases_change), cases_change_lb=quantile(cases_change, 0.025),
            cases_change_ub=quantile(cases_change, 0.975),
            cohort=mean(cohort),
            cases=mean(cases),
            secondary=mean(secondary)
  )

#cumulative incidence reduction
out_sum %>% filter(year==10.5) %>% select(scenario_lab, starts_with("inc_change"))

write.csv(out_comb, file=paste0(path_out, "incidence_out_all", "_", analysis, ".csv"), row.names=F)
write.csv(out_sum, file=paste0(path_out, "incidence_out_sum", "_", analysis, ".csv"), row.names=F)

#also save secondary cases and cohort sizes to use in final cost, death, DALY calculations
#change back to annual timestep for consistency with other output
out_yr <- out_comb %>% 
  select(country, scenario, scenario_lab, sim, year, cohort, secondary) %>%
  mutate(year_int=floor(year)) %>%
  group_by(country, scenario, scenario_lab, sim, year_int) %>%
  select(-year) %>%
  summarise(cohort=sum(cohort), secondary=sum(secondary))
write.csv(out_yr, file=paste0(path_out, "secondary_out_all", "_", analysis, ".csv"), row.names=F)


