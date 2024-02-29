#load packages
library(abind)
library(dplyr)
library(tictoc) 

#load parameters
load("data/params.Rda")
load("data/regimen_params.Rda")
load("data/cost_params.Rda")
load("data/cost_params_LAI.Rda")

#load functions
source("code/model_functions.R")

#set up model
#country <- "india"
analyses <- c("main", "resistN", "dstN", "worse_adhereLTFU", "better_adhereLTFU", "worse_Pan", "short_RS")  
country <- Sys.getenv('country')
analyses <- analyses[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(analyses)

out_all <- list()
cures_reg_all <- list()
cures_reg_all_init <- list()

#maximum duration across regimens - used as an input throughout
max_dur <- max(c(reg_params[[1]][[1]]$dur_RS, reg_params[[1]][[1]]$dur_RR,
                 reg_params[[1]][[1]]$dur_P, reg_params[[1]][[1]]$dur_I))

#specify which costs to exclude in HS perspective
costs_exclude_hs <- c("oop", "indirect")

#specify indices of each regimen in the matrix parameters
ind_RS <- seq(from=1, to=16, by=4)
ind_RR <- seq(from=2, to=16, by=4)
ind_P <- seq(from=3, to=16, by=4)
ind_I <- seq(from=4, to=16, by=4)

tic()

for(analysis in analyses) {
  print(analysis)
  for(scenario in names(reg_params[[country]])) {
    print(scenario)
    if(scenario=="Pan_inject") {
      pars <- c(params[[country]], reg_params[[country]][[scenario]], cost_params_lai[[country]])
    } else {
      pars <- c(params[[country]], reg_params[[country]][[scenario]], cost_params[[country]])
    }
    #adjust cost schedules for shorter RS and RR regimen scenarios (use Pan schedule)
    if(scenario %in% c("RS", "RS_RR")) {
      pars[["outpatient_sched"]][ind_RS, ] <- pars[["outpatient_sched"]][ind_P, ]
      pars[["smear_sched"]][ind_RS, ] <- pars[["smear_sched"]][ind_P, ]
      pars[["culture_sched"]][ind_RS, ] <- pars[["culture_sched"]][ind_P, ]
      pars[["cxr_sched"]][ind_RS, ] <- pars[["cxr_sched"]][ind_P, ]
      pars[["livertest_sched"]][ind_RS, ] <- pars[["livertest_sched"]][ind_P, ]
      pars[["fullblood_sched"]][ind_RS, ] <- pars[["fullblood_sched"]][ind_P, ]
      pars[["ecg_sched"]][ind_RS, ] <- pars[["ecg_sched"]][ind_P, ]
      pars[["neuroscreen_sched"]][ind_RS, ] <- pars[["neuroscreen_sched"]][ind_P, ]
      
      pars[["support_sched"]]  <- lapply(1:length(pars$support_sched), function(x)
        pars$support_sched[[x]][c(3,2,3,4,7,6,7,8,11,10,11,12,15,14,15,16),])
      pars[["ooptx_sched"]]  <- lapply(1:length(pars$ooptx_sched), function(x)
        pars$ooptx_sched[[x]][c(3,2,3,4,7,6,7,8,11,10,11,12,15,14,15,16),])
      pars[["indirecttx_sched"]]  <- lapply(1:length(pars$indirecttx_sched), function(x)
        pars$indirecttx_sched[[x]][c(3,2,3,4,7,6,7,8,11,10,11,12,15,14,15,16),])
      
      pars[["p_liverAE"]][,ind_RS] <- pars[["p_liverAE"]][,ind_P]
      pars[["p_pancreasAE"]][,ind_RS] <- pars[["p_pancreasAE"]][,ind_P]
      pars[["p_anemiaAE"]][,ind_RS] <- pars[["p_anemiaAE"]][,ind_P]
      pars[["p_neutropeniaAE"]][,ind_RS] <- pars[["p_neutropeniaAE"]][,ind_P]
      pars[["p_qtcfAE"]][,ind_RS] <- pars[["p_qtcfAE"]][,ind_P]
      pars[["p_renalAE"]][,ind_RS] <- pars[["p_renalAE"]][,ind_P]
      pars[["p_visionAE"]][,ind_RS] <- pars[["p_visionAE"]][,ind_P]
      pars[["p_arthralgiaAE"]][,ind_RS] <- pars[["p_arthralgiaAE"]][,ind_P]
      pars[["p_neuroshortAE"]][,ind_RS] <- pars[["p_neuroshortAE"]][,ind_P]
      pars[["p_neurolongAE"]][,ind_RS] <- pars[["p_neurolongAE"]][,ind_P]
      
      pars[["p_hosp_RS"]] <- pars[["p_hosp_P"]]
      
    }
    if(scenario %in% c("RR", "RS_RR")) {
      pars[["outpatient_sched"]][ind_RR, ] <- pars[["outpatient_sched"]][ind_P, ]
      pars[["smear_sched"]][ind_RR, ] <- pars[["smear_sched"]][ind_P, ]
      pars[["culture_sched"]][ind_RR, ] <- pars[["culture_sched"]][ind_P, ]
      pars[["cxr_sched"]][ind_RR, ] <- pars[["cxr_sched"]][ind_P, ]
      pars[["livertest_sched"]][ind_RR, ] <- pars[["livertest_sched"]][ind_P, ]
      pars[["fullblood_sched"]][ind_RR, ] <- pars[["fullblood_sched"]][ind_P, ]
      pars[["ecg_sched"]][ind_RR, ] <- pars[["ecg_sched"]][ind_P, ]
      pars[["neuroscreen_sched"]][ind_RR, ] <- pars[["neuroscreen_sched"]][ind_P, ]
      
      pars[["support_sched"]]  <- lapply(1:length(pars$support_sched), function(x)
        pars$support_sched[[x]][c(1,3,3,4,5,7,7,8,9,11,11,12,13,15,15,16),])
      pars[["ooptx_sched"]]  <- lapply(1:length(pars$ooptx_sched), function(x)
        pars$ooptx_sched[[x]][c(1,3,3,4,5,7,7,8,9,11,11,12,13,15,15,16),])
      pars[["indirecttx_sched"]]  <- lapply(1:length(pars$indirecttx_sched), function(x)
        pars$indirecttx_sched[[x]][c(1,3,3,4,5,7,7,8,9,11,11,12,13,15,15,16),])
      
      pars[["p_liverAE"]][,ind_RR] <- pars[["p_liverAE"]][,ind_P]
      pars[["p_pancreasAE"]][,ind_RR] <- pars[["p_pancreasAE"]][,ind_P]
      pars[["p_anemiaAE"]][,ind_RR] <- pars[["p_anemiaAE"]][,ind_P]
      pars[["p_neutropeniaAE"]][,ind_RR] <- pars[["p_neutropeniaAE"]][,ind_P]
      pars[["p_qtcfAE"]][,ind_RR] <- pars[["p_qtcfAE"]][,ind_P]
      pars[["p_renalAE"]][,ind_RR] <- pars[["p_renalAE"]][,ind_P]
      pars[["p_visionAE"]][,ind_RR] <- pars[["p_visionAE"]][,ind_P]
      pars[["p_arthralgiaAE"]][,ind_RR] <- pars[["p_arthralgiaAE"]][,ind_P]
      pars[["p_neuroshortAE"]][,ind_RR] <- pars[["p_neuroshortAE"]][,ind_P]
      pars[["p_neurolongAE"]][,ind_RR] <- pars[["p_neurolongAE"]][,ind_P]
      
      pars[["p_hosp_RR"]] <- pars[["p_hosp_P"]]
      
    }
    
    #adjust parameters for scenario analyses
    if(analysis=="resistN") { #faster trend in novel resistance emergence & higher baseline prevalence
      pars$resist_trend[,c(2,3)] <- 2*pars$resist_trend[,c(2,3)]
      pars$pop_start <- pars$pop_start_alt 
    } else if(analysis=="dstN" & scenario %in% c("Pan", "Pan_inject")) { #introduce DST for the Pan regimen in yr 6
      pars$dstN_RS <- reg_params[[country]][["SQ"]][["dstN_RR"]] #only adjust RS because no upfront RR testing in this scenario
      pars$dstR <- matrix(1, nrow=10000, ncol=10) #everyone who gets novel DST gets rif DST 
      pars$dstR_retx <- rep(1, 10000) #everyone who gets novel DST gets rif DST 
      pars$dstR_after_N <- TRUE #include rif DST after novel to assign to either RS SOC or individualized
    } else if(analysis=="worse_adhereLTFU") { #worse adherence, LTFU, and pre-tx LTFU under all scenarios (to the extent possible - RR and Ind stay with poor adherence)
      if(scenario!="Pan_inject") {
        pars$adherence_P <- pars$adherence_P + rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_P))
      }
      pars$adherence_RS <- pars$adherence_RS + rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_RS))
      pars$adherence_RR <- pars$adherence_RR + rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_RR))
      pars$adherence_I <- pars$adherence_I + rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_I))
      pars$ltfu <- pars$ltfu*1.5
      pars$initialloss_base <- pars$initialloss_base*1.5
      pars$initialloss_extra <- pars$initialloss_extra*1.5
    } else if(analysis=="better_adhereLTFU") { #better adherence, LTFU, and pre-tx LTFU under all scenarios
      if(scenario!="Pan_inject") {
        pars$adherence_P <- pars$adherence_P - rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_P))
      }
      pars$adherence_RS <- pars$adherence_RS - rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_RS))
      pars$adherence_RR <- pars$adherence_RR - rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_RR))
      pars$adherence_I <- pars$adherence_I - rep(c(0, 0.1, 0, 0, -0.1), each = nrow(pars$adherence_I))
      if(sum(pars$adherence_P<0)>1) {
        flags <- which(pars$adherence_P<0, arr.ind=TRUE)[,"row"]
        pars$adherence_P[flags,3] <- pars$adherence_P[flags,3] + pars$adherence_P[flags,2]
        pars$adherence_P[flags,2] <- 0
      }
      pars$ltfu <- pars$ltfu*0.5
      pars$initialloss_base <- pars$initialloss_base*0.5
      pars$initialloss_extra <- pars$initialloss_extra*0.5
    } else if(analysis=="worse_Pan") { #Pan regimen is less of an improvement over the RS SOC
      if(scenario!="Pan_inject") {
        pars$adherence_P <- (reg_params[[country]][["SQ"]][["adherence_RS"]] + 
                             reg_params[[country]][["SQ"]][["adherence_P"]])/2
        pars$forgive_P <- rep(0.1, length(pars$forgive_P))
      }
      if(scenario %in% c("RS", "RS_RR")) {
        pars$adherence_RS <- pars$adherence_P
        pars$forgive_RS <- pars$forgive_P
      }
      if(scenario %in% c("RR", "RS_RR")) {
        pars$adherence_RR <- pars$adherence_P
        pars$forgive_RR <- pars$forgive_P
      }
    } else if(analysis=="short_RS" & !(scenario %in% c("RS", "RS_RR"))) {
      pars$dur_RS <- 17
    }
    
    out_all_scen <- data.frame()
    cures_reg_all_scen <- data.frame() #track initial cures over time by regimen assignment
    cures_reg_all_init_scen <- data.frame() #track initial cures over time by regimen initiation
    
    for(j in 1:10000) {
      print(j)
      #cat('\r', paste(round(j/10000 * 100), "% done", sep = " ")) # display the progress of the simulation
      pop_all <- pars$pop_start[j,] #track pop by resistance over time
      deaths_all <- c() #track deaths over time
      cures_all <- c() #track cures over time - this is % initiating in a given yr that eventually get cured (either initially or upon retx)
      cures_init_all <- c() #track initial cures over time too
      fail_all <- c() #track total failures (from relapse, failure, initial loss) left alive over time
      costs_all <- c() #track costs by category over time
      dalys_all <- c() #track DALYs by category over time
      startmat_all <- c() #track % initiating each reg (including retx & secondary) over time
      
      pop <- pop_all #set starting pop
      retx_tx <- list("tx_reg"=0, "tx"=0) #set to 0 for 1st year since no retx yet
      
      for(i in 1:10) { #loop over 10 years time horizon
        #1. PRETREATMENT LTFU, DST, REGIMEN INITIATION, COSTS
        pretx <- calc_pretx(pop, pars, i, j, type="new")
        startmat <- pretx[[1]] #regimen initiation - used in the model
        assignmat <- pretx[[2]] #regimen assignment - only used to calculate "effectiveness" values of the regimens (for Aim 3)
        costs_pretx <- pretx[[3]] #depends on DST etc. so easier to wrap into calc_pretx function
        
        #calculate tx %s (for resistance trends) - weighted avg of prop new this year plus prop retx the previous year
        new_tx <- calc_num_tx(startmat, pars)
        tx_props <- (new_tx$tx_reg + retx_tx$tx_reg)/(new_tx$tx + retx_tx$tx) #% treated w/ R and N 
        
        #outcomes for pre-treatment LTFU
        initialloss <- pop - colSums(startmat)
        
        #2. OUTCOMES & COSTS WHILE ON TX
        cures <- treatment_model(startmat=startmat, pars=pars, index=j)
        
        #calculate costs while on treatment
        costs_dalys_tx <- calc_costs_dalys_tx(startmat=startmat, pars=pars, index=j, max_dur=max_dur)
        costs_tx <- costs_dalys_tx[[1]]
        dalys_ae <- costs_dalys_tx[[2]]
        
        #calculate costs per person by resistance type
        costs_tx_sum <- rbind(colSums(costs_tx[1:4, ]),
                              colSums(costs_tx[5:8, ]),
                              colSums(costs_tx[9:12, ]),
                              colSums(costs_tx[13:16, ]))
        costs <- costs_pretx + costs_tx_sum
        costs_pp <- rowSums(costs)/pop
        costs_pp_hs <- rowSums(costs[, !names(costs) %in% costs_exclude_hs])/pop
        #these costs will accrue to retreatments, with arbitrary 10% inflation (for now)
        
        #other outcomes while on treatment, deaths/retreatments among those with pre-tx LTFU
        deaths_tx <- colSums(startmat - cures)*(pars$cfr_tx[j])
        #adjust CFR based on % of deaths that are on tx, so we don't double count them. use this across scenarios
        if(i==1 & scenario=="SQ") {
          params[[country]][["cfr_notx"]][j] <- pars$cfr[j]*(1-pars$cdr[j]*sum(startmat)*sum(deaths_tx))
          pars$cfr_notx[j] <- pars$cfr[j]*(1-pars$cdr[j]*sum(startmat)*sum(deaths_tx))
        }
        deaths_primary <- deaths_tx + initialloss*pars$cfr_notx[j] #apply CFR to pre-tx LTFU
        retx_initialloss <- initialloss*(1-pars$cfr_notx[j]) #those that remain alive get treated later (1 year later for now)
        failures <- colSums(startmat - cures)*(1-pars$cfr_tx[j])*(1-pars$prop_relapse[j]) #rest are non-cures
        relapses <- colSums(startmat - cures)*(1-pars$cfr_tx[j])*pars$prop_relapse[j] #rest are non-cures
        fail_all <- c(fail_all, sum(failures + relapses + retx_initialloss)) #total failures/relapses/non-initiators alive (and generating 2ndary cases)
        
        #3. POST-TREATMENT OUTCOMES
        #first, retreatments
        deaths_primary <- deaths_primary + relapses*pars$cfr_notx[j] #apply CFR before retreatment for relapses only
        retx_fail <- failures #for now, assume no mortality and all still linked to care
        retx_relapse <- relapses*(1-pars$cfr_notx[j])
        costs_retx <- (retx_fail + retx_relapse)*costs_pp*1.1 + retx_initialloss*costs_pp #don't inflate costs for those w/ pre-tx LTFU
        costs_retx_hs <- (retx_fail+retx_relapse)*costs_pp_hs*1.1 + retx_initialloss*costs_pp_hs 
        
        #calculate tx assignment of retreatments - just used for resistance trends the following year 
        retx_initialloss_startmat <- calc_pretx(retx_initialloss, pars, i, j, type="initialloss")[[1]]
        retx_fail_startmat <- calc_pretx(retx_fail, pars, i, j, type="failure")[[1]]
        retx_relapse_startmat <- calc_pretx(retx_relapse, pars, i, j, type="relapse")[[1]]
        retx_startmat <- retx_initialloss_startmat + retx_fail_startmat + retx_relapse_startmat
        retx_tx <- calc_num_tx(retx_startmat, pars)
        
        #calculate outcomes of retreatments
        retx <- retx_initialloss + retx_fail + retx_relapse
        cures_retx <- retx_initialloss*pars$RR_cure_initialloss[j]*(colSums(cures)/pop) + #same outcomes as primary patients
          retx_fail*pars$RR_cure_retx[j]*(colSums(cures)/pop) + #worse outcomes than primary patients
          retx_relapse*pars$RR_cure_retx[j]*(colSums(cures)/pop)  #worse outcomes than primary patients
        cures_retx[pop==0] <- 0 #to address 0/0 error above if scenario with 0 resistance type
        deaths_retx <- (retx - cures_retx)*pars$cfr_tx[j] #add to deaths the following year
        non_cures_retx <- (retx-cures_retx)*(1-pars$cfr_tx[j]) #retx that weren't cured but remain alive
        
        #add on more costs, dalys, and deaths for those that still aren't cured after 2nd tx attempt
        deaths_retx <- deaths_retx + non_cures_retx*pars$cfr_notx[j]
        non_cures_retx <- non_cures_retx*(1-pars$cfr_notx[j])
        costs_retx <- costs_retx + non_cures_retx*costs_pp*1.1
        retx <- retx + non_cures_retx 
        
        cures_reg_all_scen <- rbind(cures_reg_all_scen, c(j, i, c(cures/startmat))) # % cured on 1st tx among those starting a regimen
        cures_reg_all_init_scen <- rbind(cures_reg_all_init_scen, c(j, i, c(cures/assignmat)))  # % cured on 1st tx among those assigned to a regimen
        cures_init <- colSums(cures)
        cures <- colSums(cures) + cures_retx #add initial & retreated cures, remove tracking by reg assignment
        cures_init_all <- rbind(cures_init_all, cures_init)
        cures_all <- rbind(cures_all, cures)
        
        #sum startmat across resistance types (for costing regimen drugs), but keep primary txs only for short-term estimates
        startmat_sum <- startmat + retx_startmat 
        startmat_sum <- rowSums(startmat_sum)
        startmat_all <- rbind(startmat_all, c(rowSums(startmat), startmat_sum))
        
        #add to costs and sum across regimens
        costs <- cbind(costs, 
                       "retx"=costs_retx,
                       "retx_hs"=costs_retx_hs) #assuming costs accrue in same year - only matters for discounting
        costs <- colSums(costs)
        costs_all <- rbind(costs_all, costs)
        
        #calculate DALYs by resistance type
        dalys <- cbind("dalys_ae"=dalys_ae,
                       "dalys_tb"=pars$dalys_tb[j]*0.5*pop, 
                       "dalys_posttb"=pars$dalys_posttb[j]*pop,
                       "dalys_deaths"=deaths_primary*pars$disc_life_exp[j]
        )
        dalys <- cbind(dalys, 
                       "dalys_retx"=retx*pars$dalys_tb[j]*0.5 +
                         deaths_retx*pars$disc_life_exp[j] +
                         retx*dalys[,"dalys_ae"]) #no post TB because they already accrued it
        #sum across resistance types and save
        dalys <- colSums(dalys)
        dalys_all <- rbind(dalys_all, dalys)
        
        deaths <- c("deaths_primary"=sum(deaths_primary), 
                    "deaths_retx"=sum(deaths_retx))
        deaths_all <- rbind(deaths_all, deaths)
        
        #4. CALCULATE PREVALENCE OF RESISTANCE FOR THE FOLLOWING YEAR
        pop <- calc_prev_resist(pop, tx_props, pars$resist_trend[j,])
        pop_all <- rbind(pop_all, pop)
        pop <- pop_all[i+1,]
      }
      
      out <- data.frame(j,
                        1:10,
                        pop_all[1:10,],
                        cures_all,
                        cures_init_all,
                        fail_all,
                        deaths_all,
                        costs_all,
                        dalys_all,
                        startmat_all,
                        cfr_notx=pars$cfr_notx[j]) 
      names(out) <- c("sim", 
                      "year",
                      paste0("pop_", colnames(pop_all)),
                      paste0("cures_", colnames(cures_all)),
                      paste0("cures_init_", colnames(cures_all)),
                      "failures",
                      colnames(deaths_all),
                      paste0("cost_", colnames(costs_all)),
                      colnames(dalys_all),
                      paste0("initiate_primary_", colnames(startmat_all)[1:4]),
                      paste0("initiate_all_", colnames(startmat_all)[1:4]),
                      "cfr_notx")
      
      
      out_all_scen <- bind_rows(out_all_scen, out)
    }
    names(cures_reg_all_scen)  <- c("sim", "year",
                                    paste(rep(colnames(startmat), 
                                              each = dim(startmat)[1]), 
                                          rownames(startmat), sep = "_")) 
    names(cures_reg_all_init_scen) <- names(cures_reg_all_scen)
    
    rownames(out_all_scen) <- NULL
    out_all[[scenario]] <- out_all_scen
    cures_reg_all[[scenario]] <- cures_reg_all_scen
    cures_reg_all_init[[scenario]] <- cures_reg_all_init_scen
  }
  
  save(out_all, file=paste0("output/out_", country, "_", analysis, ".Rda"))
  
  if(analysis=="main") {
    save(cures_reg_all, file=paste0("output/cures_reg_", country, ".Rda"))
    save(cures_reg_all_init, file=paste0("output/cures_reg_init_", country, ".Rda"))
  }
}

toc()