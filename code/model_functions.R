#calculate prev. of resistance among new cases in the subsequent yr, based on tx among cases this year
calc_prev_resist <- function(pop, prop_tx, trend) {
  #calc % N resistant by RS vs. RR in 
  resist <- c(pop[2] + pop[4], #% RR (regardless of novel resistance status)
              pop[3]/(pop[1] + pop[3]), #% novel resistant among RS
              pop[4]/(pop[2] + pop[4]) #% novel resistant among RR)
  )
  resist_new <- resist + trend*prop_tx
  #change back to none, R only, N only, RN
  resist_new <- c(resist_new[1]*(1-resist_new[3]), # % N susceptible, R resistant
                  resist_new[2]*(1-resist_new[1]), # % N resistant, R susceptible
                  resist_new[3]*resist_new[1] # % N resistant and R resistant
                  )
  resist_new <- c(1-sum(resist_new), resist_new)
  return(resist_new)
}

#calculate pre-treatment outcomes (DST, assignment, LTFU) - varies by primary vs. retx
calc_pretx <- function(pop, pars, time, index=1, type="new") {
  #initial diagnosis - assume everyone has Xpert, 1 outpatient visit, and predx indirect costs
  outpatient_costs <- pop*pars$outpatient_c[index]
  lab_costs <- pop*pars$xpert_c[index]
  oop_costs <- pop*pars$oopdx_base_c[index]
  indirect_costs <- pop*pars$indirectdx_base_c[index]
  
  if(!(pars$dstR_after_N)) {
    #RR DST step first
    pop_dstR <- pop*pars$dstR[index, time]*(type %in% c("new", "initialloss")) + 
      pop*pars$dstR_retx[index]*(type %in% c("failure", "relapse"))
    detectRR <- pop_dstR*c(0, 1, 0, 1)
    presumeRS <- pop - detectRR
    lab_costs <- lab_costs + pop_dstR*pars$dstR_c[index]
    
    #Novel DST step next
    detectRR_dstN <- detectRR*pars$dstN_RR[index, time]*(type %in% c("new", "initialloss")) + 
      detectRR*pars$dstN_RR_retx[index]*(type %in% c("failure", "relapse"))
    presumeRS_dstN <- presumeRS*pars$dstN_RS[index, time]*(type %in% c("new", "initialloss")) + 
      presumeRS*pars$dstN_RS_retx[index]*(type %in% c("failure", "relapse"))
    lab_costs <- lab_costs + (detectRR_dstN + presumeRS_dstN)*pars$dstN_c[index]
    
    #INH and FQ DST steps
    pop_dstH <- presumeRS*pars$dstH*pars$p_dstH[index] #whether there is INH DST and % that receive it
    pop_dstFQ <- (detectRR - detectRR_dstN*c(0, 0, 1, 1))*pars$dstFQ*pars$p_dstFQ[index]
    lab_costs <- lab_costs + pop_dstH*pars$dstH_c[index] + pop_dstFQ*pars$dstFQ_c[index]
    
    #resulting regimen assignment
    assign_RS <- (presumeRS - presumeRS_dstN*pars$RS_contains_N*c(0, 0, 1, 1))*(1-pars$pan) +
      presumeRS_dstN*(1-pars$RS_contains_N)*c(0, 0, 1, 1)*pars$pan
    assign_RR <- (detectRR - detectRR_dstN*c(0, 0, 1, 1))*(1-pars$pan)
    assign_P <- (presumeRS - presumeRS_dstN*(1-pars$RS_contains_N)*c(0, 0, 1, 1))*pars$pan +
      (detectRR - detectRR_dstN*c(0, 0, 1, 1))*pars$pan
    assign_I <- presumeRS_dstN*pars$RS_contains_N*c(0, 0, 1, 1) + detectRR_dstN*c(0, 0, 1, 1)
    
    #pre-treatment LTFU
    initiate_RS <- (presumeRS - presumeRS_dstN*pars$RS_contains_N*c(0, 0, 1, 1))*(1-pars$pan)*(1-pars$initialloss_base[[index]]) +
      presumeRS_dstN*(1-pars$RS_contains_N)*c(0, 0, 1, 1)*pars$pan*(1-(pars$initialloss_base[[index]] + pars$initialloss_extra[[index]]))
    initiate_RR <- assign_RR*(1-(pars$initialloss_base[[index]] + pars$initialloss_extra[[index]]))
    initiate_P <- assign_P*(1-pars$initialloss_base[[index]])
    initiate_I <- assign_I*(1-(pars$initialloss_base[[index]] + pars$initialloss_extra[[index]]))
  }
  
  if(pars$dstR_after_N) { #only used in Pan scenarios 
    #Novel DST step first
    pop_dstN <- pop*pars$dstN_RS[index, time]*(type %in% c("new", "initialloss")) + 
      pop*pars$dstN_RS_retx[index]*(type %in% c("failure", "relapse"))
    detectNR <- pop_dstN*c(0, 0, 1, 1)
    presumeNS <- pop - detectNR
    lab_costs <- lab_costs + pop_dstN*pars$dstN_c[index]
    
    #RR DST step next
    detectNR_dstR <- detectNR*pars$dstR[index, time]*(type %in% c("new", "initialloss")) + 
      detectNR*pars$dstR_retx[index]*(type %in% c("failure", "relapse"))
    lab_costs <- lab_costs + detectNR_dstR*pars$dstR_c[index]
    detectRR <- detectNR_dstR*c(0, 1, 0, 1)
    
    #INH and FQ DST steps
    pop_dstH <- (detectNR - detectRR)*pars$dstH*pars$p_dstH[index] #whether there is INH DST and % that receive it
    pop_dstFQ <- presumeNS*pars$dstFQ*pars$p_dstFQ[index]
    lab_costs <- lab_costs + pop_dstH*pars$dstH_c[index] + pop_dstFQ*pars$dstFQ_c[index]
    
    #resulting regimen assignment
    assign_RS <- (detectNR - detectRR)*(!pars$RS_contains_N)
    assign_RR <- 0
    assign_P <- presumeNS
    assign_I <- detectRR + (detectNR - detectRR)*(pars$RS_contains_N)
    
    #pre-treatment LTFU
    initiate_RS <- assign_RS*(1-(pars$initialloss_base[[index]] + pars$initialloss_extra[[index]]))
    initiate_RR <- assign_RR*(1-(pars$initialloss_base[[index]] + pars$initialloss_extra[[index]]))
    initiate_P <- assign_P*(1-pars$initialloss_base[[index]])
    initiate_I <- assign_I*(1-(pars$initialloss_base[[index]] + pars$initialloss_extra[[index]]))
  }
  
  #those not LTFU assigned to separate regimen pathway accrue additional costs
  oop_costs <- oop_costs + initiate_RR*pars$oopdx_extra_c[index]
  oop_costs <- oop_costs + initiate_I*pars$oopdx_extra_c[index]
  oop_costs <- oop_costs + pars$pan*initiate_RS*pars$oopdx_extra_c[index] #if pan scenario, RS SOC is a separate pathway too
  indirect_costs <- indirect_costs + initiate_RR*pars$indirectdx_extra_c[index]
  indirect_costs <- indirect_costs + initiate_I*pars$indirectdx_extra_c[index]
  indirect_costs <- indirect_costs + pars$pan*initiate_RS*pars$indirectdx_extra_c[index] #if pan scenario, RS SOC is a separate pathway too
  
  #of those that initiate treatment, a portion accrue up-front hospitalization costs
  inpatient_costs <- initiate_RS*pars$p_hosp_RS[index]*pars$hosp_RS_c[index] +
    initiate_RR*pars$p_hosp_RR[index]*pars$hosp_RR_c[index] +
    initiate_P*pars$p_hosp_P[index]*pars$hosp_P_c[index] +
    initiate_I*pars$p_hosp_I[index]*pars$hosp_I_c[index]
  indirect_costs <- indirect_costs +
    initiate_RS*pars$p_hosp_RS[index]*pars$indirecthosp_RS_c[index] +
    initiate_RR*pars$p_hosp_RR[index]*pars$indirecthosp_RR_c[index] +
    initiate_P*pars$p_hosp_P[index]*pars$indirecthosp_P_c[index] +
    initiate_I*pars$p_hosp_I[index]*pars$indirecthosp_I_c[index] 

  assignmat <- rbind(assign_RS, assign_RR, assign_P, assign_I)
  rownames(assignmat) <- c("RS", "RR", "P", "I")
  startmat <- rbind(initiate_RS, initiate_RR, initiate_P, initiate_I)
  rownames(startmat) <- c("RS", "RR", "P", "I")
  
  #add 0s for costs that don't accrue prediagnosis
  ae_costs <- rep(0, 4)
  support_costs <- rep(0, 4)
  drug_costs <- rep(0, 4)
  costs_out <- data.frame("outpatient"=outpatient_costs,
                          "inpatient"=inpatient_costs,
                          "labs"=lab_costs,
                          "aes"=ae_costs,
                          "support"=support_costs,
                          "drugs"=drug_costs,
                          "oop"=oop_costs,
                          "indirect"=indirect_costs)
  
  return(list(startmat, assignmat, costs_out))
}

#calculate tx initiators by regimen - used in resistance trends
calc_num_tx <- function(startmat, pars) {
  
  #number treated, by resistance and regimen (numerator)
  tx_R <- sum(startmat["RS", ]) #everyone treated getting rif
  tx_N_RS <- sum(startmat[c("RR", "P"), c("none", "N")]) + 
    pars$RS_contains_N*sum(startmat["RS", c("none", "N")]) #treated RS getting novel
  tx_N_RR <- sum(startmat[c("RR", "P"), c("R", "RN")]) + 
    pars$RS_contains_N*sum(startmat["RS", c("R", "RN")]) #treated RR getting novel
  tx_reg <- c(tx_R, tx_N_RS, tx_N_RR)
  
  #number treated, by resistance (denominator)
  tx_all <- sum(startmat) #all treated
  tx_RS <- sum(startmat[,c("none", "N")]) #rif susceptible treated
  tx_RR <- sum(startmat[, c("R", "RN")]) #rif resistant treated
  tx <- c(tx_all, tx_RS, tx_RR)
  
  return(list("tx_reg"=tx_reg, "tx"=tx))
}

#main treatment model
treatment_model <- function(startmat, pars, index=1) {
  with(pars, {
    # start from assignment to treatment regimens (after pre-tx LTFU and DST)
    
    # efficacy matrix based on resistance and regimen assignment
    efficacy_mat <- array(data=c(efficacy_RS[index], efficacy_RR[index], efficacy_P[index], efficacy_I[index],
                                 efficacy_RS[index]*RR_cure_R[index], efficacy_RR[index], efficacy_P[index], efficacy_I[index],
                                 efficacy_RS[index], efficacy_RR[index]*RR_cure_N[index], efficacy_P[index]*RR_cure_N[index], efficacy_I[index],
                                 efficacy_RS[index]*RR_cure_R[index], efficacy_RR[index]*RR_cure_N[index], efficacy_P[index]*RR_cure_N[index], efficacy_I[index]),
                          dim=c(4,4))
    if(RS_contains_N) {
      efficacy_mat[1,c(3,4)] <- efficacy_mat[1, c(3,4)]*RR_cure_N[index]
    }

    #multiply by regimen efficacy (= % cured)
    cured <- startmat*efficacy_mat
    
    #model LTFU (doesn't vary by regimen) over regimen duration
    #propagate LTFU by week across each regimen and resistance type
    durs <- c(dur_RS, dur_RR, dur_P, dur_I)
    cured <- outer(cured, ((1-ltfu[index])^(0:(max_dur-1)))) #cure prob multiplied by timing of LTFU
    #zero out rows that exceed the duration of the regimen
    durs_zero <- abind(array(rep(c(rep(1, durs[[1]]), rep(0, max_dur-durs[[1]])), 4), dim=c(max_dur, 4)),
                       array(rep(c(rep(1, durs[[2]]), rep(0, max_dur-durs[[2]])), 4), dim=c(max_dur, 4)),
                       array(rep(c(rep(1, durs[[3]]), rep(0, max_dur-durs[[3]])), 4), dim=c(max_dur, 4)),
                       array(rep(c(rep(1, durs[[4]]), rep(0, max_dur-durs[[4]])), 4), dim=c(max_dur, 4)),
                       along=3)
    durs_zero <- aperm(durs_zero, c(3,2,1))
    cured <- cured*durs_zero
    #subtract so that cured represents % that only took regimen for that long
    cured <- cured-abind(cured[,,2:max_dur], array(0, dim=c(4,4)), along=3)
    if(ltfu0_pan) { #no LTFU on pan regimen if this is turned on
      cured[3,,durs[[3]]+1] <- (startmat*efficacy_mat)[3,]
      cured[3,,1:durs[[3]]] <- 0 
    }
    #multiply by relative efficacy to calc % cured by time of LTFU
    cured <- model_ltfu(cured, durs, max_dur, as.vector(ltfu_rel_eff[index,]), ltfu_durs)
    cured <- rowSums(cured, dims=2) #no need to track dropout over time anymore
    
    #model adherence (varies by regimen)
    #adherence-efficacy relationship determined by forgiveness (varies by regimen)
    #first calculate % by # doses taken per week (while still on tx)
    #propogate adherence by group over each resistance type
    adherence <- array(data=c(rep(adherence_RS[index,], 4),
                              rep(adherence_RR[index,], 4), 
                              rep(adherence_P[index,], 4),
                              rep(adherence_I[index,], 4)),
                       dim=c(5, 4, 4))
    adherence <- aperm(adherence, c(3, 2, 1))
    cured <- abind(cured, cured, cured, cured, cured, along=3) * adherence
    #determine which adherence categories are poor vs. adequate - based on forgiveness
    forgive <- c(forgive_RS[index], forgive_RR[index], forgive_P[index], forgive_I[index])
    adequate_adhere <- match((1-forgive), adhere_grps)
    #next multiply by relative efficacy
    rel_eff <- array(data=rep(cbind(c(0, rep(rel_eff_forgive[index], adequate_adhere[1]-1), 
                                      rep(1, length(adhere_grps)-adequate_adhere[1])),
                                    c(0, rep(rel_eff_forgive[index], adequate_adhere[2]-1),
                                      rep(1, length(adhere_grps)-adequate_adhere[2])),
                                    c(0, rep(rel_eff_forgive[index], adequate_adhere[3]-1), 
                                      rep(1, length(adhere_grps)-adequate_adhere[3])),
                                    c(0, rep(rel_eff_forgive[index], adequate_adhere[4]-1),
                                      rep(1, length(adhere_grps)-adequate_adhere[4]))), 4),
                     dim=c(5, 4, 4))
    rel_eff <- aperm(rel_eff, c(2, 3, 1))
    cured <- cured*rel_eff
    #sum over adherence categories
    cured <- rowSums(cured, dim=2)

    return(cured)
  })
}

#model impact of LTFU over regimen duration - used within treatment_model function
model_ltfu <- function(cured, durs, max_dur, ltfu_rel_eff, ltfu_durs) { 
  # use relapse %s at 2, 4, and 6 months (defined in pars) for a 6 month regimen, and 100% relapse at 0, 
  # interpolate linearly 0-2, 2-4, and 4-6 to get relapse % after any fraction of treatment course.
  
  # first calculate % of duration completed corresponding to each timestep
  fracs <- cbind(c((1:(durs[1]-1)), durs[1], 
                   rep(durs[1], max_dur-durs[1])),
                 c((1:(durs[2]-1)), durs[2], 
                   rep(durs[2], max_dur-durs[2])),
                 c((1:(durs[3]-1)), durs[3], 
                   rep(durs[3], max_dur-durs[3])),
                 c((1:(durs[4]-1)), durs[4], 
                   rep(durs[4], max_dur-durs[4])))%*%diag(1/durs)
  fracs <- abind(fracs, fracs, fracs, fracs, along=3)
  fracs <- aperm(fracs, perm=c(2,3,1))
  
  # next calculate % efficacy achieved (and resulting % cured) based on % duration completed
  loops <- length(ltfu_rel_eff)-1
  for(i in 1:loops) {
    cured[fracs >= ltfu_durs[i] & fracs < ltfu_durs[i+1]] <-
      ((cured*ltfu_rel_eff[i])*(ltfu_durs[i+1] - fracs))[fracs >= ltfu_durs[i] & fracs < ltfu_durs[i+1]]/(1/(length(ltfu_durs)-1)) +
      ((cured*ltfu_rel_eff[i+1])*(fracs-ltfu_durs[i]))[fracs >= ltfu_durs[i] & fracs < ltfu_durs[i+1]]/(1/(length(ltfu_durs)-1)) 
  }
  return(cured)
}

#calculate costs and DALYs during tx based on params and health outcomes
calc_costs_dalys_tx <- function(startmat, pars, index=1, max_dur=76) { 
  with(pars, {
    #calculate tx time - same as in tx model but without the efficacy adjustment
    durs <- c(dur_RS, dur_RR, dur_P, dur_I)
    tx_time <- outer(startmat, ((1-ltfu[index])^(0:(max_dur-1)))) #cure prob multiplied by timing of LTFU
    if(ltfu0_pan) { #no LTFU on pan regimen if this is turned on
      tx_time[3,,1:durs[[3]]] <- startmat[3,]
    }
    #zero out rows that exceed the duration of the regimen
    durs_zero <- abind(array(rep(c(rep(1, durs[[1]]), rep(0, max_dur-durs[[1]])), 4), dim=c(max_dur, 4)),
                       array(rep(c(rep(1, durs[[2]]), rep(0, max_dur-durs[[2]])), 4), dim=c(max_dur, 4)),
                       array(rep(c(rep(1, durs[[3]]), rep(0, max_dur-durs[[3]])), 4), dim=c(max_dur, 4)),
                       array(rep(c(rep(1, durs[[4]]), rep(0, max_dur-durs[[4]])), 4), dim=c(max_dur, 4)),
                       along=3)
    durs_zero <- aperm(durs_zero, c(3,2,1))
    tx_time <- tx_time*durs_zero
    #manipulate tx_time to be a 16 x max_dur matrix (rows = regimens x resistance types: panS on RS, panS on RR, panS on P, etc.)
    dim(tx_time) <- c(16, max_dur)
    
    #costs with a schedule (depend on ppl still on tx each)
    outpatient_costs <- tx_time*outpatient_sched*outpatient_c[index]
    lab_costs <- tx_time*(smear_sched*smear_c[index] +
                            culture_sched*culture_c[index] + 
                            cxr_sched*cxr_c[index] + 
                            livertest_sched*livertest_c[index] +
                            fullblood_sched*fullblood_c[index] + 
                            ecg_sched*ecg_c[index] +
                            neuroscreen_sched*neuroscreen_c[index])
    support_costs <- tx_time*support_sched[[index]]
    oop_costs <- tx_time*ooptx_sched[[index]]
    indirect_costs <- tx_time*indirecttx_sched[[index]]
    drug_costs <- tx_time*drugs_sched[[index]]*(1+wastage[index]) #add wastage to all drug components
    
    #once costs are calculated, sum over weeks (no need to track time anymore)
    outpatient_costs <- rowSums(outpatient_costs)
    lab_costs <- rowSums(lab_costs)
    support_costs <- rowSums(support_costs)
    drug_costs <- rowSums(drug_costs)
    oop_costs <- rowSums(oop_costs)
    indirect_costs <- rowSums(indirect_costs)
    
    #AEs have weekly incidence - can just sum person time (no need to do this by week)
    ae_costs <- rowSums(tx_time)*(p_liverAE[index,]*liverAE_c[index] + 
                                    p_pancreasAE[index,]*pancreasAE_c[index] +
                                    p_anemiaAE[index,]*anemiaAE_c[index] +
                                    p_neutropeniaAE[index,]*neutropeniaAE_c[index] +
                                    p_qtcfAE[index,]*qtcfAE_c[index] +
                                    p_renalAE[index,]*renalAE_c[index] +
                                    p_visionAE[index,]*visionAE_c[index] +
                                    p_arthralgiaAE[index,]*arthralgiaAE_c[index] +
                                    p_neuroshortAE[index,]*neuroAE_c[index])
    
    inpatient_costs <- rep(0, 16) #inpatient costs are all up front, but include so we can add later
    
    costs_out <- data.frame("outpatient"=outpatient_costs,
                            "inpatient"=inpatient_costs,
                            "labs"=lab_costs,
                            "aes"=ae_costs,
                            "support"=support_costs,
                            "drugs"=drug_costs,
                            "oop"=oop_costs,
                            "indirect"=indirect_costs)
    
    #also calculate DALYs incurred during treatment (just AEs)
    #all AEs last for 1 month except short-term neuropathy (3 months) - happening up front so no need to discount
    dalys_ae <- rowSums(tx_time)*(p_liverAE[index,]*dalys_liver[index]/12 + 
                                    p_pancreasAE[index,]*dalys_pancreas[index]/12 +
                                    p_anemiaAE[index,]*dalys_anemia[index]/12 +
                                    p_neutropeniaAE[index,]*dalys_neutropenia[index]/12 +
                                    p_qtcfAE[index,]*dalys_qtcf[index]/12 +
                                    p_renalAE[index,]*dalys_renal[index]/12 +
                                    p_visionAE[index,]*dalys_vision[index]/12 +
                                    p_arthralgiaAE[index,]*dalys_arthralgia[index]/12 +
                                    p_neuroshortAE[index,]*dalys_neuroshort[index]/4)
    
    #assume only those on tx for at least a month incur long-term disability from it
    dalys_ae <- dalys_ae + tx_time[,5]*p_neurolongAE[index,]*dalys_neurolong[index] #already discounted
    
    #only need DALYs by resistance type, not regimen assignment
    dalys_ae_sum <- c(sum(dalys_ae[1:4]), sum(dalys_ae[5:8]), sum(dalys_ae[9:12]), sum(dalys_ae[13:16]))
    
    return(list("costs_out"=costs_out, "dalys_ae"=dalys_ae_sum))
  })
}


