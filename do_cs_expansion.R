#inputs: df (ob_cs data frame), bycatchsp (bycatch species of interest), bycatchunit(column name with bycatch unit, eg gstg_count)

# df<-ob_cs
# strata<-c("year","r_state")
# bycatchspp <-c("GSTG", "EULC")
# save(df, file = "df.Rdata")

do_cs_expansion <- function(df, strata, bycatchspp)
  {
  #Identify hauls with NIFQ
  df_nifq <- df %>% 
    dplyr::filter(datatype == "Unsampled IFQ"
               & catch_category_code == "NIFQ"
               & catch_disposition == "D")
  
  hauls_nifq <- df[df$haul_id %in% unique(df_nifq$haul_id),]
  
  #Make an empty list that will get filled with expansions for each species
  cs_sp_expansions <- list()
    
  #Loop through species we're expanding. Doing this all in tidyverse was slow, let's see if this is faster
  for (s in 1:length(bycatchspp))
  {
    bycatchsp <- bycatchspp[s]
    
    #For salmon and green sturgeon, identify NIFQ hauls that need expansion. These are hauls where species is present in strata, but not listed as species specific discard in a given haul. 
    if(bycatchsp %in% c("CHNK", "CHUM", "COHO", "PINK", "STLH", "SOCK", "USMN", "GSTG"))
    {
       nifq_hauls_to_expand <- hauls_nifq %>% 
         distinct(haul_id, spid_eqv, catch_disposition) %>% 
         group_by(haul_id) %>% 
         mutate(prot_in_haul = ifelse(paste0(bycatchsp, "D") %in% paste0(spid_eqv, catch_disposition), TRUE, FALSE)) %>%  #This identifies hauls within the hauls with NIFQ that have discareded protected species
         filter(!prot_in_haul) #filter to only hauls without discarded protected species
    }
    
    # For eulachon, always expand if EULC occurs in the strata, even if EULC is also listed in DISCARD of a given haul. Note that you should include this explanation in any reports for EULC and highlight that this estimate is conservative due to the potential over expansion.
    
    if(bycatchspp[s] == "EULC")
    {
      nifq_hauls_to_expand <- hauls_nifq 
    }
    
    #identify failed trips
    failed_trips <- unique(df$trip_id[df$datatype == "Failed Data"])
    
    #Now create DF with numerators, denominators, and expansion factors. Note conversion to MT
    unsamp_expansion <- df %>% 
      group_by_at(strata) %>% 
      summarise(nifq_num_ct = sum(exp_sp_ct[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]),
                
                nifq_num_wt = sum(exp_sp_wt[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]/2204.6),
                
                nifq_denom = sum(mt[ifq == 0 &
                                      catch_disposition == "D"]),
                
                nifq_expf = sum(mt[ifq == 0 &
                                     datatype == "Unsampled IFQ" &
                                     catch_category_code == "NIFQ" &
                                     catch_disposition == "D" &
                                     haul_id %in% nifq_hauls_to_expand$haul_id]),
                
                zmis_num_ct = sum(exp_sp_ct[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]),
                
                zmis_num_wt = sum(exp_sp_wt[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]/2204.6),
                
                zmis_denom = sum(mt[datatype == "Analysis Data" &
                                      catch_disposition == "D"]),
                
                zmis_expf = sum(mt[datatype == "Unsampled ZMIS" &
                                     catch_category_code == "ZMIS" &
                                     species_composition_id == 0 &
                                     catch_disposition == "D"]),
                
                unst_num_ct = sum(exp_sp_ct[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]),
                
                unst_num_wt = sum(exp_sp_wt[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]/2204.6),
                
                unst_denom = sum(mt[datatype == "Analysis Data"]),
                
                unst_expf = sum(mt[datatype == "Unsampled ZMIS" &
                                     catch_category_code == "UNST" &
                                     species_composition_id == 0 &
                                     catch_disposition == "D"]),
                
                fail_num_ct = sum(exp_sp_ct[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]),
                
                fail_num_wt = sum(exp_sp_wt[spid_eqv == bycatchsp &
                                              datatype == "Analysis Data" &
                                              catch_disposition == "D"]/2204.6),
                
                fail_denom = sum(tgt_mt[datatype == "Analysis Data"]),
                
                fail_expf = sum(tgt_mt[trip_id %in% failed_trips])
                
                ) %>% 
      mutate(est_nifq_ct = (nifq_num_ct / nifq_denom) * nifq_expf,
             est_nifq_wt = (nifq_num_wt / nifq_denom) * nifq_expf,
             est_zmis_ct = (zmis_num_ct / zmis_denom) * zmis_expf,
             est_zmis_wt = (zmis_num_wt / zmis_denom) * zmis_expf,
             est_unst_ct = (unst_num_ct / unst_denom) * unst_expf,
             est_unst_wt = (unst_num_wt / unst_denom) * unst_expf,
             est_fail_ct = (fail_num_ct / fail_denom) * fail_expf,
             est_fail_wt = (fail_num_wt / fail_denom) * fail_expf) %>% 
      mutate_at(vars(-group_cols()), ~replace(., is.nan(.), 0)) %>% #when dividing by 0 we get NANs that need to be replaced
      mutate(est_unsamp_ct = est_nifq_ct + est_zmis_ct + est_unst_ct + est_fail_ct,
             est_unsamp_wt = est_nifq_wt + est_zmis_wt + est_unst_wt + est_fail_wt,
             bycsp = bycatchsp)
    
    cs_sp_expansions[[s]] <- unsamp_expansion
  
    
  } #end species loop
  
  out <- bind_rows(cs_sp_expansions)
  
  return(out)
}