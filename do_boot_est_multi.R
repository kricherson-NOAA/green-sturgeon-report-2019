do_boot_multi <- function(ob_dat, ft_dat, strata, expfactor = tgt_mt, bycatchspp, bycatchunit = dis_mt, seed = 42)
{
  require(dplyr)
  require(tidyr)
  require(rlang)
  require(janitor)
  require(boot)
  
  set.seed(seed)
  
  #Identify whether the input data is ASHOP or WCGOP (won't ever be bootstrapping ASHOP but leaving this in for now)
  if(any(ob_dat$sector %in% c("CATCHER-PROCESSOR", "MOTHERSHIP", "ATSEA TRIBAL"))){program <- "ASHOP"} else {
    program <- "WCGOP"
  }
  
  #Do a little bit of data prep. Maybe would be better NOT to have clean_names in here?
  if(program == "WCGOP")
  {
    ob_data <- clean_names(ob_dat) %>% 
      filter(datatype == "Analysis Data")
  }
  
  #For ASHOP, we put in NAs for trip_id and just remove later. 
  if(program == "ASHOP")
  {
    ob_data <- clean_names(ob_dat) %>% 
      mutate(tgt_mt = pwht_mt) %>% 
      mutate(trip_id = NA)
  }
  
  #coverage and ratio
  #First, summarise the expansion factor data by strata.
  obdf_expf <- ob_data %>% 
    filter(!is.na(haul_id)) %>% #removes any trip without catch
    #Haul-level summary
    group_by_at(c(strata, "haul_id", "trip_id", "drvid")) %>% 
    summarise(expf_haul = sum(!!sym(expfactor), na.rm=T)) %>%
    ungroup() %>%
    #Get the expansion factor by strata
    group_by_at(strata) %>% 
    summarise(total_expf = sum(expf_haul, na.rm=T),
              mean_expf = mean(expf_haul, na.rm=T),
              se_expf = sqrt(var(expf_haul, na.rm=T)/length(expf_haul)),
              n_obs_ves = n_distinct(drvid),
              n_obs_trips = n_distinct(trip_id), 
              n_obs_hauls = n_distinct(haul_id)) %>% 
    #Now cobble together a column that contains the bycatch species so that we can join bycatch info without losing any strata with 0 bycatch (maybe not needed? but for completeness...). Clumsy but ok I guess.
    mutate(all_byc_species = paste(bycatchspp, sep="", collapse=",")) %>% 
    tidyr::separate(all_byc_species, bycatchspp, sep = ",") %>% 
    tidyr::gather(species, whatev, all_of(bycatchspp)) %>% 
    dplyr::select(-whatev)
  
  
  #Get unique combination of bycatch species, hauls, to make sure we account for hauls with 0 bycatch below. This seemed simpler than dealing with nested tidyr::expand.
  unique_hauls_sp <- ob_data %>% 
    filter(!is.na(haul_id)) %>% #removes any trip without catch
    dplyr::select(!!!syms(strata), haul_id, trip_id, drvid) %>% 
    distinct() %>% 
    #This is a hacky way of making sure each row contains strata, trip, haul, drvid, and bycatch species
    mutate(all_byc_species = paste(bycatchspp, sep="", collapse=",")) %>% 
    tidyr::separate(all_byc_species, bycatchspp, sep = ",") %>% 
    tidyr::gather(species, whatev, all_of(bycatchspp)) %>% 
    dplyr::select(-whatev)
  
  #Now summarise bycatch data by strata
  obdf_byc <- ob_data %>% 
    #Haul-level summary  
    #Filter to only species of interest
    filter(species %in% bycatchspp) %>% #
    #Joining all the unique haul/bycatch species combinations ensures that hauls with 0 catch are preserved.
    right_join(unique_hauls_sp) %>% 
    group_by_at(c(strata, "species", "haul_id", "trip_id", "drvid")) %>% 
    ###summarise(byc_haul = sum(!!byc_unit, na.rm=T)) %>% 
    summarise(byc_haul = sum(!!sym(bycatchunit), na.rm=T)) %>% 
    #byc = sum(dis_mt, na.rm = T)) %>% #!! unquotes the quosure
    ungroup() %>% 
    #Fill in 0s for strata without observed bycatch, for completeness. "Nesting" should return only strata combinations present in the data. 
    #tidyr::complete(tidyr::nesting(!!!syms(strata)), spid_eqv, fill = list(byc = 0)) %>% 
    #Summarise the observed bycatch by strata, including the number of hauls that encountered bycatch
    group_by_at(c(strata, "species")) %>% 
    summarise(total_byc = sum(byc_haul),
              mean_byc = mean(byc_haul),
              se_byc = sqrt(var(byc_haul)/length(byc_haul)),
              n_hauls_byc = n_distinct(haul_id[byc_haul > 0]))
  
  obdf <- obdf_expf %>% 
    full_join(obdf_byc) %>% 
    ungroup() %>% 
    mutate(byc_ratio = total_byc / total_expf,
           #byc_se = ratio.se(mean_byc, mean_expf, se_byc, se_expf),
           #byc_ratio_lower = byc_ratio - 1.96 * byc_se, 
           #byc_ratio_upper = byc_ratio + 1.96 * byc_se,
           pct_hauls_byc = round((n_hauls_byc/n_obs_hauls)*100, 2))
  
  #Bootstrap portion
  ob_byvessel_expf <- ob_data %>% 
    group_by_at(c(strata, "drvid")) %>% 
    summarise(ves_expf = sum(!!sym(expfactor), na.rm=T)) %>% 
    #Now cobble together a column that contains the bycatch species so that we can join bycatch info without losing any strata with 0 bycatch (maybe not needed? but for completeness...). Clumsy but ok I guess.
    mutate(all_byc_species = paste(bycatchspp, sep="", collapse=",")) %>% 
    tidyr::separate(all_byc_species, bycatchspp, sep = ",") %>% 
    tidyr::gather(species, whatev, all_of(bycatchspp)) %>% 
    dplyr::select(-whatev)
  
  ob_byvessel_byc <- ob_data %>% 
    filter(species %in% bycatchspp) %>% 
    group_by_at(c(strata, "drvid", "species")) %>% 
    summarise(ves_byc = sum(!!sym(bycatchunit), na.rm=T))
  
  ob_byvessel <- full_join(ob_byvessel_expf, ob_byvessel_byc, by = c(strata, "drvid", "species")) %>% 
    mutate(ves_byc = ifelse(is.na(ves_byc), 0, ves_byc))

  #Get FT data by vessel
  ft_byvessel <- ft_dat %>%
    # mutate(tgt_mt = case_when(sector %in% c("Limited Entry Trawl", "Midwater Rockfish", "Catch Shares", "OA Fixed Gear", "LE Fixed Gear DTL") ~ gfr_mt,
    #                           sector ==  "Limited Entry Sablefish" ~ sabl_mt,
    #                           sector %in% c("LE CA Halibut", "OA CA Halibut") ~ chlb_mt,
    #                           sector == "Nearshore" ~ ns_mt,
    #                           sector == "Pink Shrimp" ~ ps_mt,
    #                           sector == "Directed P Halibut" ~ phlb_mt,
    #                           sector == "Sea Cucumber" ~ cuke_mt,
    #                           sector == "Ridgeback Prawn" ~ prwn_mt,
    #                           sector %in% c("Shoreside Hake", "Midwater Hake") ~ pwht_mt)) %>% 
    group_by_at(c(strata, "drvid")) %>%  
    summarise(ves_ft_expf = sum(!!sym(expfactor), na.rm=T)) %>% 
    group_by_at(strata) %>%  
    mutate(fleet_ft_expf = sum(ves_ft_expf, na.rm=T))
    
  #This function feeds into boot::boot as the statistic argument. It takes a data frame of vessel-level data and calculates the bycatch ratio according to the vessels sampled by the bootstrapping routine
  byc_ratio_fun <- function(vdata, indices){ #vdata is the vessel-level data, indices will be the rows sampled by the boot function
    
    sum(vdata$ves_byc[indices]) / sum(vdata$ves_expf[indices])
    
  }
  
  booted <- ob_byvessel %>% 
    left_join(ft_byvessel, by=c(strata, "drvid")) %>% 
    filter(ves_expf>0) %>% #So that we don't end up dividing by 0 -- there are two vessel/season/state combos with 0 retained groundfish (and 0 gstg)
    filter(!is.na(ves_ft_expf)) %>% #Because of missing LE FTs, there are some NAs. Ignore these few cases (3) where a vessel apparently had observed landings, but no FT landings. There could be a better way of doing things?
    group_by_at(c(strata, "species"))%>% 
    nest() %>% #This creates a row for each year/state/season, with a data column that contains a df with all the data for that year/state/season
    mutate(booted = map(.x = data, # The list-column with the vessel-level data
                        ~ boot(data = .x, # The <S3 tibble> being sampled
                               statistic = byc_ratio_fun, # The user-defined function
                               R = 10000, # The number of replicates
                               stype = "i"))) %>% 
    rowwise() %>% #"rowwise() is used for the results of do() when you create list-variables." (presumably similar for map()-created list variables)
    mutate(point_est_ratio = booted$t0,  #Get point estimate, median, lower 2.5th percentile, upper 97.5th percentile of bootstrapped ratios
           median_ratio = median(booted$t),
           lower_ci = quantile(booted$t, 0.025),
           upper_ci = quantile(booted$t, 0.975)) %>% 
    # Drop the list-columns (no longer needed)
    dplyr::select(-data, -booted) %>%
    # Unnest the dataframe
    tidyr::unnest()
  
  #create expansion factors data frame
  exp_factors <- ft_dat %>% 
    mutate(r_state = state) %>% 
    group_by_at(strata) %>% 
    summarise(fleet_expf = sum(!!sym(expfactor), na.rm=T))
  
  out <- exp_factors %>% 
    full_join(booted, by = strata) %>% 
    full_join(obdf_expf, by = c(strata, "species")) %>% 
    full_join(obdf_byc, by = c(strata, "species")) %>% 
    mutate(pct_cvg = round((total_expf / fleet_expf) * 100, 2),
           pct_hauls_byc = round((n_hauls_byc/n_obs_hauls)*100, 2),
           byc_ratio = total_byc / total_expf,
           est_byc = byc_ratio * fleet_expf,
           est_byc_lower = lower_ci* fleet_expf,
           est_byc_upper = upper_ci * fleet_expf) %>% 
    ungroup()
  
 return(out) 
  
}