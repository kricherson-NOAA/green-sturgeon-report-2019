#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Green sturgeon report code
#------------------------------------------------------------
# Date          Modification                       Author
#------------------------------------------------------------
#  2018 Oct   Original Coding & Development       Kate Richerson
#             (Drawing from old code by YWL)

#  2019 Nov   Updating with new data              Kate Richerson
#------------------------------------------------------------

##%######################################################%##
#                                                          #
####   Load data and libraries, set drives, date etc    ####
#                                                          #
##%######################################################%##

rm(list=ls())

library(tidyverse)
library(janitor)
library(boot)
library(viridis)
#library(broom)

#For environmental data/gams:
library(raster) #Extracting gridded data at given locations
library(ncdf4) # .nc SST files
library(mgcv) #GAMs
library(mgcViz) #Plotting GAMs. Note that this package isn't compatible with version of R currently on Tantalus...
library(lubridate) #Julian date etc
library(marmap) #Bathymetry

source("~/observer/Input/load_data_2020-06-24.R") #

source("do_ratio_multi.R") #does bycatch ratios
source("do_ratio_est_multi.R") #does ratio expansion

load_data(c("WCGOP","ASHOP_proc","FT_proc", "EM", "BIO"))

tday <- gsub("-", "", Sys.Date()) #Current date

out_drive <- "~/observer/Output/Richerson reports/2019/Green sturgeon/"

in_drive <- "~/observer/Input/Richerson/" #This only has the GSI data for now (no other inputs)

set.seed(42) #This SHOULD make bootstrap results reproducible, I think

##%######################################################%##
#                                                          #
####             Initial data manipulation              ####
#                                                          #
##%######################################################%##
#Modify names of data frames for ease of use, add some relevant columns, do some initial data exploration (e.g. which sectors have GSTG and when?)
ob <- OBOrig_Proc %>% 
  clean_names() %>% 
  mutate(gstg_count = ceiling(ifelse((spid_eqv == "GSTG" | species == "Green Sturgeon") & catch_disposition == "D",
                                     exp_sp_ct, 0)), #Add column for expanded discard count. Rounding up following YWL.
         season = ifelse(rmonth %in% 5:10,
                         "summer", "winter"),
         state = r_state)  #Add column for season, following YWL
#We check for retained GSTG below

ft <- FTOrig_Proc %>% 
  clean_names() %>% 
  mutate(state = case_when(agency_code == "C" ~ "CA",
                           agency_code == "O" ~ "OR",
                           agency_code == "W" ~ "WA"), 
         season = ifelse(landing_month %in% 5:10,
                         "summer", "winter")) #Add column for season, followingt YWL

ashop <- ASOrig_Proc %>% 
  clean_names() %>% 
  mutate(gstg_count = ceiling(ifelse((spid_eqv == "GSTG" | species == "Green Sturgeon") & !is.na(spid_eqv), #Note there are some entries with NA for species
                                     expanded_num_2sector_level, 0))) %>%  #Add column for expanded discard count.
  filter(sector != "ATSEA TRIBAL")

#Check whether there are any EM landings. I also checked with PFMC and they have not seen any on video review
"Green Sturgeon" %in% EMOrig_Proc$spc.name | "Acipenser medirostris" %in% EMOrig_Proc$sci.name | "Green Sturgeon" %in% EMOrig_Proc$sci.name
#FALSE

#Which sectors have GSTG observed?
ob %>% 
  filter(gstg_count>0) %>% 
  group_by(sector, gear) %>% 
  summarise(total_gstg=sum(gstg_count),
            min_year=min(year),
            max_year=max(year))

# sector              gear         total_gstg min_year max_year
# <chr>               <chr>             <dbl>    <dbl>    <dbl>
# 1 Catch Shares        Bottom Trawl        140     2011     2017
# 2 Directed P Halibut  Hook & Line           1     2019     2019
# 3 LE CA Halibut       Bottom Trawl        443     2002     2013
# 4 Limited Entry Trawl Bottom Trawl         19     2002     2010
# 5 Nearshore           Fixed Gears           1     2017     2017
# 6 OA CA Halibut       Bottom Trawl        842     2003     2019

gstg_sectors <- unique(filter(ob, spid_eqv == "GSTG")$sector)

#Check for retained (dockside discards OR retained for human consumption back when it was legal. Supposedly this stopped in 2006, but we still have 2 in 2009?)
ob %>% 
  filter(spid_eqv == "GSTG" & catch_disposition == "R") %>% 
  group_by(sector, year, r_state) %>% 
  summarise(n=n())
# sector               year r_state     n
# <chr>               <int> <chr>   <int>
# 1 Limited Entry Trawl  2002 OR         29
# 2 Limited Entry Trawl  2002 WA          1
# 3 Limited Entry Trawl  2003 OR          2
# 4 Limited Entry Trawl  2004 OR         27
# 5 Limited Entry Trawl  2004 WA          1
# 6 Limited Entry Trawl  2005 OR          1
# 7 Limited Entry Trawl  2005 WA          2
# 8 Limited Entry Trawl  2006 OR          2
# 9 Limited Entry Trawl  2009 OR          2

#Pull out Catch Shares OB data. Starting with OBOrig_Pre and processing it to get the unsampled data. Following KS's methods from GM_2017 and YWL's old GSTG code. 

ob_cs_pre <- OBOrig_Pre %>% 
  filter(PROGRAM_ID == 14 &
           DATATYPE %in% c('Analysis Data', 'Unsampled IFQ', 'Unsampled ZMIS', 'Failed Data') &
           MMSBT == 0 &
           sector == "Catch Shares") #Since we already know GSTG was only observed in the CS trawl fishery

ob_ifq <- OB.processing(ob_cs_pre, 'Catch Shares', SPC)

######## LE California Halibut 2011-2013
#ob_chlb_cs <- OB.processing(ob_cs_pre, 'LE CA Halibut', SPC)
#Note: It looks like YWL *did not* include LE CHLB with CS. Maybe because it is not covered under the BiOp? Leave out for now to be consistent.

#Note: we don't need the following sectors (SSH, MWH, MWR) since no GSTG has been observed in them, but leaving the commented-out code in case we need it later
######## Shoreside Hake 2011-2014
#note from KS's GM code: Note that this fishery "ends" in 2014 because the definition changed from being based on logbook target
#   to being based on landings. 2015-forward, the fishery is Midwater Hake and Midwater Rockfish
#ob_ssh <- OB.processing(ob_cs_pre, 'Shoreside Hake', SPC)

######## Midwater Hake 2015-forward
#ob_mwh <- OB.processing(ob_cs_pre, 'Midwater Hake', SPC)

######## Midwater Rockfish 2015-forward
#ob_mwr <- OB.processing(ob_cs_pre, 'Midwater Rockfish', SPC)

#### combine all processed CS data
ob_cs <- ob_ifq %>% #bind_rows(ob_ifq, ob_ssh, ob_mwh, ob_mwr)
  mutate(cs_gear = gear.type(.)) %>% 
  clean_names() %>% 
  filter(catch_disposition != '') %>%  #KS code says: remove any hauls w/no catch to avoid issues with stratification - ie, no depth
  mutate(cs_sector = paste0("CS - ",cs_gear)) %>%  #Note this makes LE CHLB trawl "CS - Bottom Trawl" if it's included
  mutate(gstg_count = ceiling(ifelse(spid_eqv == "GSTG" & catch_disposition == "D", #Add in column for GSTG counts
                                     exp_sp_ct, 0)))

#save(ob_cs, file = "ob_cs.Rdata")
#load("ob_cs.Rdata")

#GSTG in fish tickets. This is pretty much all for human consumption, so doesn't count as bycatch, I believe
ft %>% 
  filter(spid == "GSTG") %>% 
  group_by(sector, year, state, product_use_code) %>% 
  summarise(n=n()) %>% 
  as.data.frame()
# sector year state product_use_code  n
# 1                  EFP 2003    OR                H  3
# 2  Limited Entry Trawl 2002    OR                H 13
# 3  Limited Entry Trawl 2002    WA                H  4
# 4  Limited Entry Trawl 2003    OR                H 13
# 5  Limited Entry Trawl 2004    OR                H 28
# 6  Limited Entry Trawl 2004    WA                H  3
# 7  Limited Entry Trawl 2005    OR                H 25
# 8  Limited Entry Trawl 2005    WA                H  6
# 9  Limited Entry Trawl 2006    OR                H 22
# 10 Limited Entry Trawl 2006    WA                H  1
# 11 Limited Entry Trawl 2007    OR                H  4
# 12 Limited Entry Trawl 2008    OR                H  2
# 13 Limited Entry Trawl 2009    OR                H  6
# 14     Other Fisheries 2004    WA                H  1
# 15            Research 2003    OR                D  1
# 16              Tribal 2003    WA                H  3
# 17              Tribal 2005    WA                H  3
# 18              Tribal 2007    WA                H  2
# 19              Tribal 2010    WA                H  2
# 20              Tribal 2011    WA                H  1
# 21              Tribal 2012    WA                H  1
# 22              Tribal 2013    WA                H  1

#And what about retained, observed?

ob %>% 
  filter(spid_eqv == "GSTG" & catch_disposition == "R") %>% 
  group_by(sector, year, r_state) %>% 
  summarise(n=n()) %>% 
  as.data.frame()

# sector year r_state  n
# 1 Limited Entry Trawl 2002      OR 29
# 2 Limited Entry Trawl 2002      WA  1
# 3 Limited Entry Trawl 2003      OR  2
# 4 Limited Entry Trawl 2004      OR 27
# 5 Limited Entry Trawl 2004      WA  1
# 6 Limited Entry Trawl 2005      OR  1
# 7 Limited Entry Trawl 2005      WA  2
# 8 Limited Entry Trawl 2006      OR  2
# 9 Limited Entry Trawl 2009      OR  2

bio <- BIOOrig_Proc %>% 
  clean_names()

#Load GSI data
# gsi <- read_csv(paste0(in_drive, "gsi_results.csv")) %>% 
#   clean_names() %>% 
#   mutate(assigned_to = tolower(assigned_to)) %>% 
#   mutate(dps = case_when(grepl("southern", assigned_to) ~ "sdps",
#                          grepl("northern", assigned_to) ~ "ndps",
#                          grepl("southern", assigned_to) & grepl("southern", assigned_to) ~ "check notes", #In case there's something funky
#                          TRUE ~ "Unknown")) #Ditto

#End inital data processing section
#Next sections will create bycatch tables by sector

##%######################################################%##
#                                                          #
####          LE, PHLB, NS expansions together          ####
#                                                          #
##%######################################################%##

#can do these three sectors together because they use the same strata and output format (could also consider doing PHLB coastwide...)

ncs_gstg <- do_boot_multi(ob_dat = filter(ob, sector %in% c("Limited Entry Trawl", "Nearshore", "Directed P Halibut")), #observer data
                         ft_dat = filter(ft, sector %in% c("Limited Entry Trawl", "Nearshore", "Directed P Halibut")), #fish tickets data
                         strata = c("sector","year", "state", "season"), #strata that we want estimates for
                         expfactor = "tgt_mt", #expansion factor
                         bycatchspp = "Green Sturgeon", #species of interest (from the "species" column)
                         bycatchunit = "gstg_count") #the unit of bycatch that is being estimated


ncs_gstg_out <- ncs_gstg %>% 
  dplyr::select(sector, state, year, season, total_byc, total_expf, fleet_expf, pct_cvg, byc_ratio, lower_ci, upper_ci, est_byc, n_obs_ves, est_byc_lower, est_byc_upper) %>%
  ungroup() %>% 
  mutate(est_byc_lower = ifelse(est_byc_lower < total_byc, total_byc, est_byc_lower)) %>% #Truncate lower CI at observed value
  #mutate(gr = ifelse(est_byc_lower < total_byc, total_byc, est_byc_lower)) %>% #Truncate lower CI at observed value
  mutate_at(c("est_byc", "est_byc_lower", "est_byc_upper"), round, 0) %>%  #Round counts to whole numbers
  mutate_at(c("total_expf", "fleet_expf", "pct_cvg"), round, 1) %>%  #Round weights/percents to 1 decimal place
  mutate_at(c("byc_ratio", "lower_ci", "upper_ci"), round, 2) %>%  #Round ratios to 2 decimal places
  mutate_at(c("lower_ci", "upper_ci", "est_byc_lower", "est_byc_upper"), list(~ifelse(total_byc == 0, NA, .))) %>% #If no bycatch, replace 0s with NAs here
  mutate_at(c("total_byc", "total_expf", "pct_cvg", "lower_ci"), as.character) %>% #Make into character so we can replace confidential data with * and make 0 lower CIs display as 0.00
  mutate_at(c("total_byc", "total_expf", "pct_cvg", "lower_ci","upper_ci", "est_byc_lower", "est_byc_upper","byc_ratio"), as.character) %>% #Make into character so we can replace confidential data with * and make 0 lower CIs display as 0.00
  mutate(lower_ci = ifelse(lower_ci == "0", "0.00", lower_ci)) %>% #For display purposes
  mutate(byc_ratio = ifelse(byc_ratio == "0" & total_byc!= "0", "0.00", byc_ratio)) %>% #For display purposes (when rounded bycatch ratio is v. small)
  mutate_at(c("byc_ratio", "est_byc"), list(~ifelse(is.na(total_byc), "-", .))) %>% #If strata not observed, replace NAs with dashes here
  mutate_at(c("lower_ci", "upper_ci", "est_byc_lower", "est_byc_upper"), list(~ifelse(total_byc == "0" | is.na(total_byc), "-", .))) %>% #If no bycatch OR strata not observed, replace 0s or NAs with dashes here
  mutate_at(c("total_byc", "total_expf", "pct_cvg"), list(~ifelse( n_obs_ves < 3, "*", .))) %>% #Replace confidential data with *
  dplyr::select(-n_obs_ves) %>% 
  arrange(sector, desc(state),year,desc(season)) %>% 
  rename(State= state, 
         Year= year, 
         Season = season, 
         `Observed bycatch` = total_byc, 
         `Observed target landings (MT)` = total_expf, 
         `Fleet-total target landings (MT)`= fleet_expf, 
         `Target landings sampled (%)` = pct_cvg, 
         `Bycatch ratio` = byc_ratio, 
         `Lower CI of ratio` = lower_ci,
         `Upper CI of ratio` = upper_ci, 
         `Fleet total bycatch` = est_byc, 
         `Lower CI of bycatch`= est_byc_lower,
         `Upper CI of bycatch` = est_byc_upper)

#Separate into different sectors for output. PHLB needs a little more editing since it was not observed until 2017 (and no CA observations in 2017)
### PHLB ###
phlb_gstg_out <- ncs_gstg_out %>% 
  filter(sector == "Directed P Halibut" & Year > 2016) %>% 
  dplyr::select(-sector) %>% 
  mutate_at(c("Observed bycatch", "Observed target landings (MT)"), list(~ifelse(State == "CA" & Year ==2017, "-", .))) %>% 
  mutate_at(c( "Target landings sampled (%)"), list(~ifelse(State == "CA" & Year ==2017, 0, .)))  

write_csv(phlb_gstg_out, paste0(out_drive, "phlb_gstg_", tday, ".csv"))

### NS ###
ns_gstg_out <- ncs_gstg_out %>% 
  filter(sector == "Nearshore" & !is.na(`Observed bycatch`) & State == "CA") %>% 
  dplyr::select(-sector)

write_csv(ns_gstg_out, paste0(out_drive, "ns_gstg_", tday, ".csv"))


### LE ###
le_gstg_out <- ncs_gstg_out %>% 
  filter(sector == "Limited Entry Trawl") %>% 
  dplyr::select(-sector)

write_csv(le_gstg_out, paste0(out_drive, "le_gstg_", tday, ".csv"))

##%######################################################%##
#                                                          #
####                 Catch shares trawl                 ####
#                                                          #
##%######################################################%##
#Subset to CS bottom trawl data (only non-hake CS sector that encounters GSTG)
ob_cs_bt <- ob_cs %>% 
  filter(cs_sector == "CS - Bottom Trawl")

###############Expand to unsampled data####################

#Which kinds of unsampled data do we have? Failed data, unsampled IFQ, Unsampled ZMIS
ob_cs_bt %>% 
  group_by(datatype) %>% 
  summarise(n=n())
# datatype             n
# <chr>            <int>
# 1 Analysis Data  1827358
# 2 Failed Data       1028
# 3 Unsampled IFQ     1997
# 4 Unsampled ZMIS    4093

##########################
####   Unsampled NIFQ  ###
##########################

#Step 1. Identify hauls with unsampled NIFQ (following analyst manual here)
#NIFQ unsampled entries
nifq_cs <- filter(ob_cs_bt, datatype == 'Unsampled IFQ' & 
                    catch_category_code == 'NIFQ' &
                    catch_disposition == 'D' &
                    ifq == 0) 

#All hauls with unsampled NIFQ 
nifq_hauls_cs <-filter(ob_cs_bt, haul_id %in% nifq_cs$haul_id)

#Step 2. Identify hauls that need to be expanded (hauls where "Species is present in strata, but not listed as species specific discard in a given haul")

nifq_gstg_cs <- filter(nifq_hauls_cs, spid_eqv != "GSTG")
#This is all hauls with unsampled NIFQ and no GS (turns out to be all unsampled NIFQ hauls in this case)

nifq_gstg_hauls_cs <- filter(ob_cs_bt, haul_id %in% nifq_gstg_cs$haul_id)

#Step 3. Get expansion factor (sum of unsampled NIFQ discarded weight in strata)
nifq_expansion_cs <- nifq_cs %>% 
  group_by(year, state = r_state) %>% 
  summarise(total_uns_nifq = sum(mt))

#Step 4. Calculate bycatch ratio (sum of count or weight of sampled protected species/sum of sampled discarded weight for all non-IFQ species in strata), then join expansion factor and calculate expanded bycatch ratio
byc_nifq_cs <- ob_cs_bt %>%
  filter(datatype == "Analysis Data" &
           catch_disposition == "D") %>% 
  group_by(year, state = r_state) %>% 
  summarise(total_gstg_count = sum(gstg_count),
            total_nifq_mt = sum(mt[ ifq == 0 ])) %>% 
  mutate(byc_ratio = total_gstg_count/total_nifq_mt) %>% 
  left_join(nifq_expansion_cs, by = c("year", "state")) %>% 
  mutate(est_nifq_gstg = total_uns_nifq * byc_ratio,
         est_nifq_gstg = ifelse(is.na(est_nifq_gstg), 0, est_nifq_gstg)) %>%  #Because there is one year/state without any NIFQ catch
  arrange(state,year)

##########################
####   Unsampled ZMIS  ###
##########################
byc_zmis_cs <- ob_cs_bt %>%
  group_by(year, state = r_state) %>% 
  summarise(total_gstg_count  = sum(gstg_count[datatype == "Analysis Data" & #Numerator of ratio (sum of sampled discarded GS)
                                                 catch_disposition == "D"]),
            total_dis_mt = sum(mt[datatype == "Analysis Data" & #Denominator of ratio (sum of sampled discarded weight)
                                    catch_disposition == "D"]),
            total_zmis_mt = sum(mt[datatype == "Unsampled ZMIS" & #Expansion factor (sum of unsampled ZMIS discarded weight in strata)
                                     catch_category_code == "ZMIS" &
                                     species_composition_id == 0 & #I'm not really sure what this field means and the ML is not very informative...
                                     catch_disposition == "D"])) %>% 
  mutate(byc_ratio = total_gstg_count/total_dis_mt,
         est_zmis_gstg = byc_ratio * total_zmis_mt) %>% 
  arrange(state,year)

##########################
####   Unsampled UNST  ###
##########################
byc_unst_cs <- ob_cs_bt %>% 
  group_by(year, state = r_state) %>% 
  summarise(total_gstg_count  = sum(gstg_count[datatype == "Analysis Data" & #Numerator of ratio (sum of sampled discarded GS)
                                                 catch_disposition == "D"]),
            total_mt = sum(mt[datatype == "Analysis Data"]), #Denominator of ratio (sum of sampled weight for all species, retained and discarded)
            total_unst_mt = sum(mt[datatype == "Unsampled ZMIS" & #Expansion factor (sum of UNST weight in strata)
                                     catch_category_code == "UNST" &
                                     species_composition_id == 0], na.rm=T)) %>% 
  mutate(byc_ratio = total_gstg_count/total_mt,
         est_unst_gstg = byc_ratio * total_unst_mt) %>% 
  arrange(state,year)

##########################
####   Failed data     ###
##########################
failed_cs_trips <- unique(filter(ob_cs_bt, datatype == "Failed Data")$trip_id)

byc_failed_cs <- ob_cs_bt %>% 
  group_by(year, state = r_state) %>% 
  summarise(total_gstg_count  = sum(gstg_count[datatype == "Analysis Data" & #Numerator of ratio (sum of sampled discarded GS)
                                                 catch_disposition == "D"]),
            total_gfr_mt = sum(gfr_mt[datatype == "Analysis Data"], na.rm = T), #Denominator of ratio (sum of targeted retained)
            total_failed_mt = sum(gfr_mt[trip_id %in% failed_cs_trips], na.rm=T)) %>% #Expansion factor (sum of targeted retained in failed trips in strata)
  mutate(byc_ratio = total_gstg_count/total_gfr_mt,
         est_failed_gstg = byc_ratio * total_failed_mt) %>% 
  arrange(state,year)

######### Now put everything together for output ###############

#Summarise sampled data (note we treat 'Unsampled IFQ' as sampled according to YWL code)
ob_cs_samp <- ob_cs_bt %>%
  filter(datatype == "Analysis Data" | 
           (datatype == "Unsampled IFQ" & 
              catch_category_code != "NIFQ")) %>% 
  group_by(year, state = r_state) %>% 
  summarise(obs_gfr_mt = sum(gfr_mt, na.rm = T),
            obs_hauls = n_distinct(haul_id),
            obs_gstg = sum(gstg_count),
            n_obs_ves = n_distinct(drvid))

range(ob_cs_samp$n_obs_ves) #Make sure we don't have <3 vessels in any strata

ob_cs_unsamp <- ob_cs_bt %>%
  filter(datatype != "Analysis Data" & 
           !(datatype == "Unsampled IFQ" & 
               catch_category_code != "NIFQ")) %>% 
  group_by(year, state = r_state) %>% 
  summarise(unobs_gfr_mt = sum(gfr_mt, na.rm = T),
            unobs_hauls = n_distinct(haul_id))

#Put observed and estimated GSTG together, calculate coverage.
cs_out_table <- ob_cs_samp %>%
  full_join(ob_cs_unsamp, by = c("year", "state")) %>% 
  left_join(dplyr::select(byc_nifq_cs, year, state, est_nifq_gstg), by = c("year", "state")) %>% 
  left_join(dplyr::select(byc_zmis_cs, year, state, est_zmis_gstg), by = c("year", "state")) %>% 
  left_join(dplyr::select(byc_unst_cs, year, state, est_unst_gstg), by = c("year", "state")) %>%
  left_join(dplyr::select(byc_failed_cs, year, state, est_failed_gstg), by = c("year", "state")) %>%
  replace(., is.na(.), 0) %>% 
  rowwise() %>% 
  mutate(exp_gstg = sum(est_nifq_gstg, est_zmis_gstg, est_unst_gstg),
         total_exp_gstg = sum(obs_gstg, est_nifq_gstg, est_zmis_gstg, est_unst_gstg, na.rm = T),
         total_gfr_mt = sum(obs_gfr_mt, unobs_gfr_mt, na.rm = T),
         pct_gfr_obs = (obs_gfr_mt / (obs_gfr_mt + unobs_gfr_mt)) * 100) %>% 
  dplyr::select(state, year, obs_gstg, obs_gfr_mt, total_gfr_mt, pct_gfr_obs, exp_gstg, total_exp_gstg) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  arrange(desc(state),year)

names(cs_out_table) <- c("State", "Year", "Observed bycatch", "Observed groundfish landings (MT)", "Fleet-total groundfish landings (MT)", "Groundfish landings sampled (%)", "Estimated bycatch from unsampled catch", "Fleet-total bycatch")

##%######################################################%##
#                                                          #
####               Spatial data for maps                ####
#                                                          #
##%######################################################%##

#We need lat, lon, target species weight, GSTG count, and vessel id (for confidentiality). This gets sent to KS for mapping magic.

######LE/CS########
ob_le <- ob %>% 
  filter(datatype == "Analysis Data" & 
           sector == "Limited Entry Trawl") #Note that following YWL methods we don't limit this to bottom trawl gear

ob_le_formap <- ob_le %>%
  filter(gear == "Bottom Trawl") %>% 
  dplyr::select(sector, drvid, haul_id, avg_lat, avg_long, gfr_mt, gstg_count)

#Subset to CS bottom trawl data (only non-hake CS sector that encounters GSTG)
ob_cs_bt <- ob_cs %>% 
  filter(cs_sector == "CS - Bottom Trawl")

ob_cs_formap <- ob_cs_bt %>% 
  dplyr::select(sector, drvid, haul_id, avg_lat, avg_long, gfr_mt, gstg_count)

ob_lecs_formap <- bind_rows(ob_le_formap)

write.csv(ob_lecs_formap, file = paste0(out_drive, "LE_CS_spatial_GSTG_data_", tday, ".csv"))
####CHLB#########
#Observed CHLB data
ob_chlb <- ob %>% 
  filter(sector == "OA CA Halibut" |
           sector == "LE CA Halibut") %>% 
  mutate(stratum = paste(sector,year,season,sep="_"))

ob_chlb_formap <- ob_chlb %>%
  dplyr::select(sector, drvid, haul_id, avg_lat, avg_long, chlb_mt, gstg_count)

write.csv(ob_chlb_formap, file = paste0(out_drive, "CHLB_spatial_GSTG_data_", tday, ".csv"))


  