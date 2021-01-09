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
library(devtools) #to source functions off my github
#library(broom)

#For environmental data/gams:
library(raster) #Extracting gridded data at given locations
library(ncdf4) # .nc SST files
library(mgcv) #GAMs
library(mgcViz) #Plotting GAMs. Note that this package isn't compatible with version of R currently on Tantalus...
library(lubridate) #Julian date etc
library(marmap) #Bathymetry

source("~/observer/Input/load_data_2020-06-24.R") #

source_url("https://raw.githubusercontent.com/kricherson-NOAA/ratio-bycatch-expansion/master/do_ratio_multi.R") #does bycatch ratios
source_url("https://raw.githubusercontent.com/kricherson-NOAA/ratio-bycatch-expansion/master/do_ratio_est_multi.R") #does ratio expansion
source_url("https://raw.githubusercontent.com/kricherson-NOAA/ratio-bycatch-expansion/master/do_boot_multi.R") #does bootstrapped expansions
source_url("https://raw.githubusercontent.com/kricherson-NOAA/ratio-bycatch-expansion/master/do_cs_expansions.R") #does CS expansions

load_data(c("WCGOP","ASHOP_proc","FT_proc", "EM", "BIO")) #load_data()

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


#CS data needs to be processed to get the unsampled data. However, may only want to do this once and save the df in order to save time. 
process_cs_pre <- FALSE #set to TRUE if want to process data and save it. Set to FALSE if want to load the df of processed data

if(process_cs_pre)
#Pull out Catch Shares OB data. Starting with OBOrig_Pre and processing it to get the unsampled data. Following KS's methods from GM_2017 and YWL's old GSTG code. 
{
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
  
  save(ob_cs, file = "ob_cs.Rdata")
  
  }else{
    load("ob_cs.Rdata")
}
#
#

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

#keep track of sectors that are done?
done_sectors <- character(0)

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

#add these sectors to the done sectors list
done_sectors <- c(done_sectors, unique(ncs_gstg_out$sector))

##%######################################################%##
#                                                          #
####                 Catch shares trawl                 ####
#                                                          #
##%######################################################%##
#Subset to CS bottom trawl data (only non-hake CS sector that encounters GSTG)
ob_cs_bt <- ob_cs %>% 
  filter(cs_sector == "CS - Bottom Trawl")

###############Expand to unsampled data####################

cs_gstg <- do_cs_expansion(ob_cs_bt, c("year", "r_state"), "Green Sturgeon")

range(cs_gstg$n_obs_ves)
#[1]  3 46

cs_gstg_out <- cs_gstg %>% 
  mutate_if(is.numeric, round, 1) %>%
  dplyr::select(Year = year, 
                State= r_state, 
                `Observed bycatch` = obs_byc_ct,
                `Observed groundfish landings (MT)` = obs_tgt_mt,
                `Fleet-total groundfish landings (MT)` = total_tgt_mt,
                `Groundfish landings sampled (%)` = pct_tgt_obs,
                `Estimated bycatch from unsampled catch` = est_unsamp_ct,
                `Fleet total bycatch` = total_byc_ct) %>% 
  arrange(State, Year)

write_csv(cs_gstg_out, paste0(out_drive, "cs_gstg_", tday, ".csv"))

done_sectors <- c(done_sectors, paste(unique(ob_cs_bt$sector), unique(ob_cs_bt$gear_2)))


##%######################################################%##
#                                                          #
####                 California Halibut                 ####
#                                                          #
##%######################################################%##

#Relevant notes from YWL code:
## NOTE 2: There are 3 final products based on stratification
#  (1) 1st product is summarized with season strata, 
#  (2) 2nd product is summarized without season strata (but only yearly) for the years 2011-2013.
#  (3) 3rd product is summarized based on combined sectors (LE+OA) at annual time step for the years 2011-2013.
#  Prior to 2011, the result should be presented with season (as done in previous report) 
#  Since 2011, the result should be presented at annual time scale to avoid confidentiality issue

#  Also, bycatch numbers should be reported as LE+OA combined, because of too few vessel numbers operating in LE sector since 2011
#  For LE portion of CA halibut fishery, landing amount cannot be seprately reported in the annual GSTG bycatch report to avoid confidentiality issue, 
#  because coverage table product (by KS) does not report LE coverage as separate, but as combined with IFQ bottom trawl fishery sector

# Note that, because of low fishing rates in 'LE CA halibut' fishery after 2010, LE and OA need to be combined for the later years for bycatch numbers #
# Plus, coverage rate calculation for 'LE CA halibut' is close to 100%, but 'Unsampled ZMIS' portion needs to be accessed separately.
#End YWL notes

#Observed CHLB data
ob_chlb <- ob %>% 
  filter(sector == "OA CA Halibut" |
           sector == "LE CA Halibut") %>% 
  mutate(stratum = paste(sector,year,season,sep="_"))

#Fish ticket CHLB data
ft_chlb <- filter(ft,sector == "OA CA Halibut" |
                  sector == "LE CA Halibut") %>% 
  mutate(stratum = paste(sector,year,season,sep="_"))

########################################################
#First: OA and LE separate, 2002-2010, winter and summer
########################################################

#This is an easy way to get coverage so we know which strata to bootstrap, though we won't actually use the estimates
chlb_gstg1 <- do_ratio_multi(ob_dat = filter(ob_chlb, year < 2011),
                                 #ft_dat = filter(ft_chlb, year < 2011),
                                 strata = c("year", "sector", "season","stratum"),
                                 expfactor = "tgt_mt",
                                 bycatchspp = "Green Sturgeon",
                                 bycatchunit = "gstg_count",
                                 management_groups = FALSE)

#as.data.frame(dplyr::select(chlb_gstg1, year, sector, season, n_obs_ves, stratum)) #To guide choice of pooling strata

#This identifies confidential strata and the adjacent strata we will pool over for bootstrapping. Easier to ID by hand because there are so many special cases.
chlb_conf_strata <- filter(chlb_gstg1, n_obs_ves < 3 & !is.na(n_obs_ves) & year < 2011) %>% 
  dplyr::select(year, season, stratum, sector, n_obs_ves) %>% 
  mutate(adjacent_stratum1 = case_when(stratum == "LE CA Halibut_2002_summer" ~ "LE CA Halibut_2003_summer", 
                                       stratum == "LE CA Halibut_2010_summer" ~ "LE CA Halibut_2009_summer",
                                       stratum == "LE CA Halibut_2010_winter" ~ "LE CA Halibut_2008_winter",
                                       stratum == "OA CA Halibut_2003_winter" ~ "OA CA Halibut_2004_winter",
                                       stratum == "OA CA Halibut_2004_summer" ~ "OA CA Halibut_2003_summer",
                                       stratum == "OA CA Halibut_2005_summer" ~ "OA CA Halibut_2004_summer",
                                       stratum == "OA CA Halibut_2009_summer" ~ "OA CA Halibut_2008_summer",
                                       stratum == "OA CA Halibut_2009_winter" ~ "OA CA Halibut_2008_winter"
  ), 
  adjacent_stratum2 = case_when(stratum == "LE CA Halibut_2002_summer" ~ "LE CA Halibut_2004_summer",
                                stratum == "LE CA Halibut_2010_summer" ~ "LE CA Halibut_2011_summer",
                                stratum == "LE CA Halibut_2010_winter" ~ "LE CA Halibut_2011_winter",
                                stratum == "OA CA Halibut_2003_winter" ~ "OA CA Halibut_2005_winter",
                                stratum == "OA CA Halibut_2004_summer" ~ "OA CA Halibut_2005_summer",
                                stratum == "OA CA Halibut_2005_summer" ~ "OA CA Halibut_2007_summer",
                                stratum == "OA CA Halibut_2009_summer" ~ "OA CA Halibut_2010_summer",
                                stratum == "OA CA Halibut_2009_winter" ~ "OA CA Halibut_2010_winter"
  )) 

#unobs_strata <- filter(chlb_gstg1, is.na(n_obs_ves))$stratum

#For bootstrapping, create a df that includes the pooled data. Do this by first getting the unpooled data, then looping through the strata that need to be pooled and assigning the adjacted strata to the pooled stratum
ob_chlb_pooled <- ob_chlb %>% 
  filter(year < 2011 &
           !(stratum %in% chlb_conf_strata$stratum) & 
           datatype =="Analysis Data") %>% 
  mutate(pooled = "no") %>%  #For checking purposes...
  ungroup()

#Loop through confidential strata, summarise by vessel across pooled strata, add to df defined above
for(i in 1: nrow(chlb_conf_strata))
{
  pooled_ob_df <- ob_chlb %>% 
    filter(stratum %in% chlb_conf_strata[i, c("stratum", "adjacent_stratum1", "adjacent_stratum2")] & 
             datatype =="Analysis Data") %>% 
    mutate(stratum = chlb_conf_strata$stratum[i]) %>% #stratum now contains adjacent stratum.
    mutate(pooled = "yes")
  
    ob_chlb_pooled <- bind_rows(ob_chlb_pooled, pooled_ob_df)
}

#make sure we have all strata (should be TRUE)
nrow(distinct(dplyr::select(ungroup(ob_chlb_pooled),stratum,pooled)))==length(unique(filter(ob_chlb, year<2011 & datatype == "Analysis Data")$stratum))
unique(sort(distinct(dplyr::select(ungroup(ob_chlb_pooled),stratum,pooled))$stratum) == sort(unique(filter(ob_chlb, year<2011 & datatype == "Analysis Data")$stratum)))

#make sure all "strata" now have >=3 vessels 
ob_chlb_pooled %>% 
  group_by(stratum) %>% 
  summarise(n_vessels = n_distinct(drvid)) %>% 
  ungroup %>% 
  summarise(min=min(n_vessels),max=max(n_vessels))

# A tibble: 1 x 2
# min   max
# <dbl> <dbl>
#   1     3    13

#Now should be able to calculate bycatch ratios and estimated bycatch using the do_boot_multi function, though note that coverage output does not account for pooling (I think that's how we've reported it in the past, so that is not a problem). Note that the strata that we report at and the strata that we use to summarise data at the vessel level are different ("strata" and "vessel_strata" in the function, respectively)

chlb_booted <- do_boot_multi(ob_dat = ob_chlb_pooled, 
                                    ft_dat = filter(ft_chlb, year < 2011), 
                                    strata = "stratum", 
                                    vessel_strata = c("stratum", "sector", "year", "season"), 
                                    expfactor = "tgt_mt", 
                                    bycatchspp = "Green Sturgeon", 
                                    bycatchunit = "gstg_count", 
                                    seed = 42)


chlb_02_to_10_out <- chlb_booted %>% 
  separate(stratum, c("Sector", "Year", "Season"), sep = "_") %>% #get sector, year, season from the stratum
  mutate_at(vars(point_est_ratio, upper_ci, lower_ci), round, 2) %>% 
  mutate_at(vars(total_expf, fleet_expf, pct_cvg), round, 1) %>% 
  mutate_at(vars(est_byc, est_byc_lower_trunc, est_byc_upper), round, 0) %>% 
  mutate(pct_cvg = ifelse(is.na(n_obs_ves), 0, pct_cvg)) %>% 
  mutate_at(vars(-n_obs_ves), as.character) %>% #So we can replace confidential data with * and unobserved years with --
  mutate_at(vars(total_byc, total_expf, pct_cvg), list(~ifelse( n_obs_ves < 3, "*", .))) %>% 
  mutate_at(vars(-fleet_expf, -n_obs_ves, -Sector, -Year, -Season, -pct_cvg), list(~ifelse(is.na(n_obs_ves), "--", .))) %>% 
  dplyr::select(-n_obs_ves) %>%
  ungroup() %>% 
  dplyr::select(Sector, 
                Year, 
                Season, 
                `Observed bycatch` = total_byc, 
                `Observed CA halibut landings (MT)` = total_expf, 
                `Fleet total CA halibut landings (MT)` = fleet_expf,
                `CA halibut landings sampled (%)` = pct_cvg, 
                `Bycatch ratio` = point_est_ratio, 
                `Lower CI of ratio` = lower_ci,
                `Upper CI of ratio` = upper_ci,
                `Fleet-total bycatch` = est_byc,
                `Lower CI of bycatch` = est_byc_lower_trunc, 
                `Upper CI of bycatch` = est_byc_upper) %>% 
  arrange(Sector, Year, desc(Season))

########################################################
#Second: OA and LE combined, 2011-2017
########################################################
#Following last year: estimate OA and LE separately, then combine
oa_chlb_booted <- do_boot_multi(ob_dat = filter(ob_chlb, sector == "OA CA Halibut" & year >=2011), 
                             ft_dat = filter(ft_chlb, sector == "OA CA Halibut" & year >= 2011), 
                             strata = "year", 
                             vessel_strata = "year", 
                             expfactor = "tgt_mt", 
                             bycatchspp = "Green Sturgeon", 
                             bycatchunit = "gstg_count", 
                             seed = 42)

le_chlb <- do_cs_expansion(filter(ob_chlb, sector == "LE CA Halibut" & year >=2011), 
                           strata = "year", 
                           bycatchspp = "Green Sturgeon") #almost no expanded GSTG (will be lost in rounding)

#combine LE and OA
chlb_11_to_19_out <- oa_chlb_booted %>% 
  mutate(year = as.numeric(year)) %>% 
  left_join(dplyr::select(le_chlb, year, obs_byc_ct, total_byc_ct, pct_obs_le = pct_tgt_obs, le_obs_tgt_mt = obs_tgt_mt),
            by = "year") %>% 
  replace(is.na(.), 0) %>% 
  mutate(obs_gstg_count_oale = obs_byc_ct + total_byc, #Sum of observed LE and OA GSTG
         exp_gstg_count_oale = total_byc_ct + est_byc, #Sum of expanded LE and OA GSTG
         lower_ci_bycatch_oale = total_byc_ct + est_byc_lower_trunc,
         upper_ci_bycatch_oale = total_byc_ct + est_byc_upper) %>% 
  mutate(lower_ci_bycatch_oale = ifelse(lower_ci_bycatch_oale < obs_gstg_count_oale, #Truncate
                                        obs_gstg_count_oale, 
                                        lower_ci_bycatch_oale),
         pct_obs_le = as.character(pct_obs_le),
         pct_obs_le = ifelse(pct_obs_le == "0", "--", pct_obs_le)) %>%
  mutate_if(is.numeric, round, 1) %>%
  mutate(year = as.character(year)) %>% 
  mutate(year = ifelse(pct_obs_le == "--",
                       paste0(year, "*"), 
                       year)) %>% #To indicate no LE landings
  dplyr::select(`Year` = year, 
                `LE + OA combined observed bycatch` = obs_gstg_count_oale, 
                `LE + OA combined fleet-total bycatch` = exp_gstg_count_oale, 
                `Lower CI of fleet-total bycatch` = lower_ci_bycatch_oale, 
                `Upper CI of fleet-total bycatch` = upper_ci_bycatch_oale, 
                `LE CA halibut landings sampled (%)` = pct_obs_le, 
                `OA observed CA halibut landings (MT)` = total_expf, 
                `OA fleet-total CA halibut landings (MT)` = fleet_expf, 
                `OA CA halibut landings sampled (%)` = pct_cvg)

##################### write both to csv chlb_02_to_10_out

write_csv(chlb_02_to_10_out, paste0(out_drive, "chlb_02_to_10_gstg_out_", tday, ".csv"))

write_csv(chlb_11_to_19_out, paste0(out_drive, "chlb_11_to_19_gstg_out_", tday, ".csv"))

##%######################################################%##
#                                                          #
####                    At-Sea Hake                     ####
#                                                          #
##%######################################################%##

ashop_gstg_out <- ashop %>%
  mutate(sector2 = case_when(sector == "ATSEA TRIBAL" ~ "Tribal Mothership",
                             sector == "MOTHERSHIP" ~ "Non-tribal Mothership",
                             sector == "CATCHER-PROCESSOR" ~ "Catcher Processor")) %>% 
  mutate(exp_gstg_count = ifelse(spid_eqv == "GSTG" & !is.na(spid_eqv), 
                                 expanded_num_2sector_level, 0)) %>% 
  group_by(sector2,year) %>% 
  summarise(total_gstg_count = sum(gstg_count),
            total_exp_gstg_count = sum(exp_gstg_count),
            n_samp_tows = n_distinct(haul_id[haul_sampled_by !=0 ]),
            n_tows = n_distinct(haul_id),
            total_pwht_mt = sum(pwht_mt),
            pct_sampled = (n_samp_tows/n_tows) * 100,
            n_obs_ves = n_distinct(drvid))%>% 
  ungroup() %>% 
  mutate_if(is.numeric, round, digits = 1) %>% 
  mutate_at(c("total_gstg_count", "n_samp_tows","total_pwht_mt", "pct_sampled"), as.character) %>% #So that we can replace confidential data with *
  mutate_at(c("total_gstg_count", "n_samp_tows","total_pwht_mt", "pct_sampled"), funs(ifelse( n_obs_ves < 3, "*", .))) %>% #Replace confidential data with *
  tidyr::complete(year, sector2) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(-n_tows, -n_obs_ves) %>% 
  filter(sector2 != "Tribal Mothership") %>% 
  arrange(sector2, year) %>% 
  rename(Year = year,
         Sector = sector2,
         `Observed bycatch` = total_gstg_count,
         `Fleetwide expanded bycatch` = total_exp_gstg_count,
         `Number sampled tows` = n_samp_tows,
         `Sampled hake landings (MT)` = total_pwht_mt,
         `% tows sampled` = pct_sampled)

write_csv(ashop_gstg_out, paste0(out_drive, "ashop_gstg_out.csv", tday, ".csv"))


############################################################################################

################## End of bycatch tables ###################################################

############################################################################################


##%######################################################%##
#                                                          #
####                Depth distributions                 ####
#                                                          #
##%######################################################%##

#It looks like we do depth by haul for bottom trawl fisheries in previous reports

#Cutting into depth bins takes a while so maybe just do it once, then save and reload

make_depth_bin_dfs <- TRUE

if (make_depth_bin_dfs)
{
  #OB/LE 
  ob_le_cs <- filter(ob, (sector == "Catch Shares" & gear %in% c("Bottom Trawl"))|
                       sector == "Limited Entry Trawl") %>% 
    filter(!is.na(avg_depth)) %>%  # 2 entries
    group_by(haul_id, avg_depth, sector) %>% 
    summarise(total_gstg = sum(gstg_count),
              drvid=drvid[1]) %>% 
    mutate(depth_bin1 = cut(avg_depth, breaks=c(seq(0, 1200, by = 100), Inf), 
                            labels=c(paste(seq(0, 1100,by= 100),seq(100, 1200, by= 100),sep="-"), ">1200")), #All depths in 100 fm bins
           depth_bin2 = cut(avg_depth, breaks=c(seq(0, 100, by = 10), Inf), 
                            labels=c(paste(seq(0, 90,by = 10),seq(10, 100 ,by = 10),sep="-"),">100"))) %>% #Depths <100 in 10 fm bins
    ungroup()
  
  save(ob_le_cs, file = "ob_le_cs.Rdata")
  
  ob_chlb <- ob %>% 
    filter(sector == "OA CA Halibut" |
             sector == "LE CA Halibut")  %>% 
    filter(!is.na(avg_depth)) %>%  # 1 entry
    group_by(haul_id, avg_depth, sector) %>% 
    summarise(total_gstg = sum(gstg_count),
              drvid = drvid[1]) %>% #Could be some bins with <3 vessels? %>% 
    mutate(depth_bin1 = cut(avg_depth, breaks=c(seq(0, 1200, by = 100), Inf), 
                            labels=c(paste(seq(0, 1100,by= 100),seq(100, 1200, by= 100),sep="-"), ">1200")), #All depths in 100 fm bins
           depth_bin2 = cut(avg_depth, breaks=c(seq(0, 100, by = 10), Inf), 
                            labels=c(paste(seq(0, 90,by = 10),seq(10, 100 ,by = 10),sep="-"),">100"))) %>% 
    ungroup()
  
  save(ob_chlb, file = "ob_chlb.Rdata")
  
} else{
  
  load("ob_le_cs.Rdata")
  load("ob_chlb.Rdata")
  
  
}

##############Make LE/CS figure

#check confidentiality (in the bins I'll be putting in the report)
ob_le_cs %>% 
  group_by(depth_bin2, sector) %>% 
  summarise(n_ves=n_distinct(drvid)) %>% 
  filter(n_ves<3)

#I am including this figure in the report 
le_cs_bins <- ob_le_cs %>% 
  mutate(sector = ifelse(sector == "Catch Shares", "Individual Fishing Quota Trawl", sector)) %>% #For consistency (use IFQ not CS in report)
  filter(!is.na(avg_depth) & depth_bin2 != ">100") %>% 
  group_by(sector) %>% 
  mutate(n_hauls = n_distinct(haul_id)) %>% 
  ungroup() %>% 
  mutate(sector_n = paste0(sector,", n=", n_hauls)) %>% #This was done so the facet titles will have the total number of observed hauls in them
  group_by(depth_bin2, sector_n) %>% 
  summarise(n_hauls = n_distinct(haul_id),
            pct_hauls_gstg = (n_distinct(haul_id[total_gstg>0]) / n_distinct(haul_id)) * 100) %>% 
  ggplot()+
  aes(x = depth_bin2, y = n_hauls, fill = pct_hauls_gstg)+
  geom_bar(stat="identity")+
  facet_wrap(~sector_n)+
  scale_fill_viridis(name="% observed\nhauls with\ngreen sturgeon")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Depth (fm)")+
  ylab("Number of hauls")

ggsave(le_cs_bins, file = paste0(out_drive, "le_cs_depthbins_", tday, ".png"), height =5, width=8)

####### Make CHLB figure
ob_chlb %>% 
  group_by(depth_bin2, sector) %>% 
  summarise(n_ves=n_distinct(drvid)) %>% 
  filter(n_ves<3)

chlb_bins <- ob_chlb %>% 
  filter(!is.na(avg_depth) & depth_bin2 != ">100") %>% 
  group_by(sector) %>% 
  mutate(n_hauls = n_distinct(haul_id)) %>% 
  ungroup() %>% 
  mutate(sector_n = paste0(sector,", n=", n_hauls)) %>% #This was done so the facet titles will have the total number of observed hauls in them
  group_by(depth_bin2, sector_n) %>% 
  summarise(n_hauls = n_distinct(haul_id),
            pct_hauls_gstg = (n_distinct(haul_id[total_gstg>0]) / n_distinct(haul_id))*100,
            n_ves = n_distinct(drvid)) %>% 
  #filter(n_ves>2) %>% 
  ggplot()+
  aes(x = depth_bin2, y = n_hauls, fill = pct_hauls_gstg)+
  geom_bar(stat="identity")+
  facet_wrap(~sector_n)+
  scale_fill_viridis(name="% observed\nhauls with\ngreen sturgeon")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Depth (fm)")+
  ylab("Number of hauls")

ggsave(chlb_bins, file = paste0(out_drive, "chlb_depthbins_", tday, ".png"), height =5, width=8)




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
  group_by(sector, drvid, haul_id, avg_lat, avg_long, haul_duration) %>% 
  summarise(gstg_count = sum(gstg_count),
            gfr_mt = sum(gfr_mt, na.rm = T))

#Subset to CS bottom trawl data (only non-hake CS sector that encounters GSTG)
ob_cs_bt <- ob_cs %>% 
  filter(datatype == "Analysis Data" & 
           cs_sector == "CS - Bottom Trawl") 

ob_cs_formap <- ob_cs_bt %>% 
  group_by(sector, drvid, haul_id, avg_lat, avg_long, haul_duration) %>% 
  summarise(gstg_count = sum(gstg_count),
            gfr_mt = sum(gfr_mt, na.rm = T))

ob_lecs_formap <- bind_rows(ob_le_formap, ob_cs_formap) %>% 
  mutate(cpue_haul_dur = gstg_count/haul_duration,
         cpue_tgt_mt = gstg_count/gfr_mt)

write.csv(ob_lecs_formap, file = paste0(out_drive, "LE_CS_spatial_GSTG_data_", tday, ".csv"))
####CHLB#########
#Observed CHLB data
ob_chlb <- ob %>% 
  filter((sector == "OA CA Halibut" |
           sector == "LE CA Halibut") & 
           datatype == "Analysis Data") 

ob_chlb_formap <- ob_chlb %>%
  group_by(sector, drvid, haul_id, avg_lat, avg_long, haul_duration) %>% 
  summarise(gstg_count = sum(gstg_count),
            chlb_mt = sum(chlb_mt, na.rm = T)) %>% 
  mutate(cpue_haul_dur = gstg_count/haul_duration,
         cpue_tgt_mt = gstg_count/chlb_mt)

write.csv(ob_chlb_formap, file = paste0(out_drive, "CHLB_spatial_GSTG_data_", tday, ".csv"))


  