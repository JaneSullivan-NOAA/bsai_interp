# Queries for sablefish and turbot longline survey RPNs and RPWs
# Contacts jane.sullivan@noaa.gov
# Updated January 2024
# R version 4.2.3

# set up ----
libs <- c("tidyverse", "RODBC", "zoo")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

theme_set(theme_minimal(base_size = 13))

username_akfin <- '' # fill in your akfin creds
password_akfin <- ''

channel_akfin <- odbcConnect("akfin", uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

# RPN queries ----

# BSAI turbot RPNs (BS = odd years, AI = even years). Full time series starts in 1997.
# The raw 1996 AI in area 15 (NE AI Slope) was lost, though we have area estimates. 
turbot <- sqlQuery(channel_akfin, query = paste0("
                select    species, year, council_management_area as strata, rpn, rpn_var, rpw, rpw_var
                from      afsc.lls_fmp_subarea_all_strata
                where     country = 'United States' and
                          species_code in ('10115') and
                          ((council_management_area in 'Bering Sea' and mod(year, 2) <> 0) or
                          (council_management_area in ('Aleutians') and mod(year, 2) = 0)) and
                          year >= 1996 
                order by  year asc
                ")) %>% 
  rename_all(tolower) %>% 
  arrange(species, strata, year)
write_csv(turbot, 'data/turbotarearpns.csv')

# sablefish RPNs for calculating new AK-wide indices
sable <- sqlQuery(channel_akfin, query = paste0("
                select    species, year, council_management_area as strata, rpn, rpn_var, rpw, rpw_var
                from      afsc.lls_fmp_subarea_3_to_7_depred
                where     country = 'United States' and
                          species_code in ('20510') and
                          ((council_management_area in 'Bering Sea' and mod(year, 2) <> 0) or
                          (council_management_area in 'Aleutians' and mod(year, 2) = 0) or 
                          council_management_area in ('Western Gulf of Alaska', 'Eastern Gulf of Alaska', 'Central Gulf of Alaska')) and
                          year >= 1995 
                order by  year asc
                ")) %>% 
  rename_all(tolower) %>% 
  arrange(species, strata, year)
write_csv(sable, 'data/sablearearpns.csv')

# current akwide sablefish
aksable <- sqlQuery(channel_akfin, query = paste0("
                select    species, year, rpn, rpn_var, rpw, rpw_var
                from      afsc.lls_ak_wide_3_to_7_depred
                where     country = 'United States' and
                          species_code in ('20510') and
                          year >= 1996 
                order by  year asc
                ")) %>% 
  rename_all(tolower) 
write_csv(aksable, 'data/sableakwide.csv')


