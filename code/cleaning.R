# load libraries
library(lubridate)
library(tidyverse)

# load raw case data
case_data_raw = read_tsv(file = "data/raw/case_data_raw.tsv")

# clean case data
case_data = case_data_raw %>%
  na.omit() %>%                               # remove NA values
  filter(year(date) == 2020) %>%              # keep data from 2020 
  group_by(fips, county, state) %>%           # group by county
  summarise(total_cases = sum(cases),         # total cases per county
            total_deaths = sum(deaths)) %>%   # total deaths per county
  ungroup() %>%
  mutate(case_fatality_rate =                 # case_fatality_rate = 
           total_deaths/total_cases*100) %>%  #  total_deaths/total_cases
  select(-total_cases, -total_deaths)         # remove intermediate variables

# load raw county health data
# (omitted from this template)

# clean county health data
# (omitted from this template, reading from file instead)
county_health_data = read_tsv("data/raw/county_health_data.tsv")

# join county health data with case data
covid_data = inner_join(county_health_data, case_data, by = "fips")

# write cleaned data to file
write_tsv(covid_data, file = "data/clean/covid_data.tsv")