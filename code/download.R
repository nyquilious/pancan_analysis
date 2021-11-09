# load libraries
library(tidyverse)

# download case data from NY Times website
url = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"
case_data_raw = read_csv(url)

# write raw data to file
write_tsv(x = case_data_raw, file = "data/raw/case_data_raw.tsv")

# download health rankings data from the web 
# (omitted from this template)