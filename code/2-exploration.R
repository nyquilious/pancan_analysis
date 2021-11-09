# load libraries
library(kableExtra)                     # for printing tables
library(cowplot)                        # for side by side plots
library(lubridate)                      # for dealing with dates
library(maps)                           # for creating maps
library(tidyverse)

# read in the cleaned data
covid_data = read_tsv("data/clean/covid_data.tsv")

# calculate median case fatality rate
median_case_fatality_rate = covid_data %>%
  summarise(median(case_fatality_rate)) %>%
  pull()

# create histogram of case fatality rate
p = covid_data %>%
  ggplot(aes(x = case_fatality_rate)) + 
  geom_histogram() +
  geom_vline(xintercept = median_case_fatality_rate,
             linetype = "dashed") +
  labs(x = "Case fatality rate (percent)", 
       y = "Number of counties") +
  theme_bw()

# save the histogram
ggsave(filename = "results/response-histogram.png", 
       plot = p, 
       device = "png", 
       width = 5, 
       height = 3)

# examine top 10 counties by case fatality rate
covid_data %>% 
  select(county, state, case_fatality_rate) %>%
  arrange(desc(case_fatality_rate)) %>%
  head(10) %>%
  write_tsv("results/top-10-counties-data.tsv")

# create a heatmap of case fatality rate across the U.S.
p = map_data("county") %>%
  as_tibble() %>% 
  left_join(case_data %>% 
              rename(region = state, 
                     subregion = county,
                     `Case Fatality Rate` = case_fatality_rate) %>% 
              mutate(region = str_to_lower(region), 
                     subregion = str_to_lower(subregion)), 
            by = c("region", "subregion")) %>%
  ggplot() + 
  geom_polygon(data=map_data("state"), 
               aes(x=long, y=lat, group=group),
               color="black", fill=NA,  size = 1, alpha = .3) + 
  geom_polygon(aes(x=long, y=lat, group=group, fill = `Case Fatality Rate`),
               color="darkblue", size = .1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_void()

ggsave(filename = "results/response-map.png", 
       plot = p, 
       device = "png", 
       width = 7, 
       height = 4)