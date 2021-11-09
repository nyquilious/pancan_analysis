# read in the cleaned data
covid_data = read_tsv("data/clean/covid_data.tsv")

# split into train and test (set seed here if applicable)
test_states = c("Alabama", "Arizona", "Arkansas", 
                "California", "Colorado", "Connecticut")
covid_train = covid_data %>% filter(!(state %in% test_states))
covid_test = covid_data %>% filter(state %in% test_states)

# save the train and test data
write_tsv(x = covid_train, file = "data/clean/covid_train.tsv")
write_tsv(x = covid_test, file = "data/clean/covid_test.tsv")