# load libraries
library(glmnetUtils)                    # to run ridge and lasso
source("code/plot_glmnet.R")            # for lasso/ridge trace plots

# read in the training data
covid_train = read_tsv("data/clean/covid_train.tsv")

# run ridge regression
set.seed(1)
ridge_fit = cv.glmnet(case_fatality_rate ~ . - state - county - fips,   
                      alpha = 0,                 
                      nfolds = 10,               
                      data = covid_train)

# save the ridge fit object
save(ridge_fit, file = "results/ridge_fit.Rda")

# run lasso regression
set.seed(1)
lasso_fit = cv.glmnet(case_fatality_rate ~ . - state - county - fips,   
                      alpha = 1,                 
                      nfolds = 10,               
                      data = covid_train)

# save the lasso fit object
save(lasso_fit, file = "results/lasso_fit.Rda")

# create lasso CV plot
png(width = 6, 
    height = 4,
    res = 300,
    units = "in", 
    filename = "results/lasso-cv-plot.png")
plot(lasso_fit)
dev.off()

# create lasso trace plot
p = plot_glmnet(lasso_fit, covid_train, features_to_plot = 6)
ggsave(filename = "results/lasso-trace-plot.png", 
       plot = p, 
       device = "png", 
       width = 6, 
       height = 4)

# extract features selected by lasso and their coefficients
beta_hat_std = extract_std_coefs(lasso_fit, covid_train)
beta_hat_std %>%
  filter(coefficient != 0) %>%
  arrange(desc(abs(coefficient))) %>% 
  write_tsv("results/lasso-features-table.tsv")
