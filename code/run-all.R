# run all steps of the analysis pipeline
source("code/0-download.R")
source("code/1-cleaning.R")
source("code/2-exploration.R")
source("code/3-train-test-split.R")
source("code/4-regression-modeling.R")
source("code/5-tree-modeling.R")
source("code/6-model-evaluation.R")