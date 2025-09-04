library(dplyr)
library(stringr)
library(tidyr)

models <- list.files("../model-fits-ord/")

all_scores <- data.frame()
for (i in 1:length(models)) {
	file_name <- paste0("../model-fits-ord/", models[i],
			    "/scores/all_scores.rds")
	if (dir.exists(paste0("../model-fits-ord/", models[i], "/scores")) == FALSE) {next}
	#print(dir.exists(paste0("../model-fits/", models[i], "/scores")))
	scoreRDS <- readRDS(file_name)
	all_scores <- rbind(all_scores, scoreRDS)
}

saveRDS(all_scores, "all_ord_scores.rds")

















