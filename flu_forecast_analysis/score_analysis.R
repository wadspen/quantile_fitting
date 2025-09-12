library(dplyr)
library(stringr)
library(tidyr)

fits_loc <- "model-fits-ord/"

models <- list.files(paste0("../", fits_loc))

all_scores <- data.frame()
for (i in 1:length(models)) {
	file_name <- paste0("../", fits_loc, models[i],
			    "/scores/all_scores.rds")
	if (dir.exists(paste0("../", fits_loc, models[i], "/scores")) == FALSE) {next}
	#print(dir.exists(paste0("../model-fits/", models[i], "/scores")))
	scoreRDS <- readRDS(file_name)
	if (sum(str_detect(colnames(scoreRDS), "cover")) == 0) {next}
	all_scores <- rbind(all_scores, scoreRDS)
}

saveRDS(all_scores, "all_scores_ord.rds")

















