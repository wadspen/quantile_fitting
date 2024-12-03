library(dplyr)
teams <- list.files("../model-fits")
print(teams)
all_scores <- data.frame()
for (i in 1:length(teams)) {
	file_name <- paste0("../model-fits/", teams[i], "/scores/all_scores.rds")
	scores <- readRDS(file_name)
	print(dim(scores))
 	all_scores <- rbind(all_scores, scores)
 }


all_scores <- all_scores %>%
	filter(location == "19")

saveRDS(all_scores, "iowa_scores.rds")

