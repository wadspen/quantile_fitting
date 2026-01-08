library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)



files <- list.files(pattern = ".rds")
files <- files[!str_detect(files, "comb_res")]


all_res <- data.frame()
for (i in 1:length(files)) {
	res <- readRDS(files[i])
	all_res <- rbind(all_res, res)	
}

saveRDS(all_res, "comb_res.rds")
