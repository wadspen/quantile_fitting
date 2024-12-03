library(stringr)
library(dplyr)

args <- commandArgs()
dist <- args[6]

loc <- paste0("./sim_coverage/", dist)
files <- list.files(loc)


for (i in 1:length(files)) {
	nums <- str_extract_all(files[i], "\\(?[0-9]+\\)?")[[1]]
	n <- nums[1]
	quants <- nums[2]

	file <- readRDS(paste0(loc, "/", files[i]))
	file <- file %>%
		mutate(n = n, quants = quants)

	saveRDS(file, paste0(loc, "/", files[i]))

}








