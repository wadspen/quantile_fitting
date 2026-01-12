

library(stringr)
files <- list.files(pattern = ".rds")

files <- files[!str_detect(files, "comb_res")]

all_res <- data.frame()
for (i in 1:length(files)) {
	res <- readRDS(files[i])
	print(dim(res)[1])
	all_res <- rbind(all_res, res)
}

saveRDS(all_res, "comb_res_nova.rds")
