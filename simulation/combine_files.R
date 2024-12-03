
args <- commandArgs()
folder <- args[6]
dist <- args[7]

print(folder)
print(dist)


loc <- paste0(folder, "/", dist, "/")
#files <- list.files(loc, pattern = "\\.csv")
files <- list.files(loc, patter = "\\.rds")
files <- paste0(loc, files)

#tables <- lapply(files, read.csv)
tables <- lapply(files, readRDS)
comb_res <- do.call(rbind, tables)

file_name <- paste0("comb_res/", gsub("sim_", "", folder), "_", dist, ".rds")
saveRDS(comb_res, file_name)





















