#### -------------------------- Script header -------------------------- ####
# Date:         04 March 2024                                               #
# Author:       Siri Naerland Skodvin                                       #
# Filename:     FRhap.R                                                     #
# Description:  This script reads in genetic data and family structure      #
#               indices and does analyses with Haplin                       #
#### --------------------------------------------------------------------####



#### INPUT, DATA AND PACKAGES/FUNCTIONS ####

library(data.table)   # To read ped and map file
library(dplyr)        # To merge (adjust) ped file with inds info
library(Haplin)       # To do analyses
library(parallel)     # To run analyses in parallel

# Set chromosome number - CHANGE HERE
chr <- "21"
hapcpus <- 4

# Design and phenotype - CHANGE HERE
design <- "triad"
dname <- "inds_tridy_art"

# Paths
cpath <- /PATH_TO_CODE_AND_RESULTS/
dpath <- /PATH_TO_TEMPORARY_STORAGE_AND_MAIN_R_SCRIPTS/

# Read file with indices to use for analyses
inds_list <- readRDS(file = /FAMILY_KEY_FILE.rds/)
inds <- inds_list[[which(names(inds_list) == dname)]]

# Identify the genetic files
filenames <- list.files(dpath)
filenames <- c(filenames[which(endsWith(filenames, ".ped"))],
               filenames[which(endsWith(filenames, ".map"))])
filenames <- paste0(dpath, "/", filenames)

# Read ped file, subset on sample inds, remove original file and save new
ped <- as.data.frame(fread(file = filenames[1]))
ped <- ped[inds$inds,]
ped$V6 <- inds$pheno
# invisible(file.remove(file = filenames[1]))
fwrite(ped, file = paste0(substr(filenames[1], 1, nchar(filenames[1]) - 4), "_samp.ped"),
       sep = " ", col.names = FALSE)

# Read map file
map <- as.data.frame(fread(file = filenames[2]))
map <- subset(map, select = c(V1, V2, V4))
colnames(map) <- c("CHR", "SNP", "BP")
fwrite(map, file = paste0(substr(filenames[1], 1, nchar(filenames[1]) - 4), "_samp.map"),
       sep = " ", col.names = TRUE)

# Read and preprocess files in Haplin
hapgen <- genDataRead(file.in = paste0(substr(filenames[1], 1, nchar(filenames[1]) - 4), "_samp.ped"),
                      dir.out = paste0(dpath, "/data"),
                      format = "ped",
                      map.file = paste0(substr(filenames[1], 1, nchar(filenames[1]) - 4), "_samp.map"),
                      overwrite = TRUE)
hapgenpres <- genDataPreprocess(hapgen, design = design,
                                dir.out = paste0(dpath, "/data"),
                                ncpu = hapcpus,
                                overwrite = TRUE)

# Save timestamp to save along with results
t <- Sys.time()
t <- gsub(" ", "_", t)

# Haplin - fetal and maternal effect
res_fetmat <- haplinSlide(data = hapgenpres, markers = "ALL", design = design, use.missing = TRUE,
                          maternal = TRUE, poo = FALSE, reference = "ref.cat", response = "mult",
                          threshold = 0.01, data.out = "no", winlength = 1, table.output = TRUE,
                          cpus = hapcpus, para.env = "parallel",
                          slaveOutfile = paste0(cpath, "/Results/Chr", chr, "/", t, "_hapout_fetmat_chr", chr, ".txt"),
                          printout = FALSE, verbose = FALSE)
res_fetmat <- do.call("rbind", res_fetmat)
markernames <- do.call("rbind", strsplit(rownames(res_fetmat), split = '.', fixed = TRUE))
res_fetmat$marker <- markernames[,1]
colnames(res_fetmat)[which(colnames(res_fetmat) == "marker")] <- "SNP"
res_fetmat <- left_join(res_fetmat, map, by = c("SNP" = "SNP"))
# res_fetmat <- right_join(map, res_fetmat, by = c("m" = "marker"))

# Haplin - maternal and PoO effect
res_matpoo <- haplinSlide(data = hapgenpres, markers = "ALL", design = design, use.missing = TRUE,
                          maternal = TRUE, poo = TRUE, reference = "ref.cat", response = "mult",
                          threshold = 0.01, data.out = "no", winlength = 1, table.output = TRUE,
                          cpus = hapcpus, para.env = "parallel",
                          slaveOutfile = paste0(cpath, "/Results/Chr", chr, "/", t, "_hapout_matpoo_chr", chr, ".txt"),
                          printout = FALSE, verbose = FALSE)
res_matpoo <- do.call("rbind", res_matpoo)
markernames <- do.call("rbind", strsplit(rownames(res_matpoo), split = '.', fixed = TRUE))
res_matpoo$marker <- markernames[,1]
colnames(res_matpoo)[which(colnames(res_matpoo) == "marker")] <- "SNP"
res_matpoo <- left_join(res_matpoo, map, by = c("SNP" = "SNP"))
# res_matpoo <- right_join(map, res_matpoo, by = c("m" = "marker"))

# Haplin - PoO effect
res_poo <- haplinSlide(data = hapgenpres, markers = "ALL", design = design, use.missing = TRUE,
                       maternal = FALSE, poo = TRUE, reference = "ref.cat", response = "mult",
                       threshold = 0.01, data.out = "no", winlength = 1, table.output = TRUE,
                       cpus = hapcpus, para.env = "parallel",
                       slaveOutfile = paste0(cpath, "/Results/Chr", chr, "/", t, "_hapout_poo_chr", chr, ".txt"),
                       printout = FALSE, verbose = FALSE)
res_poo <- do.call("rbind", res_poo)
markernames <- do.call("rbind", strsplit(rownames(res_poo), split = '.', fixed = TRUE))
res_poo$marker <- markernames[,1]
colnames(res_poo)[which(colnames(res_poo) == "marker")] <- "SNP"
res_poo <- left_join(res_poo, map, by = c("SNP" = "SNP"))
# res_poo <- right_join(map, res_poo, by = c("m" = "marker"))

res <- list(res_fetmat, res_matpoo, res_poo)
names(res) <- c("res_fetmat", "res_matpoo", "res_poo")

# Save results
saveRDS(res, file = paste0(cpath, "/Results/Chr", chr, "/", t, "_hapout_chr", chr, ".rds"))


