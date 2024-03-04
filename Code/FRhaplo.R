
#### INPUT, DATA AND PACKAGES/FUNCTIONS ####

library(data.table)   # To read ped and map file
library(dplyr)        # To merge (adjust) ped file with inds info
library(Haplin)       # To do analyses
library(parallel)     # To run analyses in parallel
library(stringr)

# Set chromosome number - CHANGE HERE
chr <- "04"
hapcpus <- 4

# Effects - CHANGE HERE
mat <- FALSE
poo <- FALSE

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

# Do analysis in Haplin with haplotypes of a given length
hapl <- 3
res <- haplinSlide(data = hapgenpres, markers = "ALL", design = design, use.missing = TRUE,
                   maternal = mat, poo = poo, reference = "ref.cat", response = "mult",
                   threshold = 0.01, data.out = "no", winlength = hapl, table.output = TRUE,
                   cpus = hapcpus, para.env = "parallel",
                   slaveOutfile = paste0(cpath, "/Results/Chr", chr, "/", t, "_hapout_chr", chr, ".txt"),
                   printout = FALSE, verbose = FALSE)

# Merge results
res <- do.call("rbind", res)

# Add SNP and haplowindow information
res$SNP <- rownames(res)
rownames(res) <- NULL
res$haplowindow <- sapply(res$SNP, function(x){sub("\\..*", "", x)})

# Add BP information for middle SNP
res$haplomiddle <- sapply(res$haplowindow, function(x){str_extract(x, "(?<=-).*(?=-)")})
res <- left_join(res, map, by = c("haplomiddle" = "SNP"))

# Save results
saveRDS(res, file = paste0(cpath, "/Results/Chr", chr, "/", t, "_hapout_hapl", hapl, "_chr", chr, ".rds"))


