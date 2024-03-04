
#### MERGE AND PLOT HAPLOTYPE ANALYSIS RESULTS FOR ONE CHROMOSOME ####

library(qqman)
library(Haplin)

# Set path for result files
chr <- "04"
path <- paste0(/PATH_TO_CODE_AND_RESULTS/, chr, "/")

# Read and merge result files
resFiles <- list.files(path)
resFiles <- resFiles[which(endsWith(resFiles, ".rds"))]
res <- NULL
for(i in 1:length(resFiles)){
  resnew <- readRDS(file = paste0(path, resFiles[i]))
  res <- rbind(res, resnew)
}

# Save merged file
saveRDS(res, file = paste0(path, "res", chr, ".rds"))

# Delete res files if merged file exists
checkres <- list.files(path)
if(length(which(startsWith(checkres, "res"))) > 0){
  invisible(file.remove(file = paste0(path, resFiles)))
}

# Prep plot object
res_pl <- as.data.frame(subset(res, select = c(CHR, haploblock, haplomiddle, BP, RR.p.value, pv.overall), subset = (!is.na(RR.p.value) & !is.na(pv.overall))))
res_pl <- res_pl[-which(duplicated(res_pl[,c("haploblock", "pv.overall")])),]

# Save plot object
saveRDS(res, file = paste0(path, "res_pl", chr, ".rds"))

# Plot
ymax <- max((-log10(min(res_pl$pv.overall)) + 1), 8)
jpeg(file = paste(path, "fet_haplo_manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res_pl, # xlim = c(xmin, xmax),
          p = "pv.overall",
          snp = "haplomiddle",
          ylim = c(0, ymax),
          col = c("#7E1182"),
          cex = 0.8,
          suggestiveline = -log10(5e-06),
          annotatePval = 5e-06, annotateTop = FALSE)
dev.off()
jpeg(file = paste0(path, "fet_haplo_qq", chr, ".jpeg"))
pQQ(res_pl$pv.overall, mark = FALSE)
dev.off()


