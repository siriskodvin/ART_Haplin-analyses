
#### MERGE AND PLOT RESULTS FOR ONE CHROMOSOME ####

library(qqman)
library(Haplin)

# Input
chr <- "21"

# Read and merge res-files and save
path <- paste0(/PATH_TO_CODE_AND_RESULTS/, chr)
resFiles <- list.files(path)
resFiles <- resFiles[endsWith(resFiles, "rds")]
res_fetmatl <- list()
res_matpool <- list()
res_pool <- list()
if(length(resFiles) != 0){
  for(j in 1:length(resFiles)){
    pathFil <- paste0(path, "/", resFiles[j])
    res <- readRDS(file = pathFil)
    res_fetmatl[[length(res_fetmatl) + 1]] <- res$res_fetmat
    res_matpool[[length(res_matpool) + 1]] <- res$res_matpoo
    res_pool[[length(res_pool) + 1]] <- res$res_poo
    invisible(file.remove(file = pathFil))
  }
}
res_fetmat <- do.call("rbind", res_fetmatl)
dinds <- which(duplicated(res_fetmat))
if(length(dinds) != 0){
  res_fetmat <- res_fetmat[-dinds,]
}
res_matpoo <- do.call("rbind", res_matpool)
dinds <- which(duplicated(res_matpoo))
if(length(dinds) != 0){
  res_matpoo <- res_matpoo[-dinds,]
}
res_poo <- do.call("rbind", res_pool)
dinds <- which(duplicated(res_poo))
if(length(dinds) != 0){
  res_poo <- res_poo[-dinds,]
}
res <- list(res_fetmat, res_matpoo, res_poo)
names(res) <- c("res_fetmat", "res_matpoo", "res_poo")

# Save results
saveRDS(res, file = paste0(path, "/res_chr", chr, ".rds"))

# Save number of variants (save 0 if the numbers differ)
numvar <- unique(c(length(unique(res$res_fetmat$SNP)), length(unique(res$res_matpoo$SNP)), length(unique(res$res_poo$SNP))))
if(length(numvar) > 1){
  numvar <- 0
}
writeLines(as.character(numvar), paste0(path, "/numvar", chr, ".txt"))

# Manhattan plot and QQ plot fetal and maternal effect
res_pl <- subset(res$res_fetmat, subset = !is.na(RR.p.value))
res_pl <- subset(res_pl, subset = RR.p.value > 0)
ymax <- max((-log10(min(res_pl$RR.p.value)) + 1), 8)
jpeg(file = paste(path, "/fetmat_manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res_pl, # xlim = c(xmin, xmax),
          p = "RR.p.value",
          ylim = c(0, ymax),
          col = c("darkred"),
          cex = 0.8,
          suggestiveline = -log10(5e-06),
          annotatePval = 5e-06, annotateTop = FALSE)
dev.off()
jpeg(file = paste0(path, "/fetmat_qq", chr, ".jpeg"))
pQQ(res_pl$RR.p.value, mark = FALSE)
dev.off()
res_pl <- subset(res$res_fetmat, subset = !is.na(RRm.p.value))
res_pl <- subset(res_pl, subset = RRm.p.value > 0)
ymax <- max((-log10(min(res_pl$RRm.p.value)) + 1), 8)
jpeg(file = paste(path, "/matfet_manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res_pl, # xlim = c(xmin, xmax),
          p = "RRm.p.value",
          ylim = c(0, ymax),
          col = c("darkred"),
          cex = 0.8,
          suggestiveline = -log10(5e-06),
          annotatePval = 5e-06, annotateTop = FALSE)
dev.off()
jpeg(file = paste0(path, "/matfet_qq", chr, ".jpeg"))
pQQ(res_pl$RRm.p.value, mark = FALSE)
dev.off()

# Manhattan plot and QQ plot maternal and PoO effect
res_pl <- subset(res$res_matpoo, subset = !is.na(RRm.p.value))
res_pl <- subset(res_pl, subset = RRm.p.value > 0)
ymax <- max((-log10(min(res_pl$RRm.p.value)) + 1), 8)
jpeg(file = paste(path, "/matpoo_manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res_pl, # xlim = c(xmin, xmax),
          p = "RRm.p.value",
          ylim = c(0, ymax),
          col = c("#7E1182"),
          cex = 0.8,
          suggestiveline = -log10(5e-06),
          annotatePval = 5e-06, annotateTop = FALSE)
dev.off()
jpeg(file = paste0(path, "/matpoo_qq", chr, ".jpeg"))
pQQ(res_pl$RRm.p.value, mark = FALSE)
dev.off()
res_pl <- subset(res$res_matpoo, subset = !is.na(RRcm_RRcf.p.value))
res_pl <- subset(res_pl, subset = RRcm_RRcf.p.value > 0)
ymax <- max((-log10(min(res_pl$RRcm_RRcf.p.value)) + 1), 8)
jpeg(file = paste(path, "/poomat_manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res_pl, # xlim = c(xmin, xmax),
          p = "RRcm_RRcf.p.value",
          ylim = c(0, ymax),
          col = c("#7E1182"),
          cex = 0.8,
          suggestiveline = -log10(5e-06),
          annotatePval = 5e-06, annotateTop = FALSE)
dev.off()
jpeg(file = paste0(path, "/poomat_qq", chr, ".jpeg"))
pQQ(res_pl$RRcm_RRcf.p.value, mark = FALSE)
dev.off()

# Manhattan plot and QQ plot PoO effect
res_pl <- subset(res$res_poo, subset = !is.na(RRcm_RRcf.p.value))
res_pl <- subset(res_pl, subset = RRcm_RRcf.p.value > 0)
ymax <- max((-log10(min(res_pl$RRcm_RRcf.p.value)) + 1), 8)
jpeg(file = paste(path, "/poo_manh", chr, ".jpeg", sep = ""),
     width = 1200, height = 600, quality = 100)
par(mar = c(4, 18, 11, 20))
manhattan(res_pl, # xlim = c(xmin, xmax),
          p = "RRcm_RRcf.p.value",
          ylim = c(0, ymax),
          col = c("Darkblue"),
          cex = 0.8,
          suggestiveline = -log10(5e-06),
          annotatePval = 5e-06, annotateTop = FALSE)
dev.off()
jpeg(file = paste0(path, "/poo_qq", chr, ".jpeg"))
pQQ(res_pl$RRcm_RRcf.p.value, mark = FALSE)
dev.off()


