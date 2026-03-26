# libpath <- ""

library(Biostrings)
library(bsseq)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get CpG sites for this chr 
args <- commandArgs(trailingOnly = TRUE)
print(args[1])
chr <- names(Hsapiens)[as.numeric(args[1])]
print(chr)
cgs <- lapply(chr, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(seq_along(chr), function(x) GRanges(chr[x], IRanges(cgs[[x]], width = 2))))

# Twist capture regions
twist_regions<-read.delim("/path/to/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed",header = F)
twist_regions<-twist_regions[,1:3]
colnames(twist_regions) <- c("chr","start","end")
twist_regions2<-twist_regions
twist_regions2$start<-twist_regions$start-200
twist_regions2$end<-twist_regions$end+200

sites_in_twist <- as.data.frame(GenomicRanges::intersect(
  x = cpgr,
  y = GRanges(twist_regions2)
))
sites_in_twist2<-subset(sites_in_twist,width==2)
sites_in_twist2$end<-sites_in_twist2$start

sites_in_twist3<-sites_in_twist2
sites_in_twist3$start<-sites_in_twist3$end<-sites_in_twist2$end+1

sites_in_twist2$strand<-'+'
sites_in_twist3$strand<-'-'
sites_in_twist<-rbind(sites_in_twist2,sites_in_twist3)
sites_in_twist<-sites_in_twist[order(sites_in_twist$seqnames,sites_in_twist$start),]
sites_in_twist$width<-1


# paths
cg_dir  <- paste0("/path/to/CG_reports/",chr)
tmp_dir <- file.path("/path/tmp", chr)
outdir  <- file.path("/path/to/bsseq_objects/bsseq", paste0("_", chr))

# Run bsseq just for this chromosome
bismarkBSseq <- read.bismark(
  files = list.files(cg_dir, full.names = TRUE, pattern = ".CG_report.txt"),
  strandCollapse = TRUE,
  loci = GRanges(sites_in_twist),
  verbose = TRUE,
  BACKEND = "HDF5Array",
  dir = tmp_dir,
  replace = TRUE,
  nThread = 8
)

# Save properly
library(HDF5Array)
saveHDF5SummarizedExperiment(
  x = bismarkBSseq,
  dir = outdir,
  replace = TRUE
)

unlink(tmp_dir, recursive = TRUE)

print("Saved")