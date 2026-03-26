.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.4-3.20")

cat("Loading packages...\n"); flush.console()
library(bsseq)
library(DelayedMatrixStats)
library(GenomicRanges)
library(BiocParallel)
library(rtracklayer)
library(methylCC)
library(tidyr)
library(readr)



# ------------------------
# Define inputs and paths
# ------------------------
hg38_ranges_rds_path <- "/home/rstudio/methylcc_test/gr_split.rds"
chain_file_path      <- "/home/rstudio/methylcc_test/hg19ToHg38.over.chain"
file_list_path       <- "/home/rstudio/methylcc_test/methylcc_13577HG/13577_file_list.txt"  # one GCS path per line

output_dir            <- "/home/rstudio/methylcc_test/methylcc_13577HG/"

# coverage_threshold <- 1
# n_groups <- 1



# ------------------------
# Liftover Reference
# ------------------------

gr_split = readRDS(hg38_ranges_rds_path)

chain_19_to_38 = import.chain(chain_file_path)
hg38_gr_split = gr_split %>%
  rtracklayer::liftOver(chain_19_to_38) %>%
  unlist()




library(BiocParallel)
library(R.utils)

# Set up parallel backend
register(MulticoreParam(workers = 4))  # adjust to number of CPUs

# ------------------------
# Start file-wise processing (10 at a time)
# ------------------------
file_list <- readLines(file_list_path)
batch_size <- 10
bs_list <- list()

for (i in seq(1, length(file_list), by = batch_size)) {
  batch_paths <- file_list[i:min(i + batch_size - 1, length(file_list))]
  batch_filenames <- basename(batch_paths)
  batch_local_paths <- file.path(output_dir, batch_filenames)
  
  # ------------------------
  # Check existing files and download only missing ones
  # ------------------------
  cat("Downloading batch", i, "to", i + length(batch_paths) - 1, "...\n"); flush.console()
  missing <- !file.exists(batch_local_paths)
  if (any(missing)) {
    writeLines(batch_paths[missing], con = "gcs_batch_list.txt")
    system(paste("gsutil -m cp -I", shQuote(output_dir), "< gcs_batch_list.txt"))
  } else {
    cat("All files already exist. Skipping download.\n")
  }
  
  # ------------------------
  # Confirm all files exist
  # ------------------------
  valid_paths <- batch_local_paths[file.exists(batch_local_paths)]
  if (length(valid_paths) == 0) {
    warning("No valid files found in this batch. Skipping.")
    next
  }
  
  # ------------------------
  # Parallel read BSseq objects
  # ------------------------
  cat("Reading batch into BSseq object...\n"); flush.console()
  bs_objects <- bplapply(valid_paths, function(p) {
    tryCatch({
      read.bismark(files = p, loci = hg38_gr_split, verbose = FALSE)
    }, error = function(e) NULL)
  })
  bs_objects <- Filter(Negate(is.null), bs_objects)
  
  if (length(bs_objects) > 0) {
    bs_list[[length(bs_list) + 1]] <- do.call(combine, bs_objects)
  }
  
  # ------------------------
  # Zip and move source files
  # ------------------------
  for (local_path in valid_paths) {
    cat("Zipping and moving", basename(local_path), "...\n"); flush.console()
    gzipped_path <- paste0(local_path, ".gz")
    gzip(local_path, destname = gzipped_path, overwrite = TRUE)
    file.rename(gzipped_path, file.path(output_dir, "source", basename(gzipped_path)))
  }
}

# ------------------------
# Combine BSseq objects
# ------------------------
cat("Merging BSseq objects...\n"); flush.console()
bs_all <- do.call(combine, bs_list)

saveRDS(bs_all, paste0(output_dir, "concat_data_13577_hg38.rds"))


# ------------------------
# Filter
# ------------------------

filter = function(input_data, coverage = 1, n_groups = 10) {
  bs = GenomeInfoDb::keepStandardChromosomes(input_data, pruning.mode = "coarse")
  n_samples = nrow(input_data@colData)
  loci.idx = which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") >= coverage) >= n_groups)
  bs.filtered = bs[loci.idx, ]
  
  return(bs.filtered)
}

#cf_data = readRDS("/home/rstudio/methylcc_test/methylcc_tmp/cf_test50_bsseq_hg38.rds")

cf_filtered = filter(bs_all)

# ------------------------
# LiftOver
# ------------------------

liftover_prep_for_methylCC = function(filtered_bs) {
  
  # Make indices
  mcols(filtered_bs)$index = 1:length(filtered_bs)
  hg38 = rowRanges(filtered_bs)
  hg38$index = 1:length(hg38)
  
  # Liftover ranges
  hg19 = hg38 %>%
    rtracklayer::liftOver(AnnotationHub::AnnotationHub()[["AH14108"]]) %>%
    unlist()
  
  # Subset based on regions that lifted over
  bs_hg19 = filtered_bs[which(hg19$index %in% hg38$index)]
  
  # Assign lifted over regions
  rowRanges(bs_hg19) = hg19
  
  glue::glue("{length(bs_hg19)} out of {length(filtered_bs)} were lifted over")
  glue::glue("{length(filtered_bs) - length(bs_hg19)} did not liftOver")
  
  return(bs_hg19)
}

cf_hg19 = liftover_prep_for_methylCC(cf_filtered)


# ------------------------
# Run methylCC
# ------------------------

cf_est = estimatecc(cf_hg19, include_cpgs = TRUE, include_dmrs = TRUE)

# ------------------------
# Save results
# ------------------------

cf_methylCC = gather(cbind("samples" = rownames(cell_counts(cf_est)), cell_counts(cf_est)), celltype, est, -samples)
write_tsv(cf_methylCC, paste0(output_dir, "13577_methylCC_counts.tsv"))
