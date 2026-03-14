# 파일 경로: src/analysis/06b_aggregate_seqviewer.R
# 모든 pair의 staging JSON을 병합해 metadata.json 생성
# Usage: Rscript 06b_aggregate_seqviewer.R [seqviewer_dir]

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 06b_aggregate_seqviewer.R [seqviewer_dir]")
}
seqviewer_dir <- args[1]
staging_dir   <- file.path(seqviewer_dir, "staging")
meta_path     <- file.path(seqviewer_dir, "metadata.json")

# staging/*.json 전체 수집
staging_files <- list.files(staging_dir, pattern = "_entries\\.json$", full.names = TRUE)
if (length(staging_files) == 0) {
  stop(paste("No staging files found in", staging_dir))
}

new_entries <- list()
for (f in staging_files) {
  entries <- fromJSON(f, simplifyVector = FALSE)
  new_entries <- c(new_entries, entries)
}
cat(paste("Collected", length(new_entries), "entries from", length(staging_files), "staging files\n"))

# 기존 metadata.json 있으면 병합 (dataset_id 기준 중복 제거, 새 항목 우선)
existing_datasets <- list()
if (file.exists(meta_path)) {
  existing <- fromJSON(meta_path, simplifyVector = FALSE)
  existing_datasets <- existing$datasets
  cat(paste("Existing metadata.json found:", length(existing_datasets), "datasets\n"))
}

new_ids    <- sapply(new_entries, `[[`, "dataset_id")
kept_old   <- Filter(function(x) !x$dataset_id %in% new_ids, existing_datasets)
all_entries <- c(kept_old, new_entries)

metadata <- list(
  version      = "1.0",
  last_updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  datasets     = all_entries
)

write_json(metadata, meta_path, pretty = TRUE, auto_unbox = TRUE)
cat(paste("metadata.json saved:", meta_path,
          "(", length(all_entries), "total datasets )\n"))

# 완료 플래그
flag_path <- file.path(seqviewer_dir, ".seqviewer_done.flag")
file.create(flag_path)
cat("seqviewer aggregation done.\n")
