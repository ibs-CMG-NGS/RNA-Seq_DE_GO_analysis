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

# 기존 metadata.json 있으면 병합 후 alias 기준으로 최신 import_date만 남긴다.
# dataset_id는 실행할 때마다 새로 생성되는 UUID라 이것만으로 중복 제거하면, 같은 항목을
# 다시 만들어내는 재실행(예: 다른 옵션 추가로 일부 rule만 다시 도는 경우)이 있을 때마다
# 예전 항목이 지워지지 않고 계속 쌓이는 문제가 실제로 있었다. alias는 재실행해도 그대로
# 유지되는 실질적 식별자이므로 이 기준으로 병합하면 예전에 이미 쌓인 stale 중복도
# 다음 실행 때 자동으로 정리된다(self-healing).
existing_datasets <- list()
if (file.exists(meta_path)) {
  existing <- fromJSON(meta_path, simplifyVector = FALSE)
  existing_datasets <- existing$datasets
  cat(paste("Existing metadata.json found:", length(existing_datasets), "datasets\n"))
}

combined <- c(existing_datasets, new_entries)
by_alias <- list()
for (e in combined) {
  a <- e$alias
  if (is.null(by_alias[[a]]) || e$import_date > by_alias[[a]]$import_date) {
    by_alias[[a]] <- e
  }
}
all_entries <- unname(by_alias)
cat(paste("After alias-based dedup:", length(all_entries), "datasets (",
          length(combined) - length(all_entries), "stale duplicates removed )\n"))

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
