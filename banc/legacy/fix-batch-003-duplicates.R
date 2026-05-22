###########################################################
### fix-batch-003-duplicates.R
###
### Removes the 3229 duplicate rows (rows 1-3229) from
### batch_000003.parquet that are already present in
### batch_000001 and batch_000002.
###
### Idempotent: skips if .bak already exists.
###########################################################

library(arrow)
library(dplyr)

batch_dir <- "/n/data1/hms/neurobio/wilson/banc/synapses_v3/processed"

batch3_path <- file.path(batch_dir, "batch_000003.parquet")
bak_path    <- paste0(batch3_path, ".bak")

# --- Idempotent guard ---
if (file.exists(bak_path)) {
  message("Backup already exists at ", bak_path)
  message("Skipping -- fix has already been applied.")
  quit(save = "no", status = 0)
}

stopifnot(file.exists(batch3_path))

# --- Read syn_ids from batch_000001 and batch_000002 ---
message("Reading syn_ids from batch_000001 and batch_000002...")
batch1 <- read_parquet(file.path(batch_dir, "batch_000001.parquet"),
                       col_select = "syn_id")
batch2 <- read_parquet(file.path(batch_dir, "batch_000002.parquet"),
                       col_select = "syn_id")
prior_syn_ids <- unique(c(batch1$syn_id, batch2$syn_id))
message(sprintf("  %s unique syn_ids in prior batches",
                format(length(prior_syn_ids), big.mark = ",")))
rm(batch1, batch2)

# --- Read batch_000003 ---
message("Reading batch_000003...")
batch3 <- read_parquet(batch3_path)
n_before <- nrow(batch3)
message(sprintf("  %s rows before cleaning", format(n_before, big.mark = ",")))

# --- Remove duplicates ---
is_dup <- batch3$syn_id %in% prior_syn_ids
n_dup <- sum(is_dup)
message(sprintf("  %d rows have syn_ids present in prior batches", n_dup))

if (n_dup == 0) {
  message("No duplicates found -- nothing to do.")
  quit(save = "no", status = 0)
}

batch3_clean <- batch3[!is_dup, , drop = FALSE]
n_after <- nrow(batch3_clean)

# --- Backup original ---
message(sprintf("Backing up original to %s", bak_path))
file.copy(batch3_path, bak_path)
stopifnot(file.exists(bak_path))

# --- Write cleaned version ---
message("Writing cleaned batch_000003.parquet...")
write_parquet(batch3_clean, batch3_path)

message(sprintf("\nDone. Removed %d duplicate rows.", n_dup))
message(sprintf("  Before: %s rows", format(n_before, big.mark = ",")))
message(sprintf("  After:  %s rows", format(n_after, big.mark = ",")))
message(sprintf("  Backup: %s", bak_path))
