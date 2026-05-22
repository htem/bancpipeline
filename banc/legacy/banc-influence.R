################################
### CALCULATE BANC INFLUENCE ###
################################
source("banc/banc-startup.R")
library(influencer)
overwrite <- TRUE
# banc.version set in banc-startup.R
n_cores <- max(1, parallel::detectCores() - 1)
banc.influence.save.path.version <- file.path(banc.influence.save.path,paste0("banc_",banc.version))
banc.influence.save.path.version.reverse <- file.path(banc.influence.save.path,paste0("banc_",paste0(banc.version,"_reverse")))
dir.create(banc.influence.save.path.version)
# NOTE: stray "in $HOME" fragment removed 2026-04-16 during v850->v888 migration
# (pre-existing parse error; not exercised because script wasn't sourcing).
# Versioned data directory (matches banc-data.R save.path pattern)
save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                        paste0("banc_", banc.version))

#################
### BANC META ###
#################
banc.meta <- arrow::read_feather(
  file.path(save.path, sprintf("banc_%s_meta.feather", banc.version)))
nrow(banc.meta)

# Read edgelist. Default source is v3; override with BANC_SYN_SOURCE or --source.
syn_source <- {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == "--source")
  if (length(i) == 1 && length(args) >= i + 1) tolower(args[i + 1])
  else {
    env <- Sys.getenv("BANC_SYN_SOURCE", unset = NA_character_)
    if (!is.na(env) && nzchar(env)) tolower(env)
    else if (exists("banc.synapse.source.default")) tolower(banc.synapse.source.default)
    else "v3"
  }
}
stopifnot(syn_source %in% c("v2", "v3"))
banc.elist <- arrow::read_feather(
  file.path(save.path, sprintf("banc_%s_edgelist_simple_%s.feather",
                                banc.version, syn_source))) %>%
  dplyr::filter(count > 3)
nrow(banc.elist)

# Reversed
banc.elist.rev <- banc.elist %>%
  dplyr::rename(pre_rev = post,
                post_rev = pre) %>%
  dplyr::rename(pre = pre_rev,
                post = post_rev)

# Seeds to read
seed.columns <- colnames(banc.meta)[grepl("seed_",colnames(banc.meta))]

###########################################
### CALCULATE INFLUENCE FOR SEED GROUPS ###
###########################################

# Create the R calculator and pre-compute cached matrices
# Using influence_calculator_r() (pure R) instead of influence_calculator_py()
# so that forked workers can safely inherit the cached matrix decomposition
cat("Creating influence calculator and pre-computing matrix decomposition...\n")
icr <- influence_calculator_r(edgelist_simple = banc.elist, meta = banc.meta)

# Warm up cache: computes and caches W_normalized, W_factorization, max_eigenvalue
# Forked workers inherit these via copy-on-write shared memory
cat("Pre-computing cached matrix decomposition (this takes time but speeds up all subsequent calculations)...\n")
icr.test <- icr$calculate_influence(banc.meta$root_id[1])
cat("Cache computation complete! All subsequent calculations will be much faster.\n")
print(head(icr.test))

cat(sprintf("\n=== PARALLEL PROCESSING (%d cores, cached decomposition) ===\n", n_cores))

for(sg in seed.columns){
  cat(sprintf("Processing seed group: %s\n", sg))

  # Get all the seed groups
  banc.meta.sg <- banc.meta[!is.na(banc.meta[[sg]]),]
  sg.entries <- unique(banc.meta.sg[[sg]])

  # Create output directory
  banc.influence.save.path.version.seed <- file.path(banc.influence.save.path.version, sg)
  dir.create(banc.influence.save.path.version.seed, showWarnings = FALSE)

  # Parallel processing using forked workers
  # Each worker inherits the cached LU factorisation (read-only, zero-copy)
  parallel::mclapply(sg.entries, function(b) {
    seed.ids <- na.omit(unique(banc.meta[banc.meta[[sg]] == b, "root_id"]))
    if(length(seed.ids) > 0) {
      icr.run <- icr$calculate_influence(seed.ids)
      readr::write_csv(icr.run, file = file.path(banc.influence.save.path.version.seed, paste0(b,".csv")))
    }
  }, mc.cores = n_cores)

  cat(sprintf("Completed seed group %s: %d calculations\n", sg, length(sg.entries)))
}

cat("\n=== ALL FORWARD PROCESSING COMPLETE ===\n")


####################################################
### CALCULATE INFLUENCE FOR SEED GROUPS REVERSED ###
####################################################

# Create a new R calculator for the reversed edgelist
cat("Creating reversed influence calculator and pre-computing matrix decomposition...\n")
icr <- influence_calculator_r(edgelist_simple = banc.elist.rev, meta = banc.meta)

# Warm up cache
cat("Pre-computing cached matrix decomposition (reversed)...\n")
icr.test <- icr$calculate_influence(banc.meta$root_id[1])
cat("Cache computation complete!\n")
print(head(icr.test))

cat(sprintf("\n=== PARALLEL PROCESSING REVERSED (%d cores, cached decomposition) ===\n", n_cores))

for(sg in seed.columns){
  cat(sprintf("Processing seed group: %s\n", sg))

  # Get all the seed groups
  banc.meta.sg <- banc.meta[!is.na(banc.meta[[sg]]),]
  sg.entries <- unique(banc.meta.sg[[sg]])

  # Create output directory
  banc.influence.save.path.version.seed <- file.path(banc.influence.save.path.version.reverse, sg)
  dir.create(banc.influence.save.path.version.seed, showWarnings = FALSE)

  # Parallel processing using forked workers
  parallel::mclapply(sg.entries, function(b) {
    seed.ids <- na.omit(unique(banc.meta[banc.meta[[sg]] == b, "root_id"]))
    if(length(seed.ids) > 0) {
      icr.run <- icr$calculate_influence(seed.ids)
      readr::write_csv(icr.run, file = file.path(banc.influence.save.path.version.seed, paste0(b,".csv")))
    }
  }, mc.cores = n_cores)

  cat(sprintf("Completed seed group %s: %d calculations\n", sg, length(sg.entries)))
}

cat("\n=== ALL PROCESSING COMPLETE ===\n")
