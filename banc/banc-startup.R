#' banc-startup — Project startup: load libraries, set options, define save paths.
#'
#' Sourced by every BANC pipeline script as `source("banc/banc-startup.R")`.
#' Detects the host machine, loads the natverse + tidyverse + arrow + parallel
#' stack, sets gargle / fafbseg / pillar options, sources `banc-functions.R`,
#' loads private credentials via `banc/load-keys.R`, and defines the
#' canonical save-path constants (`banc.save.path`, `banc.connectivity.save.path`,
#' `banc.meta.save.path`, `banc.l2split.save.path`, `banc.versioned.save.path`,
#' etc.) plus the global `banc.version` and `banc.synapse.source.default`.
#'
#' @section Reads:
#'   - env vars `BANC_VERSION` (overrides default 888),
#'     `BANC_SYN_SOURCE` (default v3), seatable / google tokens
#'   - `banc/banc-functions.R` (helper functions)
#'   - `banc/load-keys.R` (private credentials)
#'
#' @section Notes:
#'   - Off-O2 the absolute paths (e.g. `/n/data1/...`) do not resolve;
#'     downstream scripts that need them must short-circuit or fetch from GCS.
#'   - `banc.version` is the single point of truth for which BANC
#'     materialisation everything downstream targets.

####################
##### IDENTITY #####
####################

# What machine are we using?
user<-Sys.info()["user"]
machine<-"o2"

###################
##### OPTIONS #####
#################### 

# google account authentication
options(gargle_oauth_email="alexander.shakeel.bates@gmail.com")
options(gargle_oob_default=TRUE)
options(pillar.sigfig=15)

# path to some FAFB synapse data
options(fafbseg.sqlitepath="/n/data1/hms/neurobio/wilson/fafbz/fafbsynapses")

#####################
##### LIBRARIES #####
#####################

# load required libraires
options(scipen = 999)
library(bancr)
hemibrainr:::suppress(library(hemibrainr))
hemibrainr:::suppress(library(malevnc))
hemibrainr:::suppress(library(nat.nblast))
hemibrainr:::suppress(library(fafbseg))
hemibrainr:::suppress(library(jsonlite))
hemibrainr:::suppress(library(foreach))
hemibrainr:::suppress(library(nat.jrcbrains))
hemibrainr:::suppress(library(doMC))
hemibrainr:::suppress(library(doParallel))
hemibrainr:::suppress(library(progressr))
hemibrainr:::suppress(library(googledrive))
hemibrainr:::suppress(library(elmr))
hemibrainr:::suppress(library(dplyr))
hemibrainr:::suppress(library(tidyverse))
hemibrainr:::suppress(library(bit64))
hemibrainr:::suppress(library(reticulate))
hemibrainr:::suppress(library(RSQLite))  # still needed for influence SQLite databases and legacy reads
hemibrainr:::suppress(library(arrow))
hemibrainr:::suppress(library(plyr))
hemibrainr:::suppress(library(slackr))
hemibrainr:::suppress(library(ggforce))
hemibrainr:::suppress(library(natcpp))
hemibrainr:::suppress(library(lubridate))
hemibrainr:::suppress(library(googlesheets4))
hemibrainr:::suppress(library(doSNOW))
hemibrainr:::suppress(library(fs))
hemibrainr:::suppress(library(purrr))
hemibrainr:::suppress(library(dplyr))
hemibrainr:::suppress(library(processx))
hemibrainr:::suppress(register_saalfeldlab_registrations())

# Source other custom functions
source("banc/banc-functions.R")

# Load private identifiers (Google Sheet / Doc / Drive IDs) into banc.keys.
# Soft fail: missing data/private/keys.csv leaves banc.keys empty so
# scripts that don't need these IDs are unaffected.
source("banc/load-keys.R")

# Session-level cache for banctable_query to avoid redundant API calls.
# banctable_query_cached() returns cached results for identical SQL within
# max_age_secs (default 10 min). Wraps banctable_query() transparently.
# Cache lives in .banctable_cache env; clear with banctable_cache_clear().
# Falls back to latest snapshot CSV if API fails or if running on a parallel worker.
# IMPORTANT: Scripts that UPDATE seatable must use banctable_query() directly,
# NOT banctable_query_cached(), to ensure fresh data before writes.
if (!exists(".banctable_cache", envir = .GlobalEnv)) {
  .banctable_cache <- new.env(parent = emptyenv())
  assign(".banctable_cache", .banctable_cache, envir = .GlobalEnv)
}

# Helper: read latest snapshot CSV from banc.meta.save.path/snapshots/
.banctable_read_snapshot <- function() {
  snapshot_dir <- file.path(banc.meta.save.path, "snapshots")
  if (!dir.exists(snapshot_dir)) return(NULL)
  snapshots <- sort(list.files(snapshot_dir, pattern = "_banc_seatable\\.csv$",
                                full.names = TRUE), decreasing = TRUE)
  if (length(snapshots) == 0) return(NULL)
  message(sprintf("  [snapshot] Reading %s", basename(snapshots[1])))
  readr::read_csv(snapshots[1], show_col_types = FALSE)
}

# Helper: read meta from GCS feather (last resort fallback)
# Maps table names to GCS paths within gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/
# Layout: compiled_data/{dataset}_{version}/{file}
.banctable_read_gcs <- function(table = "banc_meta") {
  gcs_files <- c(
    banc_meta   = "banc_888/banc_888_meta.feather",
    franken_meta = "manc_121/manc_121_meta.feather",
    malecns     = "malecns_09/malecns_09_meta.feather"
  )
  if (!table %in% names(gcs_files)) return(NULL)
  tryCatch({
    gcs_base <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
    gcs_path <- paste0(gcs_base, "/", gcs_files[[table]])
    cache_dir <- file.path(tempdir(), "gcs_cache")
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    local_file <- file.path(cache_dir, basename(gcs_files[[table]]))
    if (!file.exists(local_file)) {
      message(sprintf("  [GCS] Downloading %s via gsutil", basename(local_file)))
      system2("gsutil", c("cp", gcs_path, local_file), stdout = FALSE, stderr = FALSE)
      if (!file.exists(local_file)) stop("gsutil cp failed")
    } else {
      message(sprintf("  [GCS] Using cached %s", basename(local_file)))
    }
    arrow::read_feather(local_file)
  }, error = function(e) {
    warning("GCS fallback failed for ", table, ": ", e$message)
    NULL
  })
}

banctable_query_cached <- function(sql = "SELECT * FROM banc_meta",
                                   max_age_secs = 600, ...) {
  cache <- get(".banctable_cache", envir = .GlobalEnv)
  key <- digest::digest(sql)
  if (exists(key, envir = cache)) {
    entry <- get(key, envir = cache)
    age <- as.numeric(difftime(Sys.time(), entry$time, units = "secs"))
    if (age < max_age_secs) {
      message(sprintf("  [cache hit] %s (%.0fs old)", substr(sql, 1, 60), age))
      return(entry$data)
    }
  }
  # In multi-core environments, prefer snapshot to avoid concurrent API calls
  use_multicore <- exists("numCores") && numCores > 1
  if (use_multicore) {
    result <- .banctable_read_snapshot()
    if (!is.null(result) && nrow(result) > 0) {
      assign(key, list(data = result, time = Sys.time()), envir = cache)
      return(result)
    }
  }
  result <- tryCatch(banctable_query(sql = sql, ...), error = function(e) {
    warning("banctable_query failed: ", e$message)
    NULL
  })
  if (!is.null(result) && nrow(result) > 0) {
    assign(key, list(data = result, time = Sys.time()), envir = cache)
  } else {
    # API failed — fall back to latest snapshot, then GCS
    result <- .banctable_read_snapshot()
    if (is.null(result) || nrow(result) == 0) {
      # Detect table name from SQL for GCS fallback
      tbl <- if (grepl("\\bbanc_meta\\b", sql)) "banc_meta"
             else if (grepl("\\bfranken_meta\\b", sql)) "franken_meta"
             else if (grepl("\\bmalecns\\b", sql)) "malecns"
             else NULL
      if (!is.null(tbl)) result <- .banctable_read_gcs(tbl)
    }
    if (!is.null(result) && nrow(result) > 0) {
      assign(key, list(data = result, time = Sys.time()), envir = cache)
    }
  }
  result
}

banctable_cache_clear <- function() {
  cache <- get(".banctable_cache", envir = .GlobalEnv)
  rm(list = ls(cache), envir = cache)
  message("banctable cache cleared")
}

# get some nifty functions for easy use
`%dopar%` <- foreach::`%dopar%`
`%:%` <- foreach::`%:%`
load_assign <- hemibrainr:::load_assign
overlap_score_delta <- hemibrainr:::overlap_score_delta
check_package_available <- hemibrainr:::check_package_available
nullToNA <- hemibrainr:::nullToNA
# plyr's setup_parallel only WARNS if no backend is registered; it never
# registers one. That meant every %dopar% in this codebase silently fell
# through to doSEQ (1 worker), regardless of SLURM allocation. Wrap it so
# that on first call we register doMC with numCores workers, then call
# plyr's check.
setup_parallel <- function() {
  if (foreach::getDoParWorkers() == 1L &&
      exists("numCores", envir = .GlobalEnv) &&
      get("numCores", envir = .GlobalEnv) > 1L) {
    nc <- get("numCores", envir = .GlobalEnv)
    doMC::registerDoMC(cores = nc)
    message(sprintf("  setup_parallel: registered doMC with %d cores", nc))
  }
  plyr:::setup_parallel()
  invisible(NULL)
}
stop_parallel <- function(cl) {
  if (!is.null(cl) && inherits(cl, "cluster")) parallel::stopCluster(cl)
  invisible(NULL)
}

################################
##### PARRALLEL PROCESSING #####
################################

# Cores: interactive sessions use 1 core; batch respects SLURM or defaults to 10
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
numCores.possible <- ceiling(parallel::detectCores())
if (interactive()) {
  numCores <- 1
  message("  Interactive session: using 1 core")
} else if (nzchar(slurm_cpus)) {
  numCores <- as.integer(slurm_cpus)
  message("  SLURM detected: using ", numCores, " cores (SLURM_CPUS_PER_TASK)")
} else if (is.na(numCores.possible)) {
  numCores <- 1
} else {
  numCores <- min(10, numCores.possible)
}

# If the env var BANC_TEST_IDS_FILE points at a file, expose its contents as
# global banc.test.ids — every per-target NBLAST script honors that variable
# already to restrict nblast.todo. This lets the proofread-redo workflow scope
# all per-target jobs to the same id list without editing each script.
.banc_tids_file <- Sys.getenv("BANC_TEST_IDS_FILE", "")
if (nzchar(.banc_tids_file) && file.exists(.banc_tids_file)) {
  banc.test.ids <- readLines(.banc_tids_file)
  banc.test.ids <- banc.test.ids[nzchar(banc.test.ids)]
  message(sprintf("  banc.test.ids: %d ids loaded from %s",
                  length(banc.test.ids), .banc_tids_file))
}

#################
##### PATHS #####
#################

# Our intended working directory
wd <- "/home/ab714/bancpipeline/"

# Current BANC dataset version
banc.version <- "888"

# Default synapse source for edgelist-building / -reading scripts ("v2" or "v3").
# v2 = synapses_v2_human_readable.csv.gz (has NT predictions).
# v3 = new synapse pipeline (no NT predictions yet).
# Individual scripts may override via --source CLI arg or BANC_SYN_SOURCE env var.
banc.synapse.source.default <- "v3"

# Data storage for BANC
banc_data_storage <- "/n/data1/hms/neurobio/wilson/banc"
banc.save.path <- "/n/data1/hms/neurobio/wilson/banc"
banc.versioned.save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                                       paste0("banc_", banc.version))
rclone.path <- file.path(banc_data_storage,"googledrive")
banc.obj.save.path <- file.path(banc.save.path,"obj")
banc.swc.save.path <- file.path(banc.save.path,"swc")
banc.l2swc.save.path <- file.path(banc.save.path,"l2")
banc.split.save.path <- file.path(banc.save.path,"split")
banc.l2split.save.path <- file.path(banc.save.path,"l2split")
banc.synapses.save.path <- file.path(banc.save.path,"synapses")
banc.metrics.save.path <- file.path(banc.save.path,"metrics")
banc.meta.save.path <- file.path(banc.save.path,"meta")
banc.nt.save.path <- file.path(banc.save.path,"nt")
banc.influence.save.path <- file.path(banc.save.path,"influence")
rda.dir <- file.path(banc.save.path,"banc","rda")
images.dir <- file.path(banc.save.path,"banc","images")
banc.connectivity.save.path <- file.path(banc.save.path,"connectivity")
banc.influence.save.path.version.seed <- file.path(banc.save.path,"influence")

# BANC NBLAST results
bancsynapses <- file.path(banc.save.path,"synapses_250226_human_readable.csv")
bancsynapses.nts <- file.path(banc.save.path,"banc_nt_prediction_w_sizethresh_5_09072025.parquet")
banc.nblast.save.path <- file.path(banc.save.path,"matching")
banc.deform.save.path <- file.path(banc.save.path,"deformetrica")
banc.nblast.mirror.save.path <- file.path(banc.nblast.save.path,"mirror")
banc.nblast.native.save.path <- file.path(banc.nblast.save.path,"native")
banc.nblast.fafb.save.path <- file.path(banc.nblast.save.path,"fafb")
banc.nblast.fafb.obj.save.path <- file.path(banc.nblast.fafb.save.path,"banc_space_obj")
banc.nblast.fafb.swc.save.path <- file.path(banc.nblast.save.path,"fafb_banc_space_swc")
banc.nblast.fanc.save.path <- file.path(banc.nblast.save.path,"fanc")
banc.nblast.fanc.obj.save.path <- file.path(banc.nblast.fanc.save.path,"banc_space_obj")
banc.nblast.fanc.swc.save.path <- file.path(banc.nblast.fanc.save.path,"fanc_banc_space_swc")
banc.nblast.hemibrain.save.path <- file.path(banc.nblast.save.path,"hemibrain")
banc.nblast.hemibrain.obj.save.path <- file.path(banc.nblast.hemibrain.save.path,"banc_space_obj")
banc.nblast.hemibrain.swc.save.path <- file.path(banc.nblast.hemibrain.save.path,"hemibrain_banc_space_swc")
banc.nblast.manc.save.path <- file.path(banc.nblast.save.path,"manc")
banc.nblast.manc.swc.save.path <- file.path(banc.nblast.manc.save.path,"manc_banc_space_swc")
banc.nblast.manc.obj.save.path <- file.path(banc.nblast.manc.save.path,"banc_space_obj")
banc.nblast.manc.split.save.path <- file.path(banc.nblast.manc.save.path,"banc_space_split")
banc.nblast.malecns.save.path <- file.path(banc.nblast.save.path,"malecns")
banc.nblast.malecns.obj.save.path <- file.path(banc.nblast.malecns.save.path,"banc_space_obj")
banc.nblast.malecns.swc.save.path <- file.path(banc.nblast.malecns.save.path,"malecns_banc_space_swc")

# NBLAST transform versions
banc.nblast.version <- "elastix_tpsreg_240721"
banc.nblast.malecns.version <- "navis_tpsreg_250206"

# Manual matching results
banc.mirror.correct.match.path <- file.path(banc.nblast.mirror.save.path,"correct")
banc.fafb.correct.match.path <- file.path(banc.nblast.fafb.save.path,"correct")
banc.manc.correct.match.path <- file.path(banc.nblast.manc.save.path,"correct")
banc.hemibrain.correct.match.path <- file.path(banc.nblast.hemibrain.save.path,"correct")
banc.fanc.correct.match.path <- file.path(banc.nblast.fanc.save.path,"correct")
banc.malecns.correct.match.path <- file.path(banc.nblast.malecns.save.path,"correct")

# Where I have some hemibrain data
fafbz <- "/n/data1/hms/neurobio/wilson/fafbz/"
hemibrain_data <- file.path(fafbz,"hemibrainr_data")
flywire.save.path <- file.path(hemibrain_data,"flywire_neurons")
hemibrain.save.path <- file.path(hemibrain_data,"hemibrain_neurons")
hemibrain.nblast.save.path <- file.path(hemibrain_data,"hemibrain_nblast/")
synister.save.path <- file.path(hemibrain_data,"synister/")
hemibrain.split.save.path <- file.path()

# Where I already have some FAFB-FlyWire data
fafbsynapses <- file.path(fafbz,"fafbsynapses")
flywire.obj.save.path <- file.path(flywire.save.path,"obj")
flywire.swc.save.path <- file.path(flywire.save.path,"783","swc")
flywire.split.save.path <- file.path(flywire.save.path,"783","split")
flywire.synapses.save.path <- file.path(flywire.save.path,"783","synapses")
flywire.metrics.save.path <- file.path(flywire.save.path,"783","metrics")
flywire.nt.save.path <- file.path(flywire.save.path,"783","nt")
flywire.l2.save.path <- file.path(flywire.save.path,"783","l2")

# flycircuit paths
lm.save.path <- file.path(hemibrain_data,"light_level")
lm.nblast.save.path <- paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/")
lm.save.path <- file.path(hemibrain_data,"light_level")
flycircuit.save.path <- file.path(hemibrain_data,"light_level/flycircuit")
lhns.save.path <- file.path(hemibrain_data,"light_level/lhns")

# temporary directory in case we need one
tmpdir <- file.path(banc_data_storage,"rtmp")

######################
##### PARAMETERS #####
######################

# What brainspaces to use
brains <- c("FlyWire", "JRC2018F", "JFRC2")

# Splitting parametersxfbi
polypre <- TRUE
split <- "synapses"
mode <- "centrifugal"
identifier <- paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")

# Neuronlist save method
dbClass <- "ZIP"
zip <- ifelse(dbClass=="ZIP",TRUE,FALSE)

# Column types for reading individual per-query NBLAST CSVs
banc.col.types <- readr::cols(
  .default = readr::col_character(),
  cleft_segid  = readr::col_character(),
  centroid_x = readr::col_number(),
  centroid_y = readr::col_number(),
  centroid_z = readr::col_number(),
  bbox_bx = readr::col_number(),
  bbox_by = readr::col_number(),
  bbox_bz = readr::col_number(),
  bbox_ex = readr::col_number(),
  bbox_ey = readr::col_number(),
  bbox_ez = readr::col_number(),
  presyn_segid = readr::col_character(),
  postsyn_segid  = readr::col_character(),
  presyn_x = readr::col_integer(),
  presyn_y = readr::col_integer(),
  presyn_z = readr::col_integer(),
  postsyn_x = readr::col_integer(),
  postsyn_y = readr::col_integer(),
  postsyn_z = readr::col_integer(),
  clefthash = readr::col_number(),
  partnerhash = readr::col_integer(),
  size = readr::col_number(),
  l2_root = readr::col_number(),
  l2_nodes = readr::col_number(),
  l2_segments = readr::col_number(),
  l2_branchpoints = readr::col_number(),
  l2_endpoints = readr::col_number(),
  l2_cable_length = readr::col_number(),
  l2_n_trees = readr::col_number(),
  volume_nm3 = readr::col_number(),
  nb = readr::col_number(),
  score = readr::col_number(),
  hemibrain_nblast = readr::col_number(),
  fafb_nblast = readr::col_number(),
  manc_nblast = readr::col_number(),
  fanc_nblast = readr::col_number(),
  malecns_nblast = readr::col_number(),
  X = readr::col_number(),
  Y = readr::col_number(),
  Z = readr::col_number(),
  x = readr::col_number(),
  y = readr::col_number(),
  z = readr::col_number(),
  dcv_density = readr::col_number(),
  dcv_count = readr::col_integer(),
  histamine = readr::col_number() , 
  tyramine = readr::col_number() , 
  acetylcholine= readr::col_number() , 
  glutamate= readr::col_number() , 
  gaba= readr::col_number() , 
  glycine= readr::col_number() , 
  dopamine= readr::col_number() , 
  serotonin= readr::col_number() , 
  octopamine= readr::col_number() , 
  tyramine= readr::col_number() , 
  histamine= readr::col_number() , 
  nitric_oxide= readr::col_number() ,
  `allatostatin-a`= readr::col_number() , 
  `allatostatin-c`= readr::col_number() , 
  amnesiac= readr::col_number() , 
  bursicon= readr::col_number() , 
  capability= readr::col_number() , 
  ccap= readr::col_number() , 
  ccha1= readr::col_number() , 
  cnma= readr::col_number() , 
  corazonin= readr::col_number() , 
  darc1= readr::col_number() , 
  dh31= readr::col_number() , 
  dh331= readr::col_number() , 
  dh44= readr::col_number() , 
  dilp2= readr::col_number() , 
  dilp3= readr::col_number() , 
  dilp5= readr::col_number() , 
  dnpf= readr::col_number() , 
  drosulfakinin= readr::col_number() , 
  eclosion_hormone= readr::col_number() , 
  fmrf= readr::col_number() , 
  fmrfa= readr::col_number() , 
  hugin= readr::col_number() ,
  itp= readr::col_number() , 
  leucokinin= readr::col_number() , 
  mip= readr::col_number() , 
  myosuppressin= readr::col_number() , 
  myosupressin= readr::col_number() , 
  natalisin= readr::col_number() , 
  negative= readr::col_number() , 
  neuropeptide= readr::col_number() , 
  neuropeptides= readr::col_number() , 
  npf= readr::col_number() , 
  nplp1= readr::col_number() , 
  orcokinin= readr::col_number() , 
  pdf= readr::col_number() , 
  proctolin= readr::col_number() , 
  sifamide= readr::col_number() , 
  snpf= readr::col_number() , 
  space_blanket= readr::col_number() , 
  tachykinin= readr::col_number() , 
  trissin= readr::col_number(),
  input_connections = readr::col_number(),
  output_connections = readr::col_number(),
  total_connections = readr::col_number(),
  prepost = readr::col_number(),
  l2_cable_length_um  = readr::col_number(),
  l2_nodes  = readr::col_number(),
  input_side_index = readr::col_number(),
  output_side_index = readr::col_number(),
  pd_width = readr::col_number(),
  segregation_index = readr::col_number(),
  neurotransmitter_score = readr::col_number(),
  neurotransmitter_probability = readr::col_number(),
  neurotransmitter_predicted = readr::col_character(),
  count = readr::col_integer(),
  count_valid = readr::col_integer()
)

###################
##### REMOTES #####
###################

# Seatable
flytable_login()

# # slack
# slackr::slackr_setup()

# google drive 
hemibrainr_google_login <- file.path(hemibrain_data,"annotations")

# Rclone
rclone <- FALSE

#########################
##### ANNOUNCEMENTS #####
#########################

# Messages for debugging
message("R_MAX_VSIZE: ", Sys.getenv("R_MAX_VSIZE"))
message(".libPaths: ", print(.libPaths()))
message("bancr: ", packageVersion("bancr"))
print("##### SESSION INFO #####")
print(sessionInfo())
print("##### SESSION INFO #####")

# # Find Blender with python
# reticulate::py_run_string("import os")
# reticulate::py_run_string("from distutils.spawn import find_executable")
# reticulate::py_run_string("blender_executable <- find_executable('blender')")
# reticulate::py_run_string("print('Blender executable:', blender_executable)")

###############
### Colours ###
###############

banc.path <- "/n/data1/hms/neurobio/wilson/banc/BANC-project/"
if(file.exists("/n/data1/hms/neurobio/wilson/banc/BANC-project/settings/paper_colours_lacroix.csv")){
  paper.cols.df <- read.csv("/n/data1/hms/neurobio/wilson/banc/BANC-project/settings/paper_colours_lacroix.csv")
}else if(file.exists("setup/paper_colours.csv")){
  paper.cols.df <- read.csv("setup/paper_colours.csv")
}else{
  paper.cols.df <- read.csv("/Users/papers/BANC-project/settings/paper_colours_lacroix.csv")
}
paper.cols <- paper.cols.df$hex
names(paper.cols) <- paper.cols.df$label



