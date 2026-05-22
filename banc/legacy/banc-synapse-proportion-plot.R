###########################################################
### SUPERSEDED 2026-05-21 — use banc/metrics/banc-calculate-completion.R
### (section: "Proportion-on-cell ECDF") instead. That version is
### banc.version-aware and reads the 10 GB v2-enriched parquet
### (banc_meta.R output) rather than the raw 30 GB v821 CSV that
### OOM-killed at 64 G memory.
###
### Kept here for reference only — hardcoded to v821.
### -------------------------------------------------------
### Cumulative share of synapses in identified neurons   ###
###                                                       ###
### Produces banc_synapse_proportion_on_cell.{png,pdf} in ###
### BANC-project/figures/figure_1/links/supplement.       ###
###                                                       ###
### Standalone wrapper around the plot-producing portion  ###
### of banc/legacy/banc-assess-synapses.R — avoids the    ###
### ~30 min of post-plot analysis when the figure itself  ###
### is what's needed urgently.                            ###
###########################################################
source("banc/banc-startup.R")
library(vroom)

desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id',
                     'ctr_x', 'ctr_y', 'ctr_z')
col_types <- cols(
  id = col_character(), size = col_double(),
  pre_root_id = col_character(), post_root_id = col_character(),
  ctr_x = col_double(), ctr_y = col_double(), ctr_z = col_double(),
  .default = col_double()
)
column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                  'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                  'pre_root_id', 'post_supervoxel_id', 'post_root_id')

message("=== Querying SeaTable for proofread/identified neurons ===")
banc.meta <- banctable_query()
proof.svids    <- banc.meta %>% dplyr::filter(proofread == "TRUE")          %>% dplyr::pull(supervoxel_id)
roughly.svids  <- banc.meta %>% dplyr::filter(roughly_proofread == "TRUE")   %>% dplyr::pull(supervoxel_id)
bad.svids      <- banc.meta %>%
  dplyr::filter(grepl("not_a_neuron|trachea|merge|glia", super_class) |
                grepl("NOT_A_NEURON|MERGE|TRACHEA|GLIA", status)) %>%
  dplyr::pull(supervoxel_id)
other.svids    <- banc.meta %>% dplyr::filter(!supervoxel_id %in% bad.svids) %>% dplyr::pull(supervoxel_id)
svids          <- unique(c(proof.svids, roughly.svids, other.svids))

message(sprintf("  proof=%d  roughly=%d  bad=%d  total identified pool=%d",
                length(proof.svids), length(roughly.svids), length(bad.svids), length(svids)))

message("=== Resolving supervoxel_id -> root_id (version 821) ===")
identified.ids <- banc_rootid(svids,     version = "821")
bad.ids        <- banc_rootid(bad.svids, version = "821")
neuron.ids     <- c(bad.ids, identified.ids)

version.path <- file.path(banc.save.path, "v821")
csv.path     <- file.path(version.path, "synapses_v2_human_readable.csv")

message(sprintf("=== Reading v821 synapse CSV (%.1f GB) ===",
                file.info(csv.path)$size / 1024^3))
t0 <- Sys.time()
banc.syns.review <- vroom::vroom(csv.path,
                                  col_names = column_names,
                                  col_select = dplyr::all_of(desired_columns),
                                  col_types = col_types,
                                  skip = 1) %>%
  dplyr::rename(X = ctr_x, Y = ctr_y, Z = ctr_z) %>%
  dplyr::mutate(pre_status  = dplyr::if_else(pre_root_id  %in% !!neuron.ids, "neuron", "fragment"),
                post_status = dplyr::if_else(post_root_id %in% !!neuron.ids, "neuron", "fragment")) %>%
  tibble::as_tibble()
message(sprintf("  read %s rows in %.1f min",
                format(nrow(banc.syns.review), big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

process_root_id <- function(df, root_col, status_col, label) {
  df %>%
    dplyr::group_by(root_id = .data[[root_col]], pre_status = .data[[status_col]]) %>%
    dplyr::mutate(n_syn = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pre_status, n_syn) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(n_syn) %>%
    tidyr::pivot_wider(names_from = pre_status, values_from = n, values_fill = 0) %>%
    dplyr::mutate(
      cum_neuron   = if ("neuron"   %in% names(.)) base::cumsum(neuron)   else rep(0, dplyr::n()),
      cum_fragment = if ("fragment" %in% names(.)) base::cumsum(fragment) else rep(0, dplyr::n()),
      cum_total    = cum_neuron + cum_fragment,
      pct_neuron   = dplyr::if_else(cum_total > 0, cum_neuron / cum_total, NA_real_),
      root_id_col  = label
    )
}

message("=== Building cumulative distributions ===")
dat_pre  <- process_root_id(banc.syns.review, "pre_root_id",  "pre_status",  "pre_root_id")
dat_post <- process_root_id(banc.syns.review, "post_root_id", "post_status", "post_root_id")
syn_comp <- dplyr::bind_rows(dat_pre, dat_post)

g.ecdf <- ggplot2::ggplot(syn_comp,
                          ggplot2::aes(x = n_syn, y = pct_neuron, color = root_id_col)) +
  ggplot2::geom_line(size = 1.2) +
  ggplot2::scale_x_log10() +
  ggplot2::annotation_logticks(sides = "b") +
  ggplot2::labs(x = "number of synapses per root_id (threshold)",
                y = "proportion of synapses in neurons (≤ threshold)",
                color = "",
                title = "cumulative share of synapses in neurons by root type") +
  ggplot2::scale_color_manual(values = c(pre_root_id  = paper.cols[["pre"]],
                                          post_root_id = paper.cols[["post"]])) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "none")

banc.fig1.supp.path <- "/n/data1/hms/neurobio/wilson/banc/BANC-project/figures/figure_1/links/supplement"
dir.create(banc.fig1.supp.path, showWarnings = FALSE, recursive = TRUE)
for (.ext in c("png", "pdf")) {
  out <- file.path(banc.fig1.supp.path,
                   sprintf("banc_synapse_proportion_on_cell.%s", .ext))
  ggplot2::ggsave(plot = g.ecdf, filename = out,
                  width = 6, height = 4, dpi = 300, bg = "transparent")
  message(sprintf("wrote %s", out))
}

message("=== done ===")
