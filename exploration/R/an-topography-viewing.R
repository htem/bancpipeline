### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

### take code from neck_connective_neuropils_an.R

t1r <- subset(combined_df, neuropil=="v_LNp_T1_R")
t2r <- subset(combined_df, neuropil=="v_LNp_T2_R")
t3r <- subset(combined_df, neuropil=="v_LNp_T3_R")

library(dplyr)

# Assuming your dataframes are named df1, df2, and df3

result <- bind_rows(
  t1r %>% mutate(source = "t1r"),
  t2r %>% mutate(source = "t2r"),
  t3r %>% mutate(source = "t3r")
) %>%
  group_by(neuron) %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  split(.$source) %>%
  lapply(function(df) df %>% select(-source))

t1r <- result$t1r
t2r <- result$t2r
t3r <- result$t3r

# set threshold
t1r_thresh <- subset(t1r, count>100)
t2r_thresh <- subset(t2r, count>100)
t3r_thresh <- subset(t3r, count>100)

#get neurons back
t1_neuron <- l2[t1r_thresh$neuron]
t2_neuron <- l2[t2r_thresh$neuron]
t3_neuron <- l2[t3r_thresh$neuron]

# plot T1 T2 and T3 neurons in different colors
clear3d()
#plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(t1_neuron, col="red")
plot3d(t2_neuron, col="blue")
plot3d(t3_neuron, col = "green")

#3d kernel density plot showing each kernel that skews towards

#alphashape3d

#try now with inbetween dfs

library(dplyr)
library(purrr)

process_dataframes <- function(df1, df2, df3) {
  combined <- bind_rows(
    t1r %>% mutate(source = "t1r"),
    t2r %>% mutate(source = "t2r"),
    t3r %>% mutate(source = "t3r")
  )
  
  result <- combined %>%
    group_by(neuron) %>%
    mutate(
      max_count = max(count),
      count_diff = (max_count - count) / max_count,
      in_multiple = sum(count_diff <= 0.5) > 1
    ) %>%
    mutate(
      new_source = case_when(
        in_multiple & any(source %in% c("t1r", "t2r") & count_diff <= 0.5) ~ "t12r",
        in_multiple & any(source %in% c("t1r", "t3r") & count_diff <= 0.5) ~ "t13r",
        in_multiple & any(source %in% c("t2r", "t3r") & count_diff <= 0.5) ~ "t23r",
        TRUE ~ source
      )
    ) %>%
    ungroup() %>%
    select(-max_count, -count_diff, -in_multiple)
  
  # Split into separate dataframes
  dfs <- result %>%
    split(.$new_source) %>%
    map(~ .x %>% select(-source, -new_source))
  
  # Ensure all six dataframes exist, even if empty
  all_dfs <- c("t1r", "t2r", "t3r", "t12r", "t13r", "t23r")
  dfs[setdiff(all_dfs, names(dfs))] <- list(tibble())
  
  return(dfs)
}

# Use the function
results <- process_dataframes(t1r, t2r, t3r)

# Now you can access the processed dataframes like this:
t1r_processed <- results$t1r
t2r_processed <- results$t2r
t3r_processed <- results$t3r


t1r_thresh <- subset(t1r_processed, count>15)
t2r_thresh <- subset(t2r_processed, count>15)
t3r_thresh <- subset(t3r_processed, count>15)


#get neurons back
t1_neuron <- l2[t1r_thresh$neuron]
t2_neuron <- l2[t2r_thresh$neuron]
t3_neuron <- l2[t3r_thresh$neuron]

# plot T1 T2 and T3 neurons in different colors
clear3d()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(t1_neuron, col="red")
plot3d(t2_neuron, col="blue")
plot3d(t3_neuron, col = "green")