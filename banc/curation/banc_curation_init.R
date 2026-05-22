# banc/curation init: shared helpers + library loads for the curation
# scripts in this directory (banc_meta_*, banc_cell_representative_point,
# make_codex_annotations_flat_table_v888). source() this file from any
# curation script with:
#     source("banc/curation/banc_curation_init.R")
# Originally authored by Helen Yang.

# load libraries
library(bancr)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(bit64)
library(rlang)
library(googlesheets4)

# path to txt file tracking banc_meta updates from cell_info
banc_meta_cell_info_update_file <- "./annotations/banc_meta_from_cell_info_update_log.txt"

cell_ids_tag_filepath <- "./annotations/cell_ids_tag.txt"

# load helper functions
is_na_or_empty <- function(x) {
  return(is.na(x) | x == "" | (is.character(x) & trimws(x) == ""))
}

is_not_empty_na <- function(x) {
  return(!is.na(x) & x != "" & (is.character(x) & trimws(x) != ""))
}

append_status <- function(status, update){
  # Handle NA cases
  if(is.na(update) && is.na(status)) {
    return(NA_character_)
  } else if(is.na(update)) {
    # Only update is NA, process status alone
    update.col <- paste(sort(unique(unlist(strsplit(status, split=",|, ")))), collapse=",")
    return(gsub("^,| ","", update.col))
  } else if(is.na(status)) {
    # Only status is NA, process update alone
    update.col <- paste(sort(unique(unlist(strsplit(update, split=",|, ")))), collapse=",")
    return(gsub("^,| ","", update.col))
  }
  
  # Neither is NA - original logic
  status <- paste(c(status, update), collapse=",")
  update.col <- paste(sort(unique(unlist(strsplit(status, split=",|, ")))), collapse=",")
  gsub("^,| ","", update.col)
}

subtract_status <- function(status, update, invert = FALSE){
  satuses <- sort(unique(unlist(strsplit(status,split=",|, "))))
  if(invert){
    satuses <- sort(unique(intersect(satuses,update)))
  }else{
    satuses <- sort(unique(setdiff(satuses,update)))
  }
  update.col<-paste0(satuses,collapse=",")
  gsub("^,| ","",update.col)
}

sort_status <- function(status) {
  return(paste0(sort(unique(unlist(strsplit(status,split=",|, ")))),collapse=","))
}

# function to standardize position
standardize_position <- function(position_vector) {
  result <- ifelse(
    is.na(position_vector) | position_vector == "",
    position_vector,
    str_replace_all(str_trim(position_vector), "\\s*,\\s*", ", ")
  )
  
  return(result)
}

# function to join other_names - handles spaces
join_other_names <- function(status, update){
  # Helper function to split and clean items
  split_and_clean <- function(text) {
    if(is.na(text) || text == "") {
      return(character(0))
    }
    # Split on comma with or without space
    items <- unlist(strsplit(text, split=",\\s*"))
    # Trim only leading/trailing whitespace, preserve internal spaces
    items <- trimws(items)
    # Remove empty items
    items <- items[items != ""]
    return(items)
  }
  
  # Handle NA cases
  if(is.na(update) && is.na(status)) {
    return(NA_character_)
  }
  
  # Get items from both inputs
  status_items <- split_and_clean(status)
  update_items <- split_and_clean(update)
  
  # Combine all items
  all_items <- c(status_items, update_items)
  
  # If no items, return NA
  if(length(all_items) == 0) {
    return(NA_character_)
  }
  
  # Get unique items (preserving internal spaces) and sort
  unique_items <- unique(all_items)
  sorted_items <- sort(unique_items)
  
  # Join with comma and space
  return(paste(sorted_items, collapse=", "))
}

# function to remove items from other_names - handles spaces
remove_other_names <- function(status, update){
  # Helper function to split and clean items
  split_and_clean <- function(text) {
    if(is.na(text) || text == "") {
      return(character(0))
    }
    # Split on comma with or without space
    items <- unlist(strsplit(text, split=",\\s*"))
    # Trim only leading/trailing whitespace, preserve internal spaces
    items <- trimws(items)
    # Remove empty items
    items <- items[items != ""]
    return(items)
  }
  
  # Handle NA cases
  if(is.na(status)) {
    return(NA_character_)
  }
  
  if(is.na(update)) {
    return(status)  # Nothing to remove, return original status
  }
  
  # Get items from both inputs
  status_items <- split_and_clean(status)
  update_items <- split_and_clean(update)
  
  # Remove update_items from status_items
  remaining_items <- setdiff(status_items, update_items)
  
  # If no items remain, return NA
  if(length(remaining_items) == 0) {
    return(NA_character_)
  }
  
  # Sort remaining items
  sorted_items <- sort(remaining_items)
  
  # Join with comma and space
  return(paste(sorted_items, collapse=", "))
}
