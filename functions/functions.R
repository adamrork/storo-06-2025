#################
### FUNCTIONS ###
#################

# Define a helpful contrast renaming function #
renameContrast <- function(name, addition = NULL, flip = FALSE) {

  # This is highly specific to these particular contrasts and cannot be widely reused as-is. #
  renamed <- gsub("_D", "_Day ", name) %>%
    gsub("_", " ", .) %>%
    gsub("vs", "vs.", .) %>%
    gsub("PRG4.VEGF", "PRG4-VEGF", .) %>%
    gsub("Day 00", "", .)

  if ( flip == TRUE ) {
    parts <- strsplit(renamed, " vs\\. ")[[1]]
    renamed <- paste0(parts[[2]], " vs. ", parts[[1]])
  }

  if ( is.null(addition) ) {
    return(renamed)

  } else {
    return(paste0(renamed, " ", addition))

  }
}
