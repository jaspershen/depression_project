library(tidyverse)
library(ggplot2)

convert2meters <- function(height) {
  if (is.na(height))
    return(NA)
  parts <- strsplit(height, "'")[[1]]
  feet <- as.numeric(parts[1])
  inches <- as.numeric(parts[2])
  total_inches <- feet * 12 + inches
  total_meters <- total_inches * 0.0254
  return(total_meters)
}