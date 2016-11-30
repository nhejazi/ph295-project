# Example preprocessing script.

"%ni%" = Negate("%in%")

proj_dir <- getwd()
data_dir <- paste0(proj_dir, "/data")
load(paste0(data_dir, "/illumina2007.rda"))

library(dplyr)
data <- illumina2007 %>%
  select(which(colnames(.) %ni% c("box", "riboaveugml", "ng",
                                  "totalrnaug", "berkeley_vial_label",
                                  "chip", "Chip.Id", "Chip.Section",
                                  "cRNA", "label.c"))) %>%
  dplyr::filter(!duplicated(.$id))

