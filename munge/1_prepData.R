# clean up the data structure for use in downstream analysis

rm(list = ls())
"%ni%" = Negate("%in%")
library(dplyr)

proj_dir <- getwd()
data_dir <- paste0(proj_dir, "/data")
load(paste0(data_dir, "/illumina2007.rda"))

data <- illumina2007 %>%
  dplyr::select(which(colnames(.) %ni% c("box", "riboaveugml", "ng", "exclude",
                                         "hyb", "totalrnaug", "chip", "Chip.Id",
                                         "Chip.Section", "label.c", "benzene",
                                         "illumID", "berkeley_vial_label",
                                         "cRNA"))) %>%
  dplyr::filter(!duplicated(.$id)) %>%
  dplyr::mutate(
    benzene = I(newbenz),
    smoking = I(current_smoking)
  ) %>%
  dplyr::select(which(colnames(.) %ni% c("newbenz", "current_smoking")))

subjIDs <- data$id

rm(list = setdiff(ls(), c("data", "data_dir", "proj_dir", "%ni%", "subjIDs")))
