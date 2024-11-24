#################################################################################################################.
#################################################################################################################.
###
###   Author: Jose Camacho Valenzuela.
###   Date: XX.
###   R version: R version 4.1.3 (2022-03-10) -- "One Push-Up"
###   Script: Determine CNV-LOH changes and integrate multidimensional data (germline, tumor, clinical).
###
###               1) This script is executed in one individual sample of the breast cohort.
###               2) Determine LOH in every germline variant identified.
###               3) Determine somatic full loss via CNV in the 6 HRD genes:
###                   - BRCA1/2, RAD51C/D, BARD1, and PALB2.
###               6) Determine somatic point mutations as second hit.
###               7) Integrate multidimensional data:
###                   - Germline data:
###                       - Info of all the variants identified in the patient in question.
###                   - Tumor data:
###                       - Second hit: LOH status of the germline variant via CNVs.
###                       - Second hit: somatic point mutations.
###                       - Somatic full loss of both copies of the 6 HRD genes.
###                       - Methylation in BRCA1 and RAD51C.
###                       - Signature 3 (SigMA) scores.
###                   - Clinical data:
###                       - Age, Sex, Subtype, among others.
###
###   Data source: TCGA-BRCA.
###   Sample: TCGA-XX
###
#################################################################################################################.
#################################################################################################################.

# 1. Load the required libraries.
library(data.table)
library(tidyr)
library(dplyr)
library(readxl)
library(openxlsx)
library(tidyverse)

# 2. Import the source ABSOLUTE CNV data of the breast cohort.
setwd("C:/Users/path/to/dir")
ABSOLUTE <- fread(file = "ABSOLUTE_cnv_data_cohort.txt") %>% setDT()

# 3. Import the master dataframe of the germline exome cohort with Sig3+.
setwd("C:/Users/path/to/dir")
mdf.germ <- fread(file = "masterdf_germline_exomes_cohort.txt") %>% setDT()

# 4. Subset the CNV data of only the Sig3+ group samples.
# ----- a) Adjust the TCGA IDs to match the CNV data IDs.
mdf.germ2 <- mdf.germ
mdf.germ2[, CNV_Sample_Barcode := str_extract(
  Tumor_Sample_Barcode,"TCGA-..-....-..")]

# ----- b) Do the subset.
ABSOLUTE.sig3 <- ABSOLUTE[Sample %in% mdf.germ2$CNV_Sample_Barcode]
unique(ABSOLUTE.sig3$Sample)

# ----- c) Rename the "NAs" in the merged object.
ABSOLUTE.sig3[is.na(ABSOLUTE.sig3)] <- "NA"

# 5. Get only the columns of interest of the resulting ABSOLUTE dataframe.
ABSOLUTE.sig3.col <- select(ABSOLUTE.sig3, 1:9, 18:20)

# 6. Merge the IDs columns to the resulting ABSOLUTE df.
mdf.germ3 <- mdf.germ2
mdf.germ3 <- select(mdf.germ3, 1:2, 51)
mdf.germ3.nodup <- distinct(mdf.germ3,
                            CNV_Sample_Barcode,
                            .keep_all = TRUE)

setnames(ABSOLUTE.sig3.col, "Sample", "CNV_Sample_Barcode")
ABSOLUTE.merged <- merge(ABSOLUTE.sig3.col, mdf.germ3.nodup,
                         by = "CNV_Sample_Barcode",
                         all.x = TRUE, all.y = TRUE)

ABSOLUTE.merged2 <- select(ABSOLUTE.merged, 13:14, 1:11)
ABSOLUTE.merged2[is.na(ABSOLUTE.merged2)] <- "NA"

# 7. Remove the samples with "NA"s (empty CNV info from ABSOLUTE).
samples.to.be.removed <- ABSOLUTE.merged2[Chromosome == "NA"]
ABSOLUTE.merged3 <- ABSOLUTE.merged2[!Chromosome == "NA"]
unique(ABSOLUTE.merged3$Tumor_Sample_Barcode)

# 8. Get only the sample of interest and curate it as required for the CNV integration.
ABSOLUTE.sample <- ABSOLUTE.merged3[Matched_Norm_Sample_Barcode == "TCGA-XX"]
setnames(ABSOLUTE.sample, "Chromosome", "Chr")
setnames(ABSOLUTE.sample, "Start", "Start.cnv")
setnames(ABSOLUTE.sample, "End", "End.cnv")
ABSOLUTE.sample2 <- ABSOLUTE.sample
ABSOLUTE.sample2$Chr <- paste("chr", ABSOLUTE.sample2$Chr, sep = "")

# 9. Get only the sample of interest of the Sig3+ group of the masterdf germline Sig3+ file.
snp <- mdf.germ2[Matched_Norm_Sample_Barcode == "TCGA-XX"]

# 10. Backup the resulting objects, convert as.numeric the Start.cnv, End.cnv and Start columns and rename headers appropriately.
cnv.abs <- ABSOLUTE.sample2
snp.abs <- snp
cnv.abs$Start.cnv <- as.numeric(cnv.abs$Start.cnv)
cnv.abs$End.cnv <- as.numeric(cnv.abs$End.cnv)
snp.abs$Start <- as.numeric(snp.abs$Start)
setnames(cnv.abs, "Modal_Total_CN", "CNt")
setnames(cnv.abs, "Modal_HSCN_1", "B_Modal_HSCN_1")
setnames(cnv.abs, "Modal_HSCN_2", "A_Modal_HSCN_2")


# ----------------------------------------------------------------------------------------------.
################################################################################################.
##### ------------------------------------------------------------------------------------ #####.
#####                  Determine second hit (LOH) via CNVs                                 #####.
##### ------------------------------------------------------------------------------------ #####.
################################################################################################.

# 11. Determine whether a variant is within a cnv segment based on variants' genomic coordinates.
# ----- a) Function to show if the genomic coordinate of a germline variant falls within the cnv in question.
is_between <- function(Start, Start.cnv, End.cnv) {
  Start >= Start.cnv & Start <= End.cnv
}

final.abs.df <- snp.abs %>%
  inner_join(cnv.abs, by = "Chr") %>% # Merge / concatenate both cnv and snp by "Chr" column
  group_by(Chr) %>% # Group/order by "Chr" column
  mutate(Test_variant_within_segment = is_between(Start, Start.cnv, End.cnv)) %>% # Indicate if the variant is within a cnv segment via is_between function
  filter(Test_variant_within_segment == "TRUE") %>%  # Select only the boolean TRUE (meaning the variant is within a cnv based on its genomic coordinates)
  ungroup()

# 12. Get only the columns of interest and rename the required headers appropriately.
final.abs.df.col <- select(final.abs.df, 1:2, 51, 3:6, 8:16, 26:28, 42, 55:64)
colnames(final.abs.df.col)
setnames(final.abs.df.col, "Tumor_Sample_Barcode.x", "Tumor_Sample_Barcode")
setnames(final.abs.df.col, "Matched_Norm_Sample_Barcode.x", "Matched_Norm_Sample_Barcode")
setnames(final.abs.df.col, "CNV_Sample_Barcode.x", "CNV_Sample_Barcode")

final.abs.df.col$Germline_data <- "Germline"
final.abs.df.col$Somatic_data <- "Somatic"

# 13. Determine the LOH status of every germline variant in the output xlsx of the WES algorithm.
final.df.abs.LOH <- final.abs.df.col %>%
  mutate(LOH_status = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-LOH", # Hemizygous LOH
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-LOH", # Copy neutral LOH
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-LOH", # Duplication LOH
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous", # Homozygous deletion
  ))

# 14. Reorder the columns appropriately.
final.df.abs.LOH.col <- select(final.df.abs.LOH, 1:3, 31, 5:20, 32,
                               21:24, 27, 26, 25, 28, 33, 29, 4)

final.df.LOH.accumulated <- final.df.abs.LOH.col


# ------------------------------------------------------------------------------------------------------.
########################################################################################################.
##### -------------------------------------------------------------------------------------------- #####.
#####               Determine loss of both copies of the 6 HRD genes via CNVs                      #####.
#####                     BARD1, BRCA1, BRCA2, RAD51C, RAD51D, PALB2                               #####.
##### -------------------------------------------------------------------------------------------- #####.
########################################################################################################.

# ---------------------------------------------------------------------------------------------.
# ----------                          BRCA2  -  Chr13                                ----------.
# ---------------------------------------------------------------------------------------------.

# 15. Select only the segments/CNVs in chromosome 13 for BRCA2.
HRDness <- cnv.abs %>% setDT()
BRCA2 <- HRDness[Chr == "chr13"]

# 16. Add as columns the start and end genomic coordinates of BRCA2.
BRCA2$Start.BRCA2 <- 32889645
BRCA2$End.BRCA2 <- 32974405

# 17. Use a function to determine if the Start position of BRCA2 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.BRCA2, End.BRCA2) {
  Start.BRCA2 >= Start.cnv & Start.BRCA2 <= End.cnv
}

df <- BRCA2 %>%
  mutate(Test_BRCA2_start = is_between(Start.cnv, End.cnv, Start.BRCA2, End.BRCA2))

# 18. Use a function to determine if the End position of BRCA2 falls within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.BRCA2, End.BRCA2) {
  End.BRCA2 >= Start.cnv & End.BRCA2 <= End.cnv
}

df2 <- df %>%
  mutate(Test_BRCA2_end = is_between(Start.cnv, End.cnv, Start.BRCA2, End.BRCA2))

df3 <- df2
df3$Test_BRCA2_start <- as.character(df3$Test_BRCA2_start)
df3$Test_BRCA2_end <- as.character(df3$Test_BRCA2_end)

# 19. Determine if the BRCA2 start and end coordinates together fall within a segment/cnv.
# ----- a) Write "1" if the BRCA2 start coordinate is within a segment, and "0" if not.
df4 <- df3 %>%
  mutate(BRCA2_start_boolean1 = case_when(
    startsWith (Test_BRCA2_start, "TRUE") ~ 1,
    startsWith (Test_BRCA2_start, "FALSE") ~ 0
  ))

# ----- b) Write "1" if the BRCA2 end coordinate is within a segment, and "0" if not.
df5 <- df4 %>%
  mutate(BRCA2_start_boolean2 = case_when(
    startsWith (Test_BRCA2_end, "TRUE") ~ 1,
    startsWith (Test_BRCA2_end, "FALSE") ~ 0
  ))

# ----- c) Sum the numbers of BRCA2_start_boolean1 and BRCA2_start_boolean2.
df6 <- df5
df6$sum_BRCA2 <- df6$BRCA2_start_boolean1 + df6$BRCA2_start_boolean2

# ----- d) In another column:
# Write "BRCA2_fully_present" if both BRCA2 coordinates fall within a segment.
# Write "BRCA2_partially_present" if only the start or end BRCA2 coordinates fall within a segment.
# Write "None" if both BRCA2 coordinates DO NOT fall within a segment.
df7 <- df6 %>%
  mutate(BRCA2_presence = case_when(
    sum_BRCA2 == 0 ~ "None",
    sum_BRCA2 == 1 ~ "BRCA2_partially_present",
    sum_BRCA2 == 2 ~ "BRCA2_fully_present",
  ))

# 20. Based on the integer copy number values, determine LOH for each segment/cnv.
df8 <- df7 %>%
  mutate(LOH_status_BRCA2 = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-CNV",
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-CNV",
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-CNV",
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous",
  ))

# 21. Get only the desired columns.
df9 <- select(df8, 1:13, 21:22)

# 22. Select only the row/rows where BRCA2 is present (if present in a segment).
df10 <- df9[BRCA2_presence %in% c("BRCA2_partially_present",
                                  "BRCA2_fully_present")]

# 23. Add as an extra column to the main dataframe the BRCA2 LOH status.
final.df.LOH.BRCA2 <- final.df.LOH.accumulated
final.df.LOH.BRCA2$BRCA2_LOH_status <- df10$LOH_status_BRCA2
final.df.LOH.BRCA2$BRCA2_in_segment_status <- df10$BRCA2_presence
final.df.LOH.accumulated <- final.df.LOH.BRCA2


# ---------------------------------------------------------------------------------------------.
# ----------                          BRCA1  -  Chr17                                ----------.
# ---------------------------------------------------------------------------------------------.

# 24. Select only the segments/CNVs in chromosome 17 for BRCA1.
HRDness <- cnv.abs %>% setDT()
BRCA1 <- HRDness[Chr == "chr17"]

# 25. Add as columns the start and end genomic coordinates of BRCA1.
BRCA1$Start.BRCA1 <- 41196312
BRCA1$End.BRCA1 <- 41277381

# 26. Use a function to determine if the Start position of BRCA1 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.BRCA1, End.BRCA1) {
  Start.BRCA1 >= Start.cnv & Start.BRCA1 <= End.cnv
}

df <- BRCA1 %>%
  mutate(Test_BRCA1_start = is_between(Start.cnv, End.cnv, Start.BRCA1, End.BRCA1))

# 27. Use a function to determine if the End position of BRCA1 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.BRCA1, End.BRCA1) {
  End.BRCA1 >= Start.cnv & End.BRCA1 <= End.cnv
}

df2 <- df %>%
  mutate(Test_BRCA1_end = is_between(Start.cnv, End.cnv, Start.BRCA1, End.BRCA1))

# 28. Convert to as.character the resulting columns.
df3 <- df2
df3$Test_BRCA1_start <- as.character(df3$Test_BRCA1_start)
df3$Test_BRCA1_end <- as.character(df3$Test_BRCA1_end)

# 29. Determine if the BRCA1 start and end coordinates together fall within a segment/cnv.
# ----- a) Write "1" if the BRCA1 start coordinate is within a segment, and "0" if not.
df4 <- df3 %>%
  mutate(BRCA1_start_boolean1 = case_when(
    startsWith (Test_BRCA1_start, "TRUE") ~ 1,
    startsWith (Test_BRCA1_start, "FALSE") ~ 0
  ))

# ----- b) Write "1" if the BRCA1 end coordinate is within a segment, and "0" if not.
df5 <- df4 %>%
  mutate(BRCA1_start_boolean2 = case_when(
    startsWith (Test_BRCA1_end, "TRUE") ~ 1,
    startsWith (Test_BRCA1_end, "FALSE") ~ 0
  ))

# ----- c) Sum the numbers of BRCA1_start_boolean1 and BRCA1_start_boolean2.
df6 <- df5
df6$sum_BRCA1 <- df6$BRCA1_start_boolean1 + df6$BRCA1_start_boolean2

# ----- d) In another column:
# Write "BRCA1_fully_present" if both BRCA1 coordinates fall within a segment.
# Write "BRCA1_partially_present" if only the start or end BRCA1 coordinates fall within a segment.
# Write "None" if both BRCA1 coordinates DO NOT fall within a segment.
df7 <- df6 %>%
  mutate(BRCA1_presence = case_when(
    sum_BRCA1 == 0 ~ "None",
    sum_BRCA1 == 1 ~ "BRCA1_partially_present",
    sum_BRCA1 == 2 ~ "BRCA1_fully_present",
  ))

# 30. Based on the integer copy number values, determine LOH for each segment/cnv.
df8 <- df7 %>%
  mutate(LOH_status_BRCA1 = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-CNV",
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-CNV",
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-CNV",
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous",
  ))

# 31. Get only the desired columns.
df9 <- select(df8, 1:13, 21:22)

# 32. Select only the row/rows where BRCA1 is present (if present in a segment).
df10 <- df9[BRCA1_presence %in% c("BRCA1_partially_present",
                                  "BRCA1_fully_present")]

# 33. Add as an extra column to the main dataframe the BRCA1 LOH status.
final.df.LOH.BRCA1 <- final.df.LOH.accumulated
final.df.LOH.BRCA1$BRCA1_LOH_status <- df10$LOH_status_BRCA1
final.df.LOH.BRCA1$BRCA1_in_segment_status <- df10$BRCA1_presence
final.df.LOH.accumulated <- final.df.LOH.BRCA1


# ---------------------------------------------------------------------------------------------.
# ----------                          BARD1  -  Chr2                                 ----------.
# ---------------------------------------------------------------------------------------------.

# 34. Select only the segments/CNVs in chromosome 2 for BARD1.
HRDness <- cnv.abs %>% setDT()
BARD1 <- HRDness[Chr == "chr2"]

# 35. Add as columns the start and end genomic coordinates of BARD1.
BARD1$Start.BARD1 <- 215590370
BARD1$End.BARD1 <- 215674407

# 36. Use a function to determine if the Start position of BARD1 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.BARD1, End.BARD1) {
  Start.BARD1 >= Start.cnv & Start.BARD1 <= End.cnv
}

df <- BARD1 %>%
  mutate(Test_BARD1_start = is_between(Start.cnv, End.cnv, Start.BARD1, End.BARD1))

# 37. Use a function to determine if the End position of BARD1 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.BARD1, End.BARD1) {
  End.BARD1 >= Start.cnv & End.BARD1 <= End.cnv
}

df2 <- df %>%
  mutate(Test_BARD1_end = is_between(Start.cnv, End.cnv, Start.BARD1, End.BARD1))

# 38. Convert to as.character the resulting columns.
df3 <- df2
df3$Test_BARD1_start <- as.character(df3$Test_BARD1_start)
df3$Test_BARD1_end <- as.character(df3$Test_BARD1_end)

# 39. Determine if the BARD1 start and end coordinates together fall within a segment/cnv.
# ----- a) Write "1" if the BARD1 start coordinate is within a segment, and "0" if not.
df4 <- df3 %>%
  mutate(BARD1_start_boolean1 = case_when(
    startsWith (Test_BARD1_start, "TRUE") ~ 1,
    startsWith (Test_BARD1_start, "FALSE") ~ 0
  ))

# ----- b) Write "1" if the BARD1 end coordinate is within a segment, and "0" if not.
df5 <- df4 %>%
  mutate(BARD1_start_boolean2 = case_when(
    startsWith (Test_BARD1_end, "TRUE") ~ 1,
    startsWith (Test_BARD1_end, "FALSE") ~ 0
  ))

# ----- c) Sum the numbers of BARD1_start_boolean1 and BARD1_start_boolean2.
df6 <- df5
df6$sum_BARD1 <- df6$BARD1_start_boolean1 + df6$BARD1_start_boolean2

# ----- d) In another column:
# Write "BARD1_fully_present" if both BARD1 coordinates fall within a segment.
# Write "BARD1_partially_present" if only the start or end BARD1 coordinates fall within a segment.
# Write "None" if both BARD1 coordinates DO NOT fall within a segment.
df7 <- df6 %>%
  mutate(BARD1_presence = case_when(
    sum_BARD1 == 0 ~ "None",
    sum_BARD1 == 1 ~ "BARD1_partially_present",
    sum_BARD1 == 2 ~ "BARD1_fully_present",
  ))

# 40. Based on the integer copy number values, determine LOH for each segment/cnv.
df8 <- df7 %>%
  mutate(LOH_status_BARD1 = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-CNV",
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-CNV",
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-CNV",
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous",
  ))

# 41. Get only the desired columns.
df9 <- select(df8, 1:13, 21:22)

# 42. Select only the row/rows where BARD1 is present (if present in a segment).
df10 <- df9[BARD1_presence %in% c("BARD1_partially_present",
                                  "BARD1_fully_present")]

# 43. Add as an extra column to the main dataframe the BARD1 LOH status.
final.df.LOH.BARD1 <- final.df.LOH.accumulated
final.df.LOH.BARD1$BARD1_LOH_status <- df10$LOH_status_BARD1
final.df.LOH.BARD1$BARD1_in_segment_status <- df10$BARD1_presence
final.df.LOH.accumulated <- final.df.LOH.BARD1


# ---------------------------------------------------------------------------------------------.
# ----------                         RAD51C  -  Chr17                                ----------.
# ---------------------------------------------------------------------------------------------.

# 44. Select only the segments/CNVs in chromosome 17 for RAD51C.
HRDness <- cnv.abs %>% setDT()
RAD51C <- HRDness[Chr == "chr17"]

# 45. Add as columns the start and end genomic coordinates of RAD51C.
RAD51C$Start.RAD51C <- 56769934
RAD51C$End.RAD51C <- 56812972

# 46. Use a function to determine if the Start position of RAD51C is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.RAD51C, End.RAD51C) {
  Start.RAD51C >= Start.cnv & Start.RAD51C <= End.cnv
}

df <- RAD51C %>%
  mutate(Test_RAD51C_start = is_between(Start.cnv, End.cnv, Start.RAD51C, End.RAD51C))

# 47. Use a function to determine if the End position of RAD51C is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.RAD51C, End.RAD51C) {
  End.RAD51C >= Start.cnv & End.RAD51C <= End.cnv
}

df2 <- df %>%
  mutate(Test_RAD51C_end = is_between(Start.cnv, End.cnv, Start.RAD51C, End.RAD51C))

# 48. Convert to as.character the resulting columns.
df3 <- df2
df3$Test_RAD51C_start <- as.character(df3$Test_RAD51C_start)
df3$Test_RAD51C_end <- as.character(df3$Test_RAD51C_end)

# 49. Determine if the RAD51C start and end coordinates together fall within a segment/cnv.
# ----- a) Write "1" if the RAD51C start coordinate is within a segment, and "0" if not.
df4 <- df3 %>%
  mutate(RAD51C_start_boolean1 = case_when(
    startsWith (Test_RAD51C_start, "TRUE") ~ 1,
    startsWith (Test_RAD51C_start, "FALSE") ~ 0
  ))

# ----- b) Write "1" if the RAD51C end coordinate is within a segment, and "0" if not.
df5 <- df4 %>%
  mutate(RAD51C_start_boolean2 = case_when(
    startsWith (Test_RAD51C_end, "TRUE") ~ 1,
    startsWith (Test_RAD51C_end, "FALSE") ~ 0
  ))

# ----- c) Sum the numbers of RAD51C_start_boolean1 and RAD51C_start_boolean2.
df6 <- df5
df6$sum_RAD51C <- df6$RAD51C_start_boolean1 + df6$RAD51C_start_boolean2

# ----- d) In another column:
# Write "RAD51C_fully_present" if both RAD51C coordinates fall within a segment.
# Write "RAD51C_partially_present" if only the start or end RAD51C coordinates fall within a segment.
# Write "None" if both RAD51C coordinates DO NOT fall within a segment.
df7 <- df6 %>%
  mutate(RAD51C_presence = case_when(
    sum_RAD51C == 0 ~ "None",
    sum_RAD51C == 1 ~ "RAD51C_partially_present",
    sum_RAD51C == 2 ~ "RAD51C_fully_present",
  ))

# 50. Based on the integer copy number values, determine LOH for each segment/cnv.
df8 <- df7 %>%
  mutate(LOH_status_RAD51C = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-CNV",
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-CNV",
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-CNV",
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous",
  ))

# 51. Get only the desired columns.
df9 <- select(df8, 1:13, 21:22)

# 52. Select only the row/rows where RAD51C is present (if present in a segment).
df10 <- df9[RAD51C_presence %in% c("RAD51C_partially_present",
                                   "RAD51C_fully_present")]

# 53. Add as an extra column to the main dataframe the RAD51C LOH status.
final.df.LOH.RAD51C <- final.df.LOH.accumulated
final.df.LOH.RAD51C$RAD51C_LOH_status <- df10$LOH_status_RAD51C
final.df.LOH.RAD51C$RAD51C_in_segment_status <- df10$RAD51C_presence
final.df.LOH.accumulated <- final.df.LOH.RAD51C


# ---------------------------------------------------------------------------------------------.
# ----------                          RAD51D  -  Chr17                                ----------.
# ---------------------------------------------------------------------------------------------.

# 54. Select only the segments/CNVs in chromosome 17 for RAD51D.
HRDness <- cnv.abs %>% setDT()
RAD51D <- HRDness[Chr == "chr17"]

# 55. Add as columns the start and end genomic coordinates of RAD51D.
RAD51D$Start.RAD51D <- 33419240
RAD51D$End.RAD51D <- 33446879

# 56. Use a function to determine if the Start position of RAD51D is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.RAD51D, End.RAD51D) {
  Start.RAD51D >= Start.cnv & Start.RAD51D <= End.cnv
}

df <- RAD51D %>%
  mutate(Test_RAD51D_start = is_between(Start.cnv, End.cnv, Start.RAD51D, End.RAD51D))

# 57. Use a function to determine if the End position of RAD51D is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.RAD51D, End.RAD51D) {
  End.RAD51D >= Start.cnv & End.RAD51D <= End.cnv
}

df2 <- df %>%
  mutate(Test_RAD51D_end = is_between(Start.cnv, End.cnv, Start.RAD51D, End.RAD51D))

# 58. Convert to as.character the resulting columns.
df3 <- df2
df3$Test_RAD51D_start <- as.character(df3$Test_RAD51D_start)
df3$Test_RAD51D_end <- as.character(df3$Test_RAD51D_end)

# 59. Determine if the RAD51D start and end coordinates together fall within a segment/cnv.
# ----- a) Write "1" if the RAD51D start coordinate is within a segment, and "0" if not.
df4 <- df3 %>%
  mutate(RAD51D_start_boolean1 = case_when(
    startsWith (Test_RAD51D_start, "TRUE") ~ 1,
    startsWith (Test_RAD51D_start, "FALSE") ~ 0
  ))

# ----- b) Write "1" if the RAD51D end coordinate is within a segment, and "0" if not.
df5 <- df4 %>%
  mutate(RAD51D_start_boolean2 = case_when(
    startsWith (Test_RAD51D_end, "TRUE") ~ 1,
    startsWith (Test_RAD51D_end, "FALSE") ~ 0
  ))

# ----- c) Sum the numbers of RAD51D_start_boolean1 and RAD51D_start_boolean2.
df6 <- df5
df6$sum_RAD51D <- df6$RAD51D_start_boolean1 + df6$RAD51D_start_boolean2

# ----- d) In another column:
# Write "RAD51D_fully_present" if both RAD51D coordinates fall within a segment.
# Write "RAD51D_partially_present" if only the start or end RAD51D coordinates fall within a segment.
# Write "None" if both RAD51D coordinates DO NOT fall within a segment.
df7 <- df6 %>%
  mutate(RAD51D_presence = case_when(
    sum_RAD51D == 0 ~ "None",
    sum_RAD51D == 1 ~ "RAD51D_partially_present",
    sum_RAD51D == 2 ~ "RAD51D_fully_present",
  ))

# 60. Based on the integer copy number values, determine LOH for each segment/cnv.
df8 <- df7 %>%
  mutate(LOH_status_RAD51D = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-CNV",
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-CNV",
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-CNV",
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous",
  ))

# 61. Get only the desired columns.
df9 <- select(df8, 1:13, 21:22)

# 62. Select only the row/rows where RAD51D is present (if present in a segment).
df10 <- df9[RAD51D_presence %in% c("RAD51D_partially_present",
                                   "RAD51D_fully_present")]

# 63. Add as an extra column to the main dataframe the RAD51D LOH status.
final.df.LOH.RAD51D <- final.df.LOH.accumulated
final.df.LOH.RAD51D$RAD51D_LOH_status <- df10$LOH_status_RAD51D
final.df.LOH.RAD51D$RAD51D_in_segment_status <- df10$RAD51D_presence
final.df.LOH.accumulated <- final.df.LOH.RAD51D


# ---------------------------------------------------------------------------------------------.
# ----------                          PALB2  -  Chr16                                ----------.
# ---------------------------------------------------------------------------------------------.

# 64. Select only the segments/CNVs in chromosome 16 for PALB2.
HRDness <- cnv.abs %>% setDT()
PALB2 <- HRDness[Chr == "chr16"]

# 65. Add as columns the start and end genomic coordinates of PALB2.
PALB2$Start.PALB2 <- 23614486
PALB2$End.PALB2 <- 23652631

# 66. Use a function to determine if the Start position of PALB2 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.PALB2, End.PALB2) {
  Start.PALB2 >= Start.cnv & Start.PALB2 <= End.cnv
}

df <- PALB2 %>%
  mutate(Test_PALB2_start = is_between(Start.cnv, End.cnv, Start.PALB2, End.PALB2))

# 67. Use a function to determine if the End position of PALB2 is within a segment/cnv.
# Print the result as another column.
is_between <- function(Start.cnv, End.cnv, Start.PALB2, End.PALB2) {
  End.PALB2 >= Start.cnv & End.PALB2 <= End.cnv
}

df2 <- df %>%
  mutate(Test_PALB2_end = is_between(Start.cnv, End.cnv, Start.PALB2, End.PALB2))

# 68. Convert to as.character the resulting columns.
df3 <- df2
df3$Test_PALB2_start <- as.character(df3$Test_PALB2_start)
df3$Test_PALB2_end <- as.character(df3$Test_PALB2_end)

# 69. Determine if the PALB2 start and end coordinates together fall within a segment/cnv.
# ----- a) Write "1" if the PALB2 start coordinate is within a segment, and "0" if not.
df4 <- df3 %>%
  mutate(PALB2_start_boolean1 = case_when(
    startsWith (Test_PALB2_start, "TRUE") ~ 1,
    startsWith (Test_PALB2_start, "FALSE") ~ 0
  ))

# ----- b) Write "1" if the PALB2 end coordinate is within a segment, and "0" if not.
df5 <- df4 %>%
  mutate(PALB2_start_boolean2 = case_when(
    startsWith (Test_PALB2_end, "TRUE") ~ 1,
    startsWith (Test_PALB2_end, "FALSE") ~ 0
  ))

# ----- c) Sum the numbers of PALB2_start_boolean1 and PALB2_start_boolean2.
df6 <- df5
df6$sum_PALB2 <- df6$PALB2_start_boolean1 + df6$PALB2_start_boolean2

# ----- d) In another column:
# Write "PALB2_fully_present" if both PALB2 coordinates fall within a segment.
# Write "PALB2_partially_present" if only the start or end PALB2 coordinates fall within a segment.
# Write "None" if both PALB2 coordinates DO NOT fall within a segment.
df7 <- df6 %>%
  mutate(PALB2_presence = case_when(
    sum_PALB2 == 0 ~ "None",
    sum_PALB2 == 1 ~ "PALB2_partially_present",
    sum_PALB2 == 2 ~ "PALB2_fully_present",
  ))

# 70. Based on the integer copy number values, determine LOH for each segment/cnv.
df8 <- df7 %>%
  mutate(LOH_status_PALB2 = case_when(
    CNt == 1 & B_Modal_HSCN_1 == 0 ~ "HEM-CNV",
    CNt == 2 & B_Modal_HSCN_1 == 0 ~ "CN-CNV",
    CNt > 2 & B_Modal_HSCN_1 == 0 ~ "DUP-CNV",
    CNt >= 1 & B_Modal_HSCN_1 > 0 ~ "None",
    CNt == 0 & B_Modal_HSCN_1 == 0 ~ "Homozygous",
  ))

# 71. Get only the desired columns.
df9 <- select(df8, 1:13, 21:22)

# 72. Select only the row/rows where PALB2 is present (if present in a segment).
df10 <- df9[PALB2_presence %in% c("PALB2_partially_present",
                                  "PALB2_fully_present")]

# 73. Add as an extra column to the main dataframe the PALB2 LOH status.
final.df.LOH.PALB2 <- final.df.LOH.accumulated
final.df.LOH.PALB2$PALB2_LOH_status <- df10$LOH_status_PALB2
final.df.LOH.PALB2$PALB2_in_segment_status <- df10$PALB2_presence
final.df.LOH.accumulated <- final.df.LOH.PALB2

# 74. Organize and retain only exonic and splicing variants based on the Func.refGene column.
setDT(final.df.LOH.accumulated)
final.df2 <- final.df.LOH.accumulated
setnames(final.df2, "LOH_status", "Variant_LOH_status")
colnames(final.df2)
final.df.col <- final.df2[, .(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, CNV_Sample_Barcode,
                              Germline_data,
                              Chr, Start, Ref, Alt,
                              Func.refGene, Gene.refGene, ExonicFunc.refGene,
                              AAChange.refGene, ExAC_Freq, gnomAD_genome_ALL,
                              ClinVar_SIG, CADD_phred,
                              dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE, GT,
                              Somatic_data,
                              Start.cnv, End.cnv, Num_Probes, Length,
                              CNt, A_Modal_HSCN_2, B_Modal_HSCN_1, LOH, Variant_LOH_status,
                              Homozygous_deletion, 
                              BRCA1_LOH_status, BRCA1_in_segment_status,
                              BRCA2_LOH_status, BRCA2_in_segment_status,
                              BARD1_LOH_status, BARD1_in_segment_status,
                              RAD51C_LOH_status, RAD51C_in_segment_status,
                              RAD51D_LOH_status, RAD51D_in_segment_status,
                              PALB2_LOH_status, PALB2_in_segment_status,
                              Signature_3_score)]

final.df.col <- final.df.col[Func.refGene %in% c("exonic", "splicing")]


# ----------------------------------------------------------------------------------------------.
################################################################################################.
##### ------------------------------------------------------------------------------------ #####.
#####         Determine second hit in the form of somatic point mutations                  #####.
##### ------------------------------------------------------------------------------------ #####.
################################################################################################.

# 75. Load the full cohort mc3 with all the mutations per sample.
setwd("C:/Users/path/to/dir")
full.cohort <- read_excel(path = "full.somatic.cohort.mc3.xlsx", sheet = "Hoja1")
setDT(full.cohort)

# 76. Select only the sample in question of the full cohort mc3.
sample.mc3 <- full.cohort[Matched_Norm_Sample_Barcode == "TCGA-XX"]

# 77. Select only the somatic mutations that have either:
# - P/LP classification in ClinVar.
# - A truncating effect.
# - Any type of mutation in the 6 HRD genes: BARD1, BRCA1/2, RAD51C/D, PALB2.
sample.mc3.2 <- sample.mc3

sample.mc3.2 <- sample.mc3.2 %>%
  mutate(New_SPM = case_when(
    startsWith(Variant_Classification, "Frame_Shift_Del") ~ 1,
    startsWith(Variant_Classification, "Frame_Shift_Ins") ~ 1,
    startsWith(Variant_Classification, "Nonsense_Mutation") ~ 1,
    startsWith(Variant_Classification, "Splice_Site") ~ 1,
    startsWith(CLIN_SIG, "pathogenic") ~ 1,
    startsWith(CLIN_SIG, "not_provided,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "likely_pathogenic") ~ 1,
    startsWith(CLIN_SIG, "likely_pathogenic,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "uncertain_significance,benign,likely_benign,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "uncertain_significance,not_provided,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "likely_pathogenic,pathogenic,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "uncertain_significance,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "uncertain_significance,likely_pathogenic") ~ 1,
    startsWith(CLIN_SIG, "pathogenic,other") ~ 1,
    startsWith(CLIN_SIG, "likely_benign,pathogenic") ~ 1,
    startsWith(CLIN_SIG, "pathogenic,uncertain_significance") ~ 1,
    startsWith(CLIN_SIG, "uncertain_significance,likely_pathogenic,pathogenic") ~ 1,
    startsWith(Hugo_Symbol, "BARD1") ~ 1,
    startsWith(Hugo_Symbol, "BRCA1") ~ 1,
    startsWith(Hugo_Symbol, "BRCA2") ~ 1,
    startsWith(Hugo_Symbol, "RAD51C") ~ 1,
    startsWith(Hugo_Symbol, "RAD51D") ~ 1,
    startsWith(Hugo_Symbol, "PALB2") ~ 1
  ))

sample.mc3.2$New_SPM[is.na(sample.mc3.2$New_SPM)] <- "No"
sample.mc3.3 <- sample.mc3.2[New_SPM == "1"]
sample.mc3.4 <- sample.mc3.3[, !c("New_SPM")]
sample.mc3 <- sample.mc3.4
rm(sample.mc3.2, sample.mc3.3, sample.mc3.4)

# 78. From such germline sample, select all the germline genes that match in its respective tumor sample.
# The resulting rows / somatic mutations will be the genes mutated in the germline that are also mutated in the tumor.
sample.germ <- final.df.col[Gene.refGene %in% sample.mc3$Hugo_Symbol]

# 79. Now get the somatic information of such matching germline gene variants (second hits - point mutations).
sample.germ.tum <- sample.mc3[Hugo_Symbol %in% sample.germ$Gene.refGene]
sample.germ.tum.concat <- sample.germ.tum
sample.germ.tum.concat$Second_hit_point_mutations <- paste(sample.germ.tum.concat$Hugo_Symbol,
                                                           "_",
                                                           sample.germ.tum.concat$Variant_Classification,
                                                           sep = "")

sample.concat.col <- sample.germ.tum.concat[,.(Second_hit_point_mutations)]
sample.concat.col.vec <- sample.concat.col[,Second_hit_point_mutations]
final.df.col$Second_hit_point_mutations <- list(sample.concat.col.vec)

colnames(final.df.col)
final.df.col <- final.df.col[, .(Tumor_Sample_Barcode,	Matched_Norm_Sample_Barcode, CNV_Sample_Barcode,
                                 Germline_data,	Chr,	Start, Ref,	Alt, Func.refGene,
                                 Gene.refGene, ExonicFunc.refGene, AAChange.refGene,
                                 ExAC_Freq, gnomAD_genome_ALL, ClinVar_SIG,
                                 CADD_phred, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE,
                                 GT, Somatic_data, Start.cnv, End.cnv, Num_Probes, Length,
                                 CNt, A_Modal_HSCN_2, B_Modal_HSCN_1, LOH, Variant_LOH_status, Homozygous_deletion,
                                 Second_hit_point_mutations,
                                 BRCA1_LOH_status, BRCA1_in_segment_status, BRCA2_LOH_status, BRCA2_in_segment_status,
                                 BARD1_LOH_status, BARD1_in_segment_status, RAD51C_LOH_status, RAD51C_in_segment_status,
                                 RAD51D_LOH_status, RAD51D_in_segment_status, PALB2_LOH_status, PALB2_in_segment_status,
                                 Signature_3_score)]


# ----------------------------------------------------------------------------------------------.
################################################################################################.
##### ------------------------------------------------------------------------------------ #####.
#####              Integrate methylation data of BRCA1 and RAD51C                          #####.
##### ------------------------------------------------------------------------------------ #####.
################################################################################################.

# 80. Read the methylation file.
dir_methylation <- "C:/Users/path/to/dir/methylation_data.tsv"
methylation <- fread(dir_methylation) %>% 
  t() %>% 
  as.data.table(.,keep.rownames = TRUE)

names(methylation) <- methylation[1,] %>% unlist
methylation_2 <- copy(methylation[-1,])
setnames(methylation_2, "Symbols", "Sample")

# 81. Adjust the IDs substring (names of each IDs) to enable merging methylation data with exome data.
methylation_2[, Methylated_Sample_Barcode := str_extract(Sample,"TCGA-..-....-...-")]

# 82. Select only methylation in BRCA1 and RAD51C.
methylation_3 <- methylation_2[, .(Sample, Methylated_Sample_Barcode, BRCA1, RAD51C)]
unique(methylation_3$BRCA1)
methylation_3$BRCA1 <- gsub(" 0", "0", methylation_3$BRCA1)
methylation_3$BRCA1 <- gsub(" 1", "1", methylation_3$BRCA1)
methylation_3$BRCA1[is.na(methylation_3$BRCA1)] <- "Unavailable"
methylation_3$RAD51C[is.na(methylation_3$RAD51C)] <- "Unavailable"
unique(methylation_3$BRCA1)
unique(methylation_3$RAD51C)
setnames(methylation_3, "Sample", "Tumor_Sample_Barcode")
setnames(methylation_3, "BRCA1", "BRCA1_methylation")
setnames(methylation_3, "RAD51C", "RAD51C_methylation")

# 83. Adjust the IDs' headers of the sample in question.
final.df.col[, Methylated_Sample_Barcode := str_extract(
  Tumor_Sample_Barcode,"TCGA-..-....-...-")]

# 84. Get only the methylation of the sample in question.
methylation_4 <- methylation_3[Methylated_Sample_Barcode %in%
                                 final.df.col$Methylated_Sample_Barcode]

# 85. Merge the methylation data with the TCGA-BRCA data.
final.df.col.methyl <- merge(final.df.col, methylation_4,
                             by = "Methylated_Sample_Barcode",
                             all.x = TRUE, all.y = TRUE)

# 86. Rename the Tumor ID header correctly and reorganize the order of the columns.
setnames(final.df.col.methyl, "Tumor_Sample_Barcode.x", "Tumor_Sample_Barcode")
colnames(final.df.col.methyl)
final.df.col.methyl <- final.df.col.methyl[, .(Tumor_Sample_Barcode,	Matched_Norm_Sample_Barcode, Methylated_Sample_Barcode, CNV_Sample_Barcode,
                                               Germline_data,	Chr,	Start, Ref,	Alt, Func.refGene,
                                               Gene.refGene, ExonicFunc.refGene, AAChange.refGene,
                                               ExAC_Freq, gnomAD_genome_ALL, ClinVar_SIG,
                                               CADD_phred, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE,
                                               GT, Somatic_data, Start.cnv, End.cnv, Num_Probes, Length,
                                               CNt, A_Modal_HSCN_2, B_Modal_HSCN_1, LOH, Variant_LOH_status, Homozygous_deletion,
                                               Second_hit_point_mutations,
                                               BRCA1_LOH_status, BRCA1_in_segment_status, BRCA2_LOH_status, BRCA2_in_segment_status,
                                               BARD1_LOH_status, BARD1_in_segment_status, RAD51C_LOH_status, RAD51C_in_segment_status,
                                               RAD51D_LOH_status, RAD51D_in_segment_status, PALB2_LOH_status, PALB2_in_segment_status,
                                               BRCA1_methylation, RAD51C_methylation, Signature_3_score)]


# ----------------------------------------------------------------------------------------------.
################################################################################################.
##### ------------------------------------------------------------------------------------ #####.
#####                  Integrate the corresponding clinical data                           #####.
##### ------------------------------------------------------------------------------------ #####.
################################################################################################.

# 87. Select the columns having the clinical data of interest of the sample in question.
sample.mc3.2 <- sample.mc3
sample.mc3.2 <- select(sample.mc3.2, 11:13, 3:10)
sample.mc3.3 <- distinct(sample.mc3.2, Tumor_Sample_Barcode, .keep_all = TRUE)

# 88. Merge the clinical data with the main dataframe of the sample in question.
final.df.col.merging <- merge(final.df.col.methyl, sample.mc3.3,
                              by = "Tumor_Sample_Barcode",
                              all.x = TRUE,
                              all.y = TRUE)

final.df.col.merging2 <- final.df.col.merging
final.df.col.merging2$Clinical_data <- "Clinical"
setnames(final.df.col.merging2, "Total_number_of_mutations", "Total_number_of_somatic_mutations")
setnames(final.df.col.merging2, "Matched_Norm_Sample_Barcode.x", "Matched_Norm_Sample_Barcode")
setnames(final.df.col.merging2, "Cancer Type Detailed", "Cancer_Type_Detailed")
setnames(final.df.col.merging2, "Overall Survival Status", "Overall_Survival_Status")
setnames(final.df.col.merging2, "Ethnicity Category", "Ethnicity_Category")
setnames(final.df.col.merging2, "Tumor Type", "Tumor_Type")
setnames(final.df.col.merging2, "Diagnosis Age", "Diagnosis_Age")
setnames(final.df.col.merging2, "Race Category", "Race_Category")

# 89. Do a final reorder of columns.
colnames(final.df.col.merging2)
final.df.col.merging2 <- final.df.col.merging2[, .(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Methylated_Sample_Barcode, CNV_Sample_Barcode,
                                                   Clinical_data,
                                                   Cancer_Type_Detailed, Diagnosis_Age, Sex, Overall_Survival_Status,
                                                   Ethnicity_Category, Race_Category, Subtype, Tumor_Type,
                                                   Germline_data,
                                                   Chr, Start, Ref, Alt,
                                                   Func.refGene, Gene.refGene, ExonicFunc.refGene,
                                                   AAChange.refGene, ExAC_Freq, gnomAD_genome_ALL,
                                                   ClinVar_SIG, CADD_phred, dbscSNV_ADA_SCORE,
                                                   dbscSNV_RF_SCORE, GT,
                                                   Somatic_data,
                                                   Start.cnv, End.cnv, Num_Probes,
                                                   Length, CNt, A_Modal_HSCN_2, B_Modal_HSCN_1, LOH, Variant_LOH_status,
                                                   Homozygous_deletion, Second_hit_point_mutations, BRCA1_LOH_status,
                                                   BRCA1_in_segment_status, BRCA2_LOH_status, BRCA2_in_segment_status,
                                                   BARD1_LOH_status, BARD1_in_segment_status, RAD51C_LOH_status,
                                                   RAD51C_in_segment_status, RAD51D_LOH_status, RAD51D_in_segment_status,
                                                   PALB2_LOH_status, PALB2_in_segment_status,
                                                   BRCA1_methylation, RAD51C_methylation, Total_number_of_somatic_mutations, Signature_3_score)]

class(final.df.col.merging2$Second_hit_point_mutations)
final.df.col.merging2$Second_hit_point_mutations <- as.character(final.df.col.merging2$Second_hit_point_mutations)

# 90. Write the final output as xlsx and txt files.
dir.create("C:/Users/path/to/file/output")
setwd("C:/Users/path/to/file/output")
openxlsx::write.xlsx(final.df.col.merging2,
                     file = "TCGA-XX_integrated_multidimensional_data.xlsx",
                     colNames = T, rowNames = F, append = F)

setwd("C:/Users/path/to/file/output")
write.table(final.df.col.merging2,
            file = "TCGA-XX_integrated_multidimensional_data.txt",
            sep = "\t", col.name = T, row.names = F, quote = F)


###############################################################################################.
###############################################################################################.
##########                                                                           ##########.
##########                           T H E     E N D                                 ##########.
##########                                                                           ##########.
###############################################################################################.
###############################################################################################.

