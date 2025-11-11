library(tidyverse)
library(stringr)
library(ggfortify)
library(cluster)
library(factoextra)
library(magrittr)
library(tibble)
library(NbClust)
library(viridis)
library(dendextend)
library(heatmaply)
library(superheat)


# ============ Dataset import and formatting ======================== # 
filepath <- "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Be2_analysis\\EpiPro_2.2.csv"
histone_data <- read.csv(filepath, skip=1) # <== read in the dataframe 
labels <- read.csv(filepath, header = FALSE)[1,2:ncol(histone_data)]   # <== extract the label cols

# Sample ID import
filepath <- "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Be2_analysis\\SampleIDs_modified_cleannames.csv"
sample_ids <- read.csv(filepath)

# selecting only area columns (we use this so we can normalize by total intensity rather than peptide family)
histone_data <- histone_data %>% 
  select(1, contains("Area")) %>% 
  filter(!is.na(Area))

# Cleaning column names
mappings <- setNames(sample_ids$clean_name, sample_ids$Sample)
assign_names <- function(x){
  # Stupid R not having dictionaries :(
  if (x %in% names(mappings)){return(mappings[[x]])}else{return(x)}}

colnames(histone_data)[-1] <- labels
colnames(histone_data) <- lapply(colnames(histone_data), function(x){sub("^[^,]*,", "",x)})
colnames(histone_data) <- lapply(colnames(histone_data),
                                 assign_names)

# Clean Peptide names 
histone_data$Peptide <- unlist(lapply(histone_data$Peptide, function(x){unlist(str_replace_all(x, "[ :().]", "_"))}))
histone_data <- as_tibble(histone_data)


rm(labels)
rm(assign_names)


histone_data <- histone_data %>% 
  gather(Sample, value, 2:ncol(histone_data)) %>% 
  spread(Peptide, value)


# ========== Removal of bad samples
bad_samples <- c(
  #"DPI_24h-B1-T1" # Seemed to be very very different to the others 
                   # BioRep 2 seemed to be very different from the others 
                   # - may add these back in later though....
  #,"NC-B2-T1"
  #,"NC-B2-T2"
)


# filter out samples we identified before. 
histone_data <- histone_data %>% 
  filter(!histone_data$Sample %in% bad_samples)


# ========== Deconvoluting modifications 
# I chose not to do this in order to try and recreate Stephanies analysis. 
deconvolution = TRUE

if (deconvolution){
  histone_data <- histone_data %>% 
    mutate(H3_9_17_K9K14ac = apply(histone_data[c("H3_9_17_K9ac", "H3_9_17_K14ac")], 1, sum, na.rm = TRUE),
           H3_18_26_K18K23ac = apply(histone_data[c("H3_18_26_K18ac", "H3_18_26_K23ac")], 1, sum, na.rm = TRUE),
           H2A1_4_11_K5K9ac = apply(histone_data[c("H2A1_4_11_K5ac", "H2A1_4_11_K9ac")], 1, sum, na.rm = TRUE),
           H2AX_4_11_K5K9ac = apply(histone_data[c("H2AX_4_11_K5ac", "H2AX_4_11_K9ac")], 1, sum, na.rm = TRUE),
           H2AJ_4_11_K5K9ac = apply(histone_data[c("H2AJ_4_11_K5ac", "H2AJ_4_11_K9ac")], 1, sum, na.rm = TRUE)) %>% 
    select(-matches("H3_18_26_K18ac$")) %>% 
    select(-matches("H3_18_26_K23ac$")) %>%
    select(-matches("H3_9_17_K9ac$")) %>%
    select(-matches("H3_9_17_K14ac$")) %>% 
    select(-matches("H2A1_4_11_K5ac$")) %>% 
    select(-matches("H2A1_4_11_K9ac$")) %>% 
    select(-matches("H2AX_4_11_K5ac$")) %>% 
    select(-matches("H2AX_4_11_K9ac$")) %>% 
    select(-matches("H2AJ_4_11_K5ac$")) %>% 
    select(-matches("H2AJ_4_11_K9ac$"))
  
  # Tranpose back
  histone_data <- histone_data %>% 
    gather(Peptide, value, 2:ncol(histone_data)) %>% 
    spread(Sample, value)
  
  # Select all peptides around H4[4-17] that are acetylated
  H4 <- histone_data %>% 
    filter(str_detect(Peptide, "H4.4.17.*ac"))
  
  
  # Create a new col in this df that counts the number of acs found... 
  H4$Acetyl_count <- str_count(H4$Peptide, "ac")
  
  # group Peptides by the number of acetyls found on each and sum them 
  H4 <- H4 %>%
    select(-Peptide) %>% 
    group_by(Acetyl_count) %>% 
    summarise_all(sum) %>%
    arrange(Acetyl_count)
  
  # recreate the Peptide column, but now with 
  H4$Peptide = (c("H4_4_17_1ac", "H4_4_17_2ac", "H4_4_17_3ac", "H4_4_17_4ac"))
  
  
  H4 <- H4 %>% 
    ungroup() %>%  # ungroups the data 
    select(-Acetyl_count) %>%  # removes the acetyl_count column. 
    select(Peptide, everything()) # puts peptide on the right and everything else after. 
  
  
  histone_data <- rbind(H4, histone_data) %>%  # similar to pd.concat - this slides our new data back into the df 
    filter(!str_detect(Peptide, "H4.4.17.K.*ac")) # this now drops the old data 

  
} else {
  # Tranpose back regardless
  histone_data <- histone_data %>% 
    gather(Peptide, value, 2:ncol(histone_data)) %>% 
    spread(Sample, value)
  
}



# ========== Normalisation 
norm_method <- "%PeptideFamily"
#norm_method <- "total_intensity"



# ========== Removal of missing Peptides 
# Instead of looking for missing mods, I instead want to only keep mods that have EITHER:
# 1. present in at least 30% of any group (present in 8 of 23 data columns)
# this is the value I set for maximum allowable zeros
max_zeros <- 8
# 2. 4 or more values in at least one group
min_in_onegroup <- 4

# 1. Find the peptides that have at least 8 valid values
# create counts and mean columns
histone_data <- histone_data %>% 
  mutate(
    Mean = rowMeans(histone_data[,-1]), 
    zero_counts = rowSums(histone_data == 0)
  )


# Get values that were excluded for being missing...
remove_toomanyzeros <- filter(histone_data, zero_counts > max_zeros)$Peptide


# 2. Find the peptides that had at least 4 valid values for one group
# reshape
histone_data <- histone_data %>%
  select(-zero_counts, -Mean)
histone_data <- histone_data %>% 
  gather(Sample, value, 2:ncol(histone_data)) %>% 
  spread(Peptide, value)

# Group by biological condition
groups <- merge(histone_data, 
                select(sample_ids,c(clean_name,Condition)),
                by.x = "Sample",
                by.y ="clean_name")
group_flags <- groups %>%
  mutate(across(
            matches("^H\\d"), # selects the columns to apply across
            ~ ifelse(.x==0, 
                    0,
                    1) ))     # puts a 1 wherever there was a valid value 
validpergroup <- group_flags %>%
  group_by(Condition) %>%
  summarise(across(matches("^H\\d"),
                   sum
                   ))

# filter for columns which have at least 1 value >=4
one_valid_group <- names(validpergroup)[sapply(validpergroup,function(x) any(x >= min_in_onegroup, na.rm = TRUE))]

# rescue these peptides from the remove list
rescued <- intersect(one_valid_group,remove_toomanyzeros)
remove_toomanyzeros <- remove_toomanyzeros[!remove_toomanyzeros %in% rescued]


# FINALLY: remove all the peptides from this list 
histone_data <- histone_data %>%
  select(-all_of(remove_toomanyzeros))


rm(validpergroup)






# ================= Now finally do the normalisation after missing values/low abundance values have been weeded out
histone_data <- histone_data %>% 
  gather(Peptide, value, 2:ncol(histone_data)) %>% 
  spread(Sample, value)

# Caculating sums for each peptide family
# To do this, first make a label for each peptide family
histone_data$peptide_id <-unlist(lapply(histone_data$Peptide, function(x){
  unlist(str_sub(x, 1, 7)) # <--- groups by the first 7 letters in their name
}))

write_csv(histone_data,"C:\\Users\\Tim\\OneDrive\\__ThesisProject\\Results\\Appendix1_PoorSamples\\Be2_filtered_ids.csv")


if (norm_method == "%PeptideFamily"){
  print("normalisation by Peptide Family")
  normalized <- histone_data %>% 
    group_by(peptide_id) %>% 
    mutate_if(is.double, list(~./sum(., na.rm = TRUE))) %>% 
    ungroup() %>% 
    select(-peptide_id)
}

if (norm_method == "total_intensity"){
  print("normalisation by total intensity")
  normalized <- histone_data %>%
    ungroup() %>%
    mutate_if(is.double, list(~./sum(., na.rm = TRUE))) %>%
    select(-peptide_id)
}



#write_csv(normalized,"C:\\Users\\Tim\\OneDrive\\__ThesisProject\\Results\\Appendix1_PoorSamples\\Be2_normdataforchecking.csv")








# ========== Removal of Problematic (too variable) peptides
# Grouping by sample condition
normalized <- normalized %>% 
  gather(Sample, value, 2:ncol(normalized)) %>% 
  spread(Peptide, value)

sample_ids <- rename(sample_ids,
                     #newname = oldname
                     sample_id=Sample)
sample_ids <- rename(sample_ids, 
                     Sample=clean_name)

# merge with the sample_ids dataframe to map to group names 
groups <- merge(normalized, sample_ids, by ="Sample")

groups <- rename(groups, Sample_Type = Condition )



# Coefficient of variation is defined as std/mean
# We calculate CV for each biological conditions and then look at the median value for each. 
# If a modification tends to have a very high CV we filter it out 
# (we will use 0.8 as a cutoff - meaning the average value is 0.8x the standard dev - quite a lot of variation)
cv_cutoff <- 10000


# calculate avg and st.dev
average <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(-c(Sample))
stdev <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  select(-c(Sample))
# calculate CV for each peptide for each group, then get the median CV for each peptide
CV <- stdev[2:ncol(stdev)]/average[2:ncol(average)]
CV <- cbind(average[,1], CV)
CV <- CV %>%
  summarise_all(median, na.rm = TRUE)
CV <- CV %>%
  gather(Peptide, value, 2:ncol(CV)) %>% 
  spread(Sample_Type, value)
colnames(CV) <- c("Peptide", "Median")

# filter out the values above the cutoff from the original dataframe
CV <- CV %>% 
  filter(Median > cv_cutoff)
histone_data <- histone_data %>% 
  filter(!Peptide %in% CV$Peptide)


# ========== Impute the missing values 
impute_val <- 1E5
histone_data[is.na(histone_data)] <- 0
histone_data[histone_data == 0] <- impute_val



# Remove peptides with only one modification identified ... nothing to compare
# creates a new column of counts per id
histone_data <- histone_data %>%
  group_by(peptide_id) %>% 
  mutate(ids = n())

print("peptides with only one form identified (nothing to compare)")
# prints out the table of peptides that are now excluded.
print(bind_cols(histone_data[histone_data$ids <2,]$peptide_id,
                histone_data[histone_data$ids <2,]$ids)
)
histone_data <- histone_data %>% 
  filter(ids >= 2) %>% 
  select(-ids)






# ==================== renormalise data after removing problematic peptides and singly identified peptides 
new <- histone_data


if (norm_method == "%PeptideFamily"){
  print("normalisation by Peptide Family")
  normalized <- histone_data %>% 
    group_by(peptide_id) %>% 
    mutate_if(is.double, list(~./sum(., na.rm = TRUE))) %>% 
    ungroup() %>% 
    select(-peptide_id)
}

if (norm_method == "total_intensity"){
  print("normalisation by total intensity")
  normalized <- histone_data %>%
    ungroup() %>%
    mutate_if(is.double, list(~./sum(., na.rm = TRUE))) %>%
    select(-peptide_id)
}



# ==================== Validate consistency within Samples

# Transpose the data
normalized <- normalized %>% 
  gather(Sample, value, 2:ncol(normalized)) %>% 
  spread(Peptide, value)

# convert to dataframe type. 
distance <- as.data.frame(normalized)
# sets the samples as row names (indexes if we were talking pandas) rather than a column
distance <- column_to_rownames(distance, var = "Sample")

# get the pairwise spearman rank correlation coefficients between samples... 
res.dist <- get_dist(distance, method = "spearman")

fviz_dist(res.dist, lab_size = 5, order = FALSE)


library(pheatmap)
res.mat <- as.matrix(res.dist)

pdf("spearman_heatmap.pdf", width = 10, height = 8)  # set size in inches
# heatmap
pheatmap(res.mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Spearman distance Between Samples")

rm(res.mat)

# ============================= End data cleaning.... 
# Re-group and calculate averages again
groups <- merge(normalized, select(sample_ids, c(Condition, Sample)), by ="Sample")
groups <- rename(groups,
                 Sample_Type = Condition)


write.csv(groups,
  "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Histone_be2IMR_comparative\\Be2groups.csv"
)

average <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(-c(Sample))

write.csv(average,
          "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Histone_be2IMR_comparative\\Be2groups_average.csv"
)

stdev <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  select(-c(Sample))
stdev$Sample_Type <- unlist(lapply(stdev$Sample_Type,function(x){unlist(str_replace_all(x, "$", "_Sd"))}))




# Calculate Fold-Changes
average <- average %>% 
  gather(Peptide, value, 2:ncol(average)) %>% 
  spread(Sample_Type, value)

fold_change <- average

fold_change <- fold_change %>% 
  select_if(is.numeric) %>% 
  map(function(x){return(x/(average %>% select_if(is.numeric # <= function to map with
  )))}) %>%
  imap(~set_names(.x, paste0(.y, "_", names(.x)))) %>%  # <-- the part that 
  bind_cols() # <--- turns it back into a tibble  

fold_change <- bind_cols(average[1], fold_change) # <- add the peptide column back 




# ============================= Student t-tests with FDR for control vs DPI 
ctrl_dpi <- groups[groups$Sample_Type %in% c("NC","DPI 24h"),]

long <- ctrl_dpi %>% 
  select(-c(Sample)) %>% 
  select(Sample_Type, everything())

long <- long %>%   
  gather(Peptide, value, 2:ncol(long)) 
colnames(long) <- c("Sample_Type", "Peptide", "value")


# We convert our data to factors 
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))


# Temporary function for running t-tests
tmpfn_tim <- function(dat) { # <== takes as input the dataframe for that modification
  # runs t-test
  t_result <- t.test(
    value ~ Sample_Type, 
    data = dat,
    var.equal = TRUE     # <=== this setting is commented if you want Welch's t-test
  )
  p_val <- t_result$p.value
  # returns as dataframe, so that bind_rows can then stack them all at the end... 
  return(data.frame(p_val))
}

# run t-tests
long_stat_tim <- long_3 %>%
  map(tmpfn_tim) %>%
  bind_rows(.id = "Peptide")

long_stat_tim$fdr_padj <- p.adjust(long_stat_tim$p_val, method = "fdr")

# long_stat_tim$fdr_padj <- long_stat_tim$p_val

write_csv(long_stat_tim, "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Be2_analysis\\Be2_t-testvscontrol_students_newscript.csv")





# ============================= ANOVA with Tukey post-hoc

long <- groups %>% 
  select(-c(Sample)) %>% 
  select(Sample_Type, everything())
long <- long %>%   
  gather(Peptide, value, 2:ncol(long)) 
colnames(long) <- c("Sample_Type", "Peptide", "value")


# We convert our data to factors 
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))

# Temp function for running anova
tmpfn <- function(dat) {
  fit <- lm(value ~ Sample_Type, dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}


long_stat <- long_3 %>%
  map(tmpfn) %>%
  bind_rows(.id = "Peptide")

long_stat <- long_stat %>%
  filter(pval < 0.05)
list <- c(long_stat[,1])

long <- filter(long, Peptide %in% list)


# Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

tmpfn <- function(dat) {
  fit <- aov(value ~ Sample_Type, dat)
  Tukey <- TukeyHSD(fit)$Sample_Type
}

long_stat <- long_3 %>%
  map(tmpfn)

p_adj <- data.frame(long_stat)
p_adj <- p_adj %>%
  select(contains("p.adj"))
names(p_adj) <- sub(".p.adj", "", names(p_adj))
p_adj <- rownames_to_column(p_adj, "Sample")

rm(list)
rm(tmpfn)
rm(long_2)
rm(long_3)
rm(long_stat)




# 1. pivot the stdev table
stdev <- stdev %>% 
  gather(Peptide, value, 2:ncol(stdev)) %>% 
  spread(Sample_Type, value)


# 2. pivot the fold_change table
fold_change <- fold_change %>% 
  gather(Sample, value, 2:ncol(fold_change)) %>% 
  spread(Peptide, value)



# Filter out self-comparisons from fold_change
list <- c(p_adj[,1])



fold_change <- dplyr::filter(
  fold_change,
  gsub("_", "-", Sample) %in% gsub("_", "-", list)
)

print("sample_names")
print(fold_change$Sample) 
# Reformat the sample names
fold_change$Sample <- unlist(lapply(fold_change$Sample, function(x){unlist(str_replace_all(x, "_", "/"))}))

# pivot back the df
fold_change <- fold_change %>%
  gather(Peptide, value, 2:ncol(fold_change)) %>% 
  spread(Sample, value)
rm(list)


# 3. pivot the p_adj dataframe. 
p_adj <- p_adj %>%
  gather(Peptide, value, 2:ncol(p_adj)) %>% 
  spread(Sample, value)
# Clean the column names in p_adj 
p_adj$Peptide <- unlist(lapply(p_adj$Peptide, function(x){unlist(str_replace_all(x, "\\.Sample_Type", ""))}))


# Join the dataframes
final <- average %>% 
  full_join(stdev, by = "Peptide") %>% 
  full_join(fold_change, by = "Peptide") %>%
  full_join(p_adj, by = "Peptide")


write.csv(final, file = "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Be2_analysis\\output_df_2_newscript.csv")




### ======================== Further Deconvolution of peptides
# From what I can see, this next section adapted from the Thomas et al script is just 

deconvoluted <- groups %>%  
  ungroup() %>% 
  mutate(H3_9_18_K9me1 = select(groups, contains("K9me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_9_18_K9me2 = select(groups, contains("K9me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_9_18_K9me3 = select(groups, contains("K9me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K27me1 = select(groups, contains("H3_27_40_K27me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K27me2 = select(groups, contains("H3_27_40_K27me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K27me3 = select(groups, contains("H3_27_40_K27me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K36me1 = select(groups, matches("H3_.*K36me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K36me2 = select(groups, matches('H3_.*K36me2')) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K36me3 = select(groups, matches("H3_.*K36me3")) %>%
           rowSums(na.rm = TRUE)) %>% 
  mutate(H33_27_41_K27me1 = select(groups, contains("H33_27_40_K27me1")) %>% 
           rowSums(na.rm = TRUE)) %>% 
  mutate(H33_27_41_K27me2 = select(groups, contains("H33_27_40_K27me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K27me3 = select(groups, contains("H33_27_40_K27me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K36me1 = select(groups, matches("H33_.*K36me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K36me2 = select(groups, matches("H33_.*K36me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K36me3 = select(groups, matches("H33_.*K36me3")) %>% 
           rowSums(na.rm = TRUE))

deconvoluted <- deconvoluted %>% 
  select(-matches("H3_9_17_K9me(1|2|3)")) %>% 
  select(-matches("H3_27_40_K(27|36)(me1|me2|me3)")) %>% 
  select(-matches("H33_27_40_K(27|36)(me1|me2|me3)"))

write.csv(deconvoluted,
          "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Histone_be2IMR_comparative\\Be2groups_deconvol.csv"
)


# We then calculate the average, stan_dev and fold-changes

average_deconv <- deconvoluted %>% 
  group_by(Sample_Type) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(-c(Sample))

write.csv(average_deconv,
          "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Histone_be2IMR_comparative\\Be2groups_deconvol_avg.csv"
)


stdev_deconv <- deconvoluted %>% 
  group_by(Sample_Type) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  select(-c(Sample))
stdev_deconv$Sample_Type <- unlist(lapply(stdev_deconv$Sample_Type, function(x){unlist(str_replace_all(x, "$", "_Sd"))}))


average_deconv <- average_deconv %>% 
  gather(Peptide, value, 2:ncol(average_deconv)) %>% 
  spread(Sample_Type, value)


# Calculate fold changes
fold_change_deconv <- average_deconv
fold_change_deconv <- fold_change_deconv %>% 
  select_if(is.numeric) %>% 
  map(~./fold_change_deconv %>% select_if(is.numeric)) %>% 
  imap(~set_names(.x, paste0(.y, "-", names(.x)))) %>% 
  bind_cols()
fold_change_deconv <- bind_cols(average_deconv[1], fold_change_deconv)


# ============================= Student t-tests with FDR for control vs DPI 
ctrl_dpi <- groups[groups$Sample_Type %in% c("NC","DPI 24h"),]

long <- ctrl_dpi %>% 
  select(-c(Sample)) %>% 
  select(Sample_Type, everything())

long <- long %>%   
  gather(Peptide, value, 2:ncol(long)) 
colnames(long) <- c("Sample_Type", "Peptide", "value")


# We convert our data to factors 
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))


# Temporary function for running t-tests
tmpfn_tim <- function(dat) { # <== takes as input the dataframe for that modification
  # runs t-test
  t_result <- t.test(
    value ~ Sample_Type, 
    data = dat,
    var.equal = TRUE     # <=== this setting is commented if you want Welch's t-test
  )
  p_val <- t_result$p.value
  # returns as dataframe, so that bind_rows can then stack them all at the end... 
  return(data.frame(p_val))
}

# run t-tests
long_stat_tim_deconv <- long_3 %>%
  map(tmpfn_tim) %>%
  bind_rows(.id = "Peptide")

long_stat_tim_deconv$fdr_padj <- p.adjust(long_stat_tim_deconv$p_val, method = "fdr")


write_csv(long_stat_tim_deconv, "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Be2_analysis\\be2t-testvscontrol_students_deconv.csv")
























# ======================= One-way Anova with Tukey post-hoc correction
long_deconv <- deconvoluted %>% 
  select(-c(Sample)) %>% 
  select(Sample_Type, everything())
long_deconv <- long_deconv %>%   
  gather(Peptide, value, 2:ncol(long_deconv)) 
colnames(long_deconv) <- c("Sample_Type", "Peptide", "value")

long_2 <- long_deconv
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

# Since some of my data is missing, i use my own version of the function..... 
tmpfn <- function(dat) { # <== takes as input the dataframe for that modification
  # Fits ANOVA model 
  peptide_anova <- aov( 
    data = dat, 
    value ~ Sample_Type
  )
  # Summaries ANOVA model 
  stat <- summary(peptide_anova)   
  
  # extracts p value 
  pval <- stat[[1]][["Pr(>F)"]][1]
  
  # returns as dataframe, so that bind_rows can then stack them all at the end... 
  return(data.frame(pval))
}

long_stat <- long_3 %>%
  map(tmpfn) %>%
  bind_rows(.id = "Peptide")


long_stat <- long_stat %>%
  filter(pval < 0.05)
list <- c(long_stat[,1])
long_deconv <- filter(long_deconv, Peptide %in% list)

# Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
long_2 <- long_deconv
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

tmpfn <- function(dat) {
  fit <- aov(value ~ Sample_Type, dat)
  Tukey <- TukeyHSD(fit)$Sample_Type
}

long_stat <- long_3 %>%
  map(tmpfn)

p_adj_deconv <- data.frame(long_stat)
p_adj_deconv <- p_adj_deconv %>% 
  select(contains("p.adj"))
names(p_adj_deconv) <- sub(".p.adj", "", names(p_adj_deconv))
p_adj_deconv <- rownames_to_column(p_adj_deconv, "Sample")


rm(list)
rm(tmpfn)
rm(long_2)
rm(long_3)
rm(long_stat)

# Combine average, standard deviation, fold change, and p-values into one table
# First, make sure all your data is transposed
stdev_deconv <- stdev_deconv %>% 
  gather(Peptide, value, 2:ncol(stdev_deconv)) %>% 
  spread(Sample_Type, value)

# Filter out the extra fold change comparisons we don't need
fold_change_deconv <- fold_change_deconv %>% 
  gather(Sample, value, 2:ncol(fold_change_deconv)) %>% 
  spread(Peptide, value)
list <- c(p_adj_deconv[,1])
fold_change_deconv <- filter(fold_change_deconv, Sample %in% list)
fold_change_deconv$Sample <- unlist(lapply(fold_change_deconv$Sample, function(x){unlist(str_replace_all(x, "-", "/"))}))
fold_change_deconv <- fold_change_deconv %>%
  gather(Peptide, value, 2:ncol(fold_change_deconv)) %>% 
  spread(Sample, value)
rm(list)
p_adj_deconv <- p_adj_deconv %>%
  gather(Peptide, value, 2:ncol(p_adj_deconv)) %>% 
  spread(Sample, value)
p_adj_deconv$Peptide <- unlist(lapply(p_adj_deconv$Peptide, function(x){unlist(str_replace_all(x, "\\.Sample_Type", ""))}))

# Then combine
final_deconv <- average_deconv %>% 
  full_join(stdev_deconv, by = "Peptide") %>% 
  full_join(fold_change_deconv, by = "Peptide") %>%
  full_join(p_adj_deconv, by = "Peptide")

# Write everything to an excel file
write.csv(final_deconv, file = "C:\\Users\\Tim\\OneDrive\\__ThesisProject\\_DryLab\\Code\\Be2_analysis\\Be2Deconv_results.csv")























































