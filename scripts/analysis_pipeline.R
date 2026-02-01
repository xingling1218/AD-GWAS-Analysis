# ==============================================================================
# Project: GWAS Visualisation of Alzheimer's Disease (Wightman et al. 2021)
# Author: Xingling Wan (MSc Statistical Data Science, UoB)
# Date: 2026-01-25
# Purpose: Generate Manhattan and Q-Q plots from 12.69M variants while 
#          maintaining low memory footprint.
# ==============================================================================

library(data.table)
library(tidyverse)

# --- Environment Setup ---
if (!dir.exists("results")) dir.create("results")
# Please ensure the raw GWAS summary statistics are placed in the 'data/raw/' directory.
# Source: Wightman et al. (2021)
input_path <- "data/raw/PGCALZ2sumstatsExcluding23andMe.txt.gz"

# 1.read the data
sumstats <- fread(input_path)

# 2.view the data
head(sumstats)
tail(sumstats)

# --- Data Analysis ---
# 1.Filter out possibly existing invalid data with NA P-values
ad_clean <- sumstats %>% filter(!is.na(p) & p>0)

# 2.release the allocated memory
rm(sumstats)
gc()

# 3.Extract the significant results (P < 0.01)
sig_data <- ad_clean %>% filter(p < 0.01)

# 4.Sample the non-significant “background” (P ≥ 0.01)
nonsig_data <- ad_clean %>% 
  filter(p >= 0.01) %>% 
  sample_frac(0.005) 

print(colnames(ad_clean))
# 5.Combine the two parts
plot_data <- bind_rows(sig_data, nonsig_data) %>%
  select(CHR = chr, BP = PosGRCh37, P = p) %>%
  mutate(SNP = paste0(CHR, ":", BP)) # Use genomic coordinates as SNP IDs

# 6.Release the allocated memory
rm(ad_clean, sig_data, nonsig_data)
gc()

# 7.Check the current row count
nrow(plot_data)

library(qqman)

# --- Visualisation ---
# 1.Draw the plot
png("results/AD_Wightman_Manhattan.png", width = 12, height = 6, units = "in", res = 300)
manhattan(plot_data, 
          main = "Manhattan Plot: Alzheimer's Disease (Wightman 2021)",
          col = c("midnightblue", "skyblue"), 
          suggestiveline = -log10(1e-5), 
          genomewideline = -log10(5e-8))
dev.off()

# 2.Reload and uniformly subsample (1%)
qq_sample <- fread("data/raw/PGCALZ2sumstatsExcluding23andMe.txt.gz") %>%
  filter(!is.na(p)) %>%
  sample_frac(0.01) 

# 3.Release the allocated memory
gc() 

# 4.Compute Lambda
chisq <- qchisq(1 - qq_sample$p, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
print(paste("Genomic Inflation Factor (lambda):", round(lambda, 3)))

# 5.Plot and save the Q-Q plot
png("results/AD_Wightman_QQ.png", width = 6, height = 6, units = "in", res = 300)
qq(qq_sample$p, 
   main = paste("Q-Q Plot: AD (Wightman 2021), lambda =", round(lambda, 3)),
   col = "steelblue", 
   cex = 0.5)
dev.off()

