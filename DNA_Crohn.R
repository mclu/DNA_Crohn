# Analysis of Methylation for Crohn's Disease
#
# Data: GSE138311_series_matrix.txt.gz
# Source: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138311/matrix/
# Description: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138311
#
# Author: Ming-Chen Lu (mingchlu@umich.edu)
# Updated: January 5, 2020
#80: ---------------------------------------------------------------------------
setwd("/Users/Amy/Desktop/DNA_Crohn")

# Libraries: -------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(future)

# Read in the data: ------------------------------------------------------------
mat = fread("matrix.txt", skip = 68)

# Filter to those probes beginning with "ch" and add "sample_group" to indicate
# Crohn's and non-Chrohn's samples
m = mat[grep("^ch", ID_REF)] %>%
  .[, !c("GSM4105199"), with = FALSE] %>%
  melt(., id = 1) %>%
  .[grep("8[7-9]|9[0-3]", variable), `:=` (sample_group = 1L)] %>%
  .[grep("9[4-9]", variable), `:=` (sample_group = 0L)]

# Function to compute intermediate step of t-statistic: ------------------------
tfun = function(mn, sd1){
  # inputs: mn - mean value; sd1 - weighted standard deviation
  # output: a list of values for mean differences and sd for pooled variance
  diff = mn[1] - mn[2]
  sp = sqrt((sd1[1] + sd1[2]) / 10)
  return(list(diff = diff, sp = sp))
}

# Compute t-statistics with pooled variance
t_mat = 
  m[, .(mn = mean(value), s = sd(value)), by = .(ID_REF, sample_group)] %>%
  .[sample_group == 1, sd1 := (6*s^2)] %>%
  .[sample_group == 0, sd1 := (4*s^2)] %>%
  .[, s := NULL] %>%
  .[, c("diff", "sp") := tfun(mn, sd1), by = ID_REF] %>%
  .[, .(tval = diff / (sp*sqrt(1/7 + 1/5))), by = ID_REF] %>%
  unique(.) %>%
  .[, `:=` (probe_group = substr(ID_REF, 1, 5))]

# Plot: ------------------------------------------------------------------------
# Proportion of probes within each probe group that are siginificant at .05 level
plot = 
  t_mat[, pval := pt(-abs(tval), 10)*2] %>%
  .[, .(pp = sum(pval <= .05) / .N), by = probe_group] %>%
  ggplot(., aes(x = probe_group, y = pp)) +
  geom_bar(stat = "identity") +
  ylab("Proportion of probes") + xlab("Probe Group") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Function for permutation tests: ----------------------------------------------
# The function assesses the statistical significance of each probe group
permute_T = 
  function(dt, type = c("two-tailed", "greater", "lesser"), 
           flag = c(TRUE, FALSE), alpha = .05) {
    # inputs:
    #   dt - a data.table in a long format
    #   type - one type of test statistic
    #   flag - a logical flag "permute"
    #   alpha - alpha level for t-distribution
    # outputs: a data.table that stores the values of appropriate t-statistics
    #          for each probe group
    
    if ( flag == FALSE ) {
      result = 
        dt %>%
        .[, .(mn = mean(value), s = sd(value)), by = .(ID_REF, sample_group)] %>%
        .[sample_group == 1, sd1 := (6*s^2)] %>%
        .[sample_group == 0, sd1 := (4*s^2)] %>%
        .[, s := NULL] %>%
        .[, c("diff", "sp") := tfun(mn, sd1), by = ID_REF] %>%
        .[, .(tval = diff / (sp*sqrt(1/7 + 1/5))), by = ID_REF] %>%
        unique(.) %>%
        .[, `:=` (probe_group = substr(ID_REF, 1, 5))] %>%
        .[, `:=` (T_abs = abs(tval) * as.numeric(abs(tval) > qt(1-alpha/2, 10)),
                  T_up = tval * as.numeric(tval > qt(1-alpha, 10)),
                  T_down = tval * as.numeric(tval < qt(alpha, 10)))] %>%
        .[, .(T_abs = mean(T_abs), T_up = mean(T_up), T_down = mean(T_down)), 
          by = probe_group]
      
      if (type == "two-tailed") {
        return(result[, .(probe_group, T_abs)])
      } else if (type == "greater") {
        return(result[, .(probe_group, T_up)])
      } else {
        return(result[, .(probe_group, T_down)])
      }
    } else {
      # permute the sample group labels
      perm_groups = dt[ , .N, .(variable, sample_group)]
      perm_groups[ , sample_group := sample(sample_group, replace = FALSE)]
      
      perm_dt = merge( dt[, !"sample_group"], perm_groups, 
                       by = 'variable')
      
      # compute t-statistics
      result = perm_dt %>%
        .[, .(mn = mean(value), s = sd(value)), by = .(ID_REF, sample_group)] %>%
        .[sample_group == 1, sd1 := (6*s^2)] %>%
        .[sample_group == 0, sd1 := (4*s^2)] %>%
        .[, s := NULL] %>%
        .[, c("diff", "sp") := tfun(mn, sd1), by = ID_REF] %>%
        .[, .(tval = diff / (sp*sqrt(1/7 + 1/5))), by = ID_REF] %>%
        unique(.) %>%
        .[, `:=` (probe_group = substr(ID_REF, 1, 5))] %>%
        .[, `:=` (T_abs = abs(tval) * as.numeric(abs(tval) > qt(1-alpha/2, 10)),
                  T_up = tval * as.numeric(tval > qt(1-alpha, 10)),
                  T_down = tval * as.numeric(tval < qt(alpha, 10)))] %>%
        .[, .(T_abs = mean(T_abs), T_up = mean(T_up), T_down = mean(T_down)), 
          by = probe_group]
      
      if (type == "two-tailed") {
        return(result[, .(probe_group, T_abs)])
      } else if (type == "greater") {
        return(result[, .(probe_group, T_up)])
      } else {
        return(result[, .(probe_group, T_down)])
      }
    }
  }

# Test the function
#permute_T(m, type="two-tailed", flag = TRUE)

# Permute for T_abs, sequentially: ---------------------------------------------
# T_abs score for the original data
tabs = permute_T(m, type = "two-tailed", flag = FALSE, alpha = .05)

# T_abs score for 1000 times permutation
nperm = 1000
t1 = system.time({
  perm_tabs = lapply(1:nperm, function(i) 
    permute_T(m, type = "two-tailed", flag = TRUE, alpha = .05))
})

# Covert list to data.table
tbl_h = as.data.table(matrix(unlist(perm_tabs), nrow = 23))

# drop repeated names of probe group
idx = seq(3, length(tbl_h), by = 2)
tbl_h = tbl_h[, !idx, with = FALSE]

# Merge observed t-score and permuted t-score, then compute p-values
tabs_h = merge(tabs, tbl_h, by.x = "probe_group", by.y = "V1") %>%
  .[, lapply(.SD, function(x) T_abs <= x), 
    by=probe_group, .SDcols = -"T_abs"] %>%
  .[, .(p = {1 + rowSums(.SD)} / {1 + ncol(.SD)}), by = probe_group]

# Permute for T_up, using mclapply: --------------------------------------------
# T_up score for the original data
tup = permute_T(m, type = "greater", flag = FALSE, alpha = .05)

# Compute T_up score using `mclapply` for parallelism
t2 = system.time({
  perm_tup = parallel::mclapply(1:nperm, function(i) 
    permute_T(m, type = "greater", flag = TRUE, alpha = .05))
})

# Covert list to data.table
tbl_i = as.data.table(matrix(unlist(perm_tup), nrow = 23))

# drop repeated names of probe group
idx = seq(3, length(tbl_i), by = 2)
tbl_i = tbl_i[, !idx, with = FALSE]

# Merge observed t-score and permuted t-score, then compute p-values
tup_i = merge(tup, tbl_i, by.x = "probe_group", by.y = "V1") %>%
  .[, lapply(.SD, function(x) T_up <= x), 
    by=probe_group, .SDcols = -"T_up"] %>%
  .[, .(p = {1 + rowSums(.SD)} / {1 + ncol(.SD)}), by = probe_group]

# Permute for T_down, using futures: -------------------------------------------
# T_down score for the original data
tdown = permute_T(m, type = "lesser", flag = FALSE, alpha = .05)

# Compute T_down score using futures for parallelisms
plan("multisession")
t3 = system.time({
  p1 %<-% {lapply(1:200, function(i) 
    permute_T(m, type = "lesser", flag = TRUE, alpha = .05))}
  p2 %<-% {lapply(1:200, function(i) 
    permute_T(m, type = "lesser", flag = TRUE, alpha = .05))}
  p3 %<-% {lapply(1:200, function(i) 
    permute_T(m, type = "lesser", flag = TRUE, alpha = .05))}
  p4 %<-% {lapply(1:200, function(i) 
    permute_T(m, type = "lesser", flag = TRUE, alpha = .05))}
  p5 %<-% {lapply(1:200, function(i) 
    permute_T(m, type = "lesser", flag = TRUE, alpha = .05))}
})
perm_tdown = cbind(p1, p2, p3, p4, p5)

# Covert list to data.table
tbl_j = as.data.table(matrix(unlist(perm_tdown), nrow = 23))

# drop repeated names of probe group
idx = seq(3, length(tbl_j), by = 2)
tbl_j = tbl_j[, !idx, with = FALSE]

# Merge observed t-score and permuted t-score, then compute p-values
tdown_j = merge(tdown, tbl_j, by.x = "probe_group", by.y = "V1") %>%
  .[, lapply(.SD, function(x) T_down <= x), 
    by = probe_group, .SDcols = -"T_down"] %>%
  .[, .(p = {1 + rowSums(.SD)} / {1 + ncol(.SD)}), by = probe_group]

# Plot the result: -------------------------------------------------------------
tbl = merge(tabs_h, tup_i, by = "probe_group") %>%
  merge(., tdown_j, by = "probe_group") %>%
  setnames(., old = c("p.x", "p.y", "p"), new = c("T_abs", "T_up", "T_down"))
