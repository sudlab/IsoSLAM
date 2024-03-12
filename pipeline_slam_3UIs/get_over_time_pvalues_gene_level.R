library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

input_folder =  args[1]

outfile = args[2]

#########################################
### Enter manually for script testing ###
#########################################
#input_folder = "/shared/sudlab1/General/projects/stem_utrons/wet-lab/SLAMseq/MAIN/slam_3UIs/3UI_counts_and_info"

#outfile = "test_output.tsv"
#########################################

files = list.files(input_folder, pattern=".tsv", full.names = TRUE)

data_list = list()

for(file in files){
  file_name = paste(str_split(basename(file), "_")[[1]][1:2], collapse="_")
  
  df = fread(file) %>% 
    group_by(Gene_id, 
             Conversions>=1) %>% 
    summarise(counts = n()) %>%
    mutate(File = file_name)
  
  data_list[[length(data_list) + 1]] = df
  
  rm(df)
}

combined_df = rbindlist(data_list)

rm(data_list)

combined_df = combined_df  %>% 
  separate(File, into=c("Day", "Misc"), sep="_") %>%
  separate(Misc, into=c("Hour", "Rep"), sep="hr") %>%
  group_by(Gene_id,
           Day,
           Hour, 
           Rep) %>% 
  summarise(Converted = `Conversions >= 1`, 
            counts = counts, 
            total_counts = sum(counts)) %>% 
  mutate(pct = counts/total_counts) %>%
  na.omit() # na.omit here removes the no4su samples

non_captured = combined_df %>% 
  group_by(Gene_id, Day, Hour, Rep) %>% 
  summarise(n=n()) %>% 
  filter(n==1) %>% 
  inner_join(combined_df) %>%
  filter(Converted==FALSE)

# We have the events where Converted==FALSE account for 100% of the counts
# Use this to create the corresponding Converted==TRUE where counts=0. Can keep total_counts the same
non_captured$counts = 0
non_captured$Converted = TRUE
non_captured$pct = 0
non_captured = non_captured %>% dplyr::select(-n)

# Add these back onto the combo df 
combined_df = rbind(combined_df, non_captured)

combined_df = combined_df %>% 
  filter(Converted) %>%
  dplyr::select(Gene_id,
                Day,
                Hour,
                Rep,
                Converted_counts = counts,
                Total_counts = total_counts,
                Pct_converted = pct)

named_combined_df = combined_df %>% 
  rowwise() %>% 
  mutate(new_col_name = paste(Hour, Rep, sep="_"))

d0_pct_matrix = named_combined_df %>%
  filter(Day=="d0") %>%
  pivot_wider(id_cols=Gene_id,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("Gene_id")

d0_count_matrix = named_combined_df %>%
  filter(Day=="d0") %>%
  pivot_wider(id_cols=Gene_id,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("Gene_id")

d2_pct_matrix = named_combined_df %>%
  filter(Day=="d2") %>%
  pivot_wider(id_cols=Gene_id,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("Gene_id")

d2_count_matrix = named_combined_df %>%
  filter(Day=="d2") %>%
  pivot_wider(id_cols=Gene_id,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("Gene_id")

d16_pct_matrix = named_combined_df %>%
  filter(Day=="d16") %>%
  pivot_wider(id_cols=Gene_id,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("Gene_id")

d16_count_matrix = named_combined_df %>%
  filter(Day=="d16") %>%
  pivot_wider(id_cols=Gene_id,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("Gene_id")

rm(named_combined_df)

write.csv(as.data.frame(d0_count_matrix), paste0(dirname(outfile),"/","d0_count_matix.csv"))
write.csv(as.data.frame(d2_count_matrix), paste0(dirname(outfile),"/","d2_count_matix.csv"))
write.csv(as.data.frame(d16_count_matrix), paste0(dirname(outfile),"/","d16_count_matix.csv"))
write.csv(as.data.frame(d0_pct_matrix), paste0(dirname(outfile),"/","d0_pct_matix.csv"))
write.csv(as.data.frame(d2_pct_matrix), paste0(dirname(outfile),"/","d2_pct_matix.csv"))
write.csv(as.data.frame(d16_pct_matrix), paste0(dirname(outfile),"/","d16_pct_matix.csv"))

########################################################
#### FUNCTION TO GET P VALUES AND ESTIMATE B VALUES ####
########################################################

get_pval_from_combo_nls_all <- function(gene_to_interogate){
  set.seed(123)
  
  d0_counts =d0_count_matrix[gene_to_interogate,] %>% unlist() 
  d0_pct = d0_pct_matrix[gene_to_interogate,] %>% unlist()
  d0_times = as.numeric(gsub("_.*", "", names(d0_counts)))
  
  d2_counts =d2_count_matrix[gene_to_interogate,] %>% unlist() 
  d2_pct = d2_pct_matrix[gene_to_interogate,] %>% unlist()
  d2_times = as.numeric(gsub("_.*", "", names(d2_counts)))
  
  d16_counts =d16_count_matrix[gene_to_interogate,] %>% unlist() 
  d16_pct = d16_pct_matrix[gene_to_interogate,] %>% unlist()
  d16_times = as.numeric(gsub("_.*", "", names(d16_counts)))
  
  d0_df = data.frame(total_counts = d0_counts,
                     value = d0_pct,
                     time = d0_times,
                     day = "d0")
  
  d2_df = data.frame(total_counts = d2_counts,
                     value = d2_pct,
                     time = d2_times,
                     day = "d2")
  
  d16_df = data.frame(total_counts = d16_counts,
                      value = d16_pct,
                      time = d16_times,
                      day = "d16")
  
  actual_data = rbind(d0_df, d2_df, d16_df)
  
  exponential_func_forced_1 <- function(b, time) {
    #as we are going to force that y intercept to be 1 we do not need an 'a' param
    predicted <- exp(b * time)
    return(predicted)
  }
  
  average0_d0 = actual_data %>% filter(time==0, day=="d0")
  average0_d0_weighted_mean = weighted.mean(average0_d0$value, average0_d0$total_counts)
  
  average0_d2 = actual_data %>% filter(time==0, day=="d2")
  average0_d2_weighted_mean = weighted.mean(average0_d2$value, average0_d2$total_counts)
  
  average0_d16 = actual_data %>% filter(time==0, day=="d16")
  average0_d16_weighted_mean = weighted.mean(average0_d16$value, average0_d16$total_counts)
  
  actual_data_forced_1 = actual_data %>% 
    mutate(value0 = ifelse(day=="d0", average0_d0_weighted_mean, ifelse(day=="d2", average0_d2_weighted_mean, ifelse(day=="d16", average0_d16_weighted_mean, NA))),
           value = value/value0)
  
  # we should also normalize the weighting factor so that it reflects the within-day differences
  # this will be important for examples where there are loads more reads in a given day 
  
  total_counts_d0 = actual_data_forced_1 %>% filter(day=="d0") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_d2 = actual_data_forced_1 %>% filter(day=="d2") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_d16 = actual_data_forced_1 %>% filter(day=="d16") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  
  actual_data_forced_1 = actual_data_forced_1 %>%
    mutate(sum_total_counts = ifelse(day=="d0", total_counts_d0, ifelse(day=="d2", total_counts_d2, ifelse(day=="d16", total_counts_d16, NA))),
           normalized_total_counts = total_counts/sum_total_counts)
  
  actual_data_forced_1$day = factor(actual_data_forced_1$day, levels=c("d0", "d2", "d16"))
  
  tryCatch(
    {  
      d0_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[actual_data_forced_1$day == "d0", ],
                             weights = normalized_total_counts)
      
      d2_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[actual_data_forced_1$day == "d2", ],
                             weights = normalized_total_counts)
      
      d16_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                              start = c(b = 0),
                              data = actual_data_forced_1[actual_data_forced_1$day == "d16", ],
                              weights = normalized_total_counts)
      
      shared_params_fit = nls(value ~ exponential_func_forced_1(b, time),
                              start = c(b=0),
                              data = actual_data_forced_1,
                              weights= normalized_total_counts)
      
      combined_params_fit <- nls(value ~ exponential_func_forced_1(b[day], time),
                                 start = list(b = c(d0 = 0, d2 = 0, d16=0)),
                                 data = actual_data_forced_1,
                                 weights = normalized_total_counts)
      
      x =anova(combined_params_fit, shared_params_fit)
      pval = x$`Pr(>F)`[2]
      d0_b = coef(d0_fit_weighted)[1]
      d2_b = coef(d2_fit_weighted)[1]
      d16_b = coef(d16_fit_weighted)[1]
      hl_d0 = log(2)/-d0_b
      hl_d2 = log(2)/-d2_b
      hl_d16 = log(2)/-d16_b
      output = c(pval, d0_b, d2_b, d16_b, hl_d0, hl_d2, hl_d16)
      names(output) = c("pval", "d0_b", "d2_b", "d16_b", "hl_d0", "hl_d2", "hl_d16")
      return(output)
    }, error = function(e){
      return(NA)
    }
  )
}

##################################################################
#### RUN IT ON THE GENES THAT CAN BE FIT AT ALL 3 TIME POINTS ####
##################################################################

pvals = data.frame(gene=c(), pval=c(), d0_b=c(), d2_b=c(), d16_b=c(), hl_d0=c(), hl_d2=c(), hl_d16=c())

in_all = intersect(intersect(rownames(d0_pct_matrix), rownames(d2_pct_matrix)), rownames(d16_pct_matrix))

for(gene in in_all){
  data = get_pval_from_combo_nls_all(gene)
  pval = data[1]
  d0_b = data[2]
  d2_b = data[3]
  d16_b = data[4]
  hl_d0 = data[5]
  hl_d2 = data[6]
  hl_d16 = data[7]
  pvals = rbind(pvals, data.frame(gene, pval, d0_b, d2_b, d16_b, hl_d0, hl_d2, hl_d16))
}

# Filter out things where the half life estimate is >24 hours, or negative
pvals = na.omit(pvals) %>% filter(hl_d0<24, hl_d2<24, hl_d16<24, hl_d0>0, hl_d2>0, hl_d16>0)

# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, outfile, sep="\t")

#############################################################
#### RUN IT IN A PAIRWISE MANNER 0 vs 2, 0 vs 16, 2 v 16 ####
#############################################################

get_pval_from_combo_nls_pairwise <- function(gene_to_interogate, pair){
  
  comp1 = pair[1]
  comp2 = pair[2]
  
  if(comp1=="d0"){
    comp1_count_matrix = d0_count_matrix
    comp1_pct_matrix = d0_pct_matrix
  }else if(comp1=="d2"){
    comp1_count_matrix = d2_count_matrix
    comp1_pct_matrix = d2_pct_matrix
  }else if(comp1=="d16"){
    comp1_count_matrix = d16_count_matrix
    comp1_pct_matrix = d16_pct_matrix
  }
  
  if(comp2=="d0"){
    comp2_count_matrix = d0_count_matrix
    comp2_pct_matrix = d0_pct_matrix
  }else if(comp2=="d2"){
    comp2_count_matrix = d2_count_matrix
    comp2_pct_matrix = d2_pct_matrix
  }else if(comp2=="d16"){
    comp2_count_matrix = d16_count_matrix
    comp2_pct_matrix = d16_pct_matrix
  }
  
  set.seed(123)
  
  comp1_counts =comp1_count_matrix[gene_to_interogate,] %>% unlist() 
  comp1_pct = comp1_pct_matrix[gene_to_interogate,] %>% unlist()
  comp1_times = as.numeric(gsub("_.*", "", names(comp1_counts)))
  
  comp2_counts =comp2_count_matrix[gene_to_interogate,] %>% unlist() 
  comp2_pct = comp2_pct_matrix[gene_to_interogate,] %>% unlist()
  comp2_times = as.numeric(gsub("_.*", "", names(comp2_counts)))

  comp1_df = data.frame(total_counts = comp1_counts,
                     value = comp1_pct,
                     time = comp1_times,
                     day = "comp1")
  
  comp2_df = data.frame(total_counts = comp2_counts,
                     value = comp2_pct,
                     time = comp2_times,
                     day = "comp2")
  
  actual_data = rbind(comp1_df, comp2_df)
  
  exponential_func_forced_1 <- function(b, time) {
    #as we are going to force that y intercept to be 1 we do not need an 'a' param
    predicted <- exp(b * time)
    return(predicted)
  }
  
  average0_comp1 = actual_data %>% filter(time==0, day=="comp1")
  average0_comp1_weighted_mean = weighted.mean(average0_comp1$value, average0_comp1$total_counts)
  
  average0_comp2 = actual_data %>% filter(time==0, day=="comp2")
  average0_comp2_weighted_mean = weighted.mean(average0_comp2$value, average0_comp2$total_counts)
  
  actual_data_forced_1 = actual_data %>% 
    mutate(value0 = ifelse(day=="comp1", average0_comp1_weighted_mean, ifelse(day=="comp2", average0_comp2_weighted_mean, NA)),
           value = value/value0)
  
  # we should also normalize the weighting factor so that it reflects the within-day differences
  # this will be important for examples where there are loads more reads in a given day 
  
  total_counts_comp1 = actual_data_forced_1 %>% filter(day=="comp1") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_comp2 = actual_data_forced_1 %>% filter(day=="comp2") %>% dplyr::select(total_counts) %>% unlist() %>% sum()

  actual_data_forced_1 = actual_data_forced_1 %>%
    mutate(sum_total_counts = ifelse(day=="comp1", total_counts_comp1, ifelse(day=="comp2", total_counts_comp2, NA)),
           normalized_total_counts = total_counts/sum_total_counts)
  
  actual_data_forced_1$day = factor(actual_data_forced_1$day, levels=c("comp1", "comp2"))
  
  tryCatch(
    {  
      comp1_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[actual_data_forced_1$day == "comp1", ],
                             weights = normalized_total_counts)
      
      comp2_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[actual_data_forced_1$day == "comp2", ],
                             weights = normalized_total_counts)
      
      shared_params_fit = nls(value ~ exponential_func_forced_1(b, time),
                              start = c(b=0),
                              data = actual_data_forced_1,
                              weights= normalized_total_counts)
      
      combined_params_fit <- nls(value ~ exponential_func_forced_1(b[day], time),
                                 start = list(b = c(comp1 = 0, comp2 = 0)),
                                 data = actual_data_forced_1,
                                 weights = normalized_total_counts)
      
      x =anova(combined_params_fit, shared_params_fit)
      pval = x$`Pr(>F)`[2]
      comp1_b = coef(comp1_fit_weighted)[1]
      comp2_b = coef(comp2_fit_weighted)[1]
      hl_comp1 = log(2)/-comp1_b
      hl_comp2 = log(2)/-comp2_b
      output = c(pval, comp1_b, comp2_b, hl_comp1, hl_comp2)
      names(output) = c("pval", "comp1_b", "comp2_b", "hl_comp1", "hl_comp2")
      return(output)
    }, error = function(e){
      return(NA)
    }
  )
}

########################
#### Day 0 vs Day 2 ####
########################

pvals = data.frame(gene=c(), pval=c(), d0_b=c(), d2_b=c(), hl_d0=c(), hl_d2=c())
intersect_0_2 = intersect(rownames(d0_pct_matrix), rownames(d2_pct_matrix))

for(gene in intersect_0_2){
  data = get_pval_from_combo_nls_pairwise(gene, c("d0", "d2"))
  pval = data[1]
  d0_b = data[2]
  d2_b = data[3]
  hl_d0 = data[4]
  hl_d2 = data[5]
  pvals = rbind(pvals, data.frame(gene, pval, d0_b, d2_b, hl_d0, hl_d2))
}
# Filter out things where the half life estimate is >24 hours, or negative
pvals = na.omit(pvals) %>% filter(hl_d0<24, hl_d2<24, hl_d0>0, hl_d2>0)
# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, paste0(dirname(outfile),"/","d0_vs_d2_pvals.tsv"), sep="\t")

########################
#### Day 0 vs Day 16 ####
########################

pvals = data.frame(gene=c(), pval=c(), d0_b=c(), d16_b=c(), hl_d0=c(), hl_d16=c())
intersect_0_16 = intersect(rownames(d0_pct_matrix), rownames(d16_pct_matrix))

for(gene in intersect_0_16){
  data = get_pval_from_combo_nls_pairwise(gene, c("d0", "d16"))
  pval = data[1]
  d0_b = data[2]
  d16_b = data[3]
  hl_d0 = data[4]
  hl_d16 = data[5]
  pvals = rbind(pvals, data.frame(gene, pval, d0_b, d16_b, hl_d0, hl_d16))
}
# Filter out things where the half life estimate is >24 hours, or negative
pvals = na.omit(pvals) %>% filter(hl_d0<24, hl_d16<24, hl_d0>0, hl_d16>0)
# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, paste0(dirname(outfile),"/","d0_vs_d16_pvals.tsv"), sep="\t")


########################
#### Day 2 vs Day 16 ####
########################

pvals = data.frame(gene=c(), pval=c(), d2_b=c(), d16_b=c(), hl_d2=c(), hl_d16=c())
intersect_2_16 = intersect(rownames(d2_pct_matrix), rownames(d16_pct_matrix))

for(gene in intersect_2_16){
  data = get_pval_from_combo_nls_pairwise(gene, c("d2", "d16"))
  pval = data[1]
  d2_b = data[2]
  d16_b = data[3]
  hl_d2 = data[4]
  hl_d16 = data[5]
  pvals = rbind(pvals, data.frame(gene, pval, d2_b, d16_b, hl_d2, hl_d16))
}
# Filter out things where the half life estimate is >24 hours, or negative
pvals = na.omit(pvals) %>% filter(hl_d2<24, hl_d16<24, hl_d2>0, hl_d16>0)
# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, paste0(dirname(outfile),"/","d2_vs_d16_pvals.tsv"), sep="\t")