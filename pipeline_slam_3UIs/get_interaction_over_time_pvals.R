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
#input_folder = "../3UI_counts_and_info"

#outfile = "test_output.tsv"
#########################################

files = list.files(input_folder, pattern=".tsv", full.names = TRUE)

data_list = list()

for(file in files){
  file_name = paste(str_split(basename(file), "_")[[1]][1:2], collapse="_")
  
  df = fread(file) %>% 
    group_by(Transcript_id, 
             Chr,
             Strand,
             Start, 
             End, 
             Assignment, 
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
  group_by(Transcript_id, 
           Chr, 
           Strand,
           Start, 
           End, 
           Assignment,
           Day,
           Hour, 
           Rep) %>% 
  summarise(Converted = `Conversions >= 1`, 
            counts = counts, 
            total_counts = sum(counts)) %>% 
  mutate(pct = counts/total_counts) %>%
  na.omit() # na.omit here removes the no4su samples

non_captured = combined_df %>% 
  group_by(Transcript_id, Chr, Strand, Start, End, Assignment, Day, Hour, Rep) %>% 
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
  dplyr::select(Transcript_id,
                Chr, 
                Strand,
                Start,
                End,
                Assignment,
                Day,
                Hour,
                Rep,
                Converted_counts = counts,
                Total_counts = total_counts,
                Pct_converted = pct)

named_combined_df = combined_df %>% 
  rowwise() %>% 
  mutate(new_col_name = paste(Hour, Rep, sep="_")) %>%
  mutate(eventName = paste(Transcript_id, Chr, Strand, Start, End, sep=":"))

rm(combined_df)

d0_ret_pct_matrix = named_combined_df %>%
  filter(Day=="d0", Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d0_ret_count_matrix = named_combined_df %>%
  filter(Day=="d0", Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d0_spl_pct_matrix = named_combined_df %>%
  filter(Day=="d0", Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d0_spl_count_matrix = named_combined_df %>%
  filter(Day=="d0", Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")
  
d2_ret_pct_matrix = named_combined_df %>%
  filter(Day=="d2", Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d2_ret_count_matrix = named_combined_df %>%
  filter(Day=="d2", Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
             names_from=new_col_name,
             values_from=Total_counts,
             values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")
  
d2_spl_pct_matrix = named_combined_df %>%
  filter(Day=="d2", Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d2_spl_count_matrix = named_combined_df %>%
  filter(Day=="d2", Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d16_ret_pct_matrix = named_combined_df %>%
  filter(Day=="d16", Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")
  
d16_ret_count_matrix = named_combined_df %>%
  filter(Day=="d16", Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d16_spl_pct_matrix = named_combined_df %>%
  filter(Day=="d16", Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

d16_spl_count_matrix = named_combined_df %>%
  filter(Day=="d16", Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

rm(named_combined_df)


########################################################
#### FUNCTION TO GET P VALUES AND ESTIMATE B VALUES ####
########################################################

get_pval_from_combo_nls <- function(event_to_interogate){
  set.seed(123)
  
  d0_ret_counts =d0_ret_count_matrix[event_to_interogate,] %>% unlist() 
  d0_ret_pct = d0_ret_pct_matrix[event_to_interogate,] %>% unlist()
  d0_ret_times = as.numeric(gsub("_.*", "", names(d0_ret_counts)))
  
  d2_ret_counts =d2_ret_count_matrix[event_to_interogate,] %>% unlist() 
  d2_ret_pct = d2_ret_pct_matrix[event_to_interogate,] %>% unlist()
  d2_ret_times = as.numeric(gsub("_.*", "", names(d2_ret_counts)))
  
  d16_ret_counts =d16_ret_count_matrix[event_to_interogate,] %>% unlist() 
  d16_ret_pct = d16_ret_pct_matrix[event_to_interogate,] %>% unlist()
  d16_ret_times = as.numeric(gsub("_.*", "", names(d16_ret_counts)))
  
  d0_spl_counts =d0_spl_count_matrix[event_to_interogate,] %>% unlist() 
  d0_spl_pct = d0_spl_pct_matrix[event_to_interogate,] %>% unlist()
  d0_spl_times = as.numeric(gsub("_.*", "", names(d0_spl_counts)))
  
  d2_spl_counts =d2_spl_count_matrix[event_to_interogate,] %>% unlist() 
  d2_spl_pct = d2_spl_pct_matrix[event_to_interogate,] %>% unlist()
  d2_spl_times = as.numeric(gsub("_.*", "", names(d2_spl_counts)))
  
  d16_spl_counts =d16_spl_count_matrix[event_to_interogate,] %>% unlist() 
  d16_spl_pct = d16_spl_pct_matrix[event_to_interogate,] %>% unlist()
  d16_spl_times = as.numeric(gsub("_.*", "", names(d16_spl_counts)))
  
  d0_ret_df = data.frame(total_counts = d0_ret_counts,
                          value = d0_ret_pct,
                          time = d0_ret_times,
                          day = "d0",
                          isoform="ret")
  
  d2_ret_df = data.frame(total_counts = d2_ret_counts,
                     value = d2_ret_pct,
                     time = d2_ret_times,
                     day = "d2",
                     isoform="ret")
  
  d16_ret_df = data.frame(total_counts = d16_ret_counts,
                     value = d16_ret_pct,
                     time = d16_ret_times,
                     day = "d16",
                     isoform="ret")
  
  d0_spl_df = data.frame(total_counts = d0_spl_counts,
                         value = d0_spl_pct,
                         time = d0_spl_times,
                         day = "d0",
                         isoform="spl")
  
  d2_spl_df = data.frame(total_counts = d2_spl_counts,
                         value = d2_spl_pct,
                         time = d2_spl_times,
                         day = "d2",
                         isoform="spl")
  
  d16_spl_df = data.frame(total_counts = d16_spl_counts,
                          value = d16_spl_pct,
                          time = d16_spl_times,
                          day = "d16",
                          isoform="spl")
  
  actual_data = rbind(d0_ret_df, d2_ret_df, d16_ret_df, d0_spl_df, d2_spl_df, d16_spl_df)
  
  exponential_func_forced_1 <- function(b, time) {
    #as we are going to force that y intercept to be 1 we do not need an 'a' param
    predicted <- exp(b * time)
    return(predicted)
  }
  
  average0_ret_d0 = actual_data %>% filter(time==0, day=="d0", isoform=="ret")
  average0_ret_d0_weighted_mean = weighted.mean(average0_ret_d0$value, average0_ret_d0$total_counts)
  
  average0_ret_d2 = actual_data %>% filter(time==0, day=="d2", isoform=="ret")
  average0_ret_d2_weighted_mean = weighted.mean(average0_ret_d2$value, average0_ret_d2$total_counts)
  
  average0_ret_d16 = actual_data %>% filter(time==0, day=="d16", isoform=="ret")
  average0_ret_d16_weighted_mean = weighted.mean(average0_ret_d16$value, average0_ret_d16$total_counts)
  
  average0_spl_d0 = actual_data %>% filter(time==0, day=="d0", isoform=="spl")
  average0_spl_d0_weighted_mean = weighted.mean(average0_spl_d0$value, average0_spl_d0$total_counts)
  
  average0_spl_d2 = actual_data %>% filter(time==0, day=="d2", isoform=="spl")
  average0_spl_d2_weighted_mean = weighted.mean(average0_spl_d2$value, average0_spl_d2$total_counts)
  
  average0_spl_d16 = actual_data %>% filter(time==0, day=="d16", isoform=="spl")
  average0_spl_d16_weighted_mean = weighted.mean(average0_spl_d16$value, average0_spl_d16$total_counts)
  
  actual_data_forced_1 = actual_data %>% 
    mutate(value0 = ifelse((day=="d0"&isoform=="ret"), average0_ret_d0_weighted_mean, 
                           ifelse((day=="d0"&isoform=="spl"), average0_spl_d0_weighted_mean,
                                  ifelse((day=="d2"&isoform=="ret"), average0_ret_d2_weighted_mean,
                                         ifelse((day=="d2"&isoform=="spl"), average0_spl_d2_weighted_mean,
                                                ifelse((day=="d16"&isoform=="ret"), average0_ret_d16_weighted_mean,
                                                       ifelse((day=="d16"&isoform=="spl"), average0_spl_d16_weighted_mean,
                                                              NA)))))),
           value = value/value0)
  
  # we should also normalize the weighting factor so that it reflects the within-day-within-isoform differences
  # this will be important for examples where there are loads more reads in a given day or a given isoform
  # so the joint model is not just skewed by a day or isoform having much higher expression
  
  total_counts_ret_d0 = actual_data_forced_1 %>% filter(day=="d0", isoform=="ret") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_ret_d2 = actual_data_forced_1 %>% filter(day=="d2", isoform=="ret") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_ret_d16 = actual_data_forced_1 %>% filter(day=="d16", isoform=="ret") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_spl_d0 = actual_data_forced_1 %>% filter(day=="d0", isoform=="spl") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_spl_d2 = actual_data_forced_1 %>% filter(day=="d2", isoform=="spl") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_spl_d16 = actual_data_forced_1 %>% filter(day=="d16", isoform=="spl") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  
  actual_data_forced_1 = actual_data_forced_1 %>%
    mutate(sum_total_counts = ifelse((day=="d0"&isoform=="ret"), total_counts_ret_d0, 
                                     ifelse((day=="d2"&isoform=="ret"), total_counts_ret_d2, 
                                            ifelse((day=="d16"&isoform=="ret"), total_counts_ret_d16, 
                                                   ifelse((day=="d0"&isoform=="spl"), total_counts_spl_d0,
                                                          ifelse((day=="d2"&isoform=="spl"), total_counts_spl_d2,
                                                                 ifelse((day=="d16"&isoform=="spl"), total_counts_spl_d16,
                                                                        NA)))))),
           normalized_total_counts = total_counts/sum_total_counts)
  
  actual_data_forced_1$day = factor(actual_data_forced_1$day, levels=c("d0", "d2", "d16"))
  
  tryCatch(
    { 
      ################################################################################
      ### INDEPENDENT FITS - I.E. D0 RET, D0 SPL, D2 RET, D2 SPL, D16 RET, D16 SPL ###
      ################################################################################
      
      d0_ret_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                  start = c(b = 0),
                                  data = actual_data_forced_1[(actual_data_forced_1$day=="d0"&actual_data_forced_1$isoform=="ret"), ],
                                  weights = normalized_total_counts)
 
      d2_ret_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[(actual_data_forced_1$day == "d2"&actual_data_forced_1$isoform=="ret"), ],
                             weights = normalized_total_counts)
      
      d16_ret_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[(actual_data_forced_1$day=="d16"&actual_data_forced_1$isoform=="ret"), ],
                             weights = normalized_total_counts)
      
      d0_spl_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                 start = c(b = 0),
                                 data = actual_data_forced_1[(actual_data_forced_1$day=="d0"&actual_data_forced_1$isoform=="spl"), ],
                                 weights = normalized_total_counts)
      
      d2_spl_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                 start = c(b = 0),
                                 data = actual_data_forced_1[(actual_data_forced_1$day == "d2"&actual_data_forced_1$isoform=="spl"), ],
                                 weights = normalized_total_counts)
      
      d16_spl_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                  start = c(b = 0),
                                  data = actual_data_forced_1[(actual_data_forced_1$day=="d16"&actual_data_forced_1$isoform=="spl"), ],
                                  weights = normalized_total_counts)
      
      ########################################
      ### COMBINED FITS - I.E. D0, D2, D16 ###
      ########################################
      
      d0_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                 start = c(b = 0),
                                 data = actual_data_forced_1[actual_data_forced_1$day=="d0", ],
                                 weights = normalized_total_counts)
      
      d2_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[actual_data_forced_1$day=="d2", ],
                             weights = normalized_total_counts)
      
      d16_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                             start = c(b = 0),
                             data = actual_data_forced_1[actual_data_forced_1$day=="d16", ],
                             weights = normalized_total_counts)
      
      ####################
      ### JOINT MODELS ###
      ####################
      
      
      ### FIXED OFFSET ###
      ### This is like plotting 2 lines for each differentiation day based on the average fit for both isoforms
      ### + a fixed offset (which is derived from the average fold difference between isoforms) across
      ### all time points. I.e. the difference between the 2 isoforms is always relatively constant (it is
      ### a fold change, not a fixed number, therefore stays relative to the half-life at that day)
      ### e.g. if the AVERAGE retained isoform is twice as stable as the spliced isoform then here:
      ### if: average d0=6hrs, d2 = 3hrs, d16=12hrs
      ### then: d0_spl=8hrs, d0_ret=4hrs, d2_spl=2hrs, d2_ret=4hrs, d16_spl=8hrs, d16_ret=16hrs
      ####################
      
      d0_fold_b_diff = coef(d0_ret_fit_weighted)/coef(d0_spl_fit_weighted)
      d2_fold_b_diff = coef(d2_ret_fit_weighted)/coef(d2_spl_fit_weighted)
      d16_fold_b_diff = coef(d16_ret_fit_weighted)/coef(d16_spl_fit_weighted)
      mean_b_diff_ratio = mean(d0_fold_b_diff, d2_fold_b_diff, d16_fold_b_diff)
      
      fixed_offset_fit = nls(value ~ as.numeric(isoform=="ret")*as.numeric(day=="d0")*exp((d0_b*mean_b_diff_ratio) * time) + 
                               as.numeric(isoform=="spl")*as.numeric(day=="d0")*exp((d0_b/mean_b_diff_ratio) * time) + 
                               as.numeric(isoform=="ret")*as.numeric(day=="d2")*exp((d2_b*mean_b_diff_ratio) * time) + 
                               as.numeric(isoform=="spl")*as.numeric(day=="d2")*exp((d2_b/mean_b_diff_ratio) * time) + 
                               as.numeric(isoform=="ret")*as.numeric(day=="d16")*exp((d16_b*mean_b_diff_ratio) * time) + 
                               as.numeric(isoform=="spl")*as.numeric(day=="d16")*exp((d16_b/mean_b_diff_ratio) * time),
                             start=list(d0_b = coef(d0_fit_weighted),
                                        d2_b = coef(d2_fit_weighted),
                                        d16_b = coef(d16_fit_weighted)),
                             data = actual_data_forced_1,
                             weights = normalized_total_counts)
      
      ### FREE OFFSET ###
      ### We then do a similar thing, except we don't take the average b difference from all days, we just use the 
      ### difference at that particular differentiation day. 
      ### Hence if there is a flip at a given day, or the difference between isoforms becomes significantly
      ### more or less profound, then an anova between the 2 models would be significant. 
      ### If the difference between isoforms is consistent between days, then the anova will be n.s. 
      ###################
      
      free_offset_fit = nls(value ~ as.numeric(isoform=="ret")*as.numeric(day=="d0")*exp((d0_b*d0_offset) * time) + 
                            as.numeric(isoform=="spl")*as.numeric(day=="d0")*exp((d0_b/d0_offset) * time) + 
                            as.numeric(isoform=="ret")*as.numeric(day=="d2")*exp((d2_b*d2_offset) * time) + 
                            as.numeric(isoform=="spl")*as.numeric(day=="d2")*exp((d2_b/d2_offset) * time) + 
                            as.numeric(isoform=="ret")*as.numeric(day=="d16")*exp((d16_b*d16_offset) * time) + 
                            as.numeric(isoform=="spl")*as.numeric(day=="d16")*exp((d16_b/d16_offset) * time),
                          start=list(d0_b = coef(d0_fit_weighted),
                                     d0_offset = coef(d0_ret_fit_weighted)/coef(d0_spl_fit_weighted),
                                     d2_b = coef(d2_fit_weighted),
                                     d2_offset = coef(d2_ret_fit_weighted)/coef(d2_spl_fit_weighted),
                                     d16_b = coef(d16_fit_weighted),
                                     d16_offset = coef(d16_ret_fit_weighted)/coef(d16_spl_fit_weighted)),
                          data = actual_data_forced_1,
                          weights = normalized_total_counts)

      
      x =anova(fixed_offset_fit, free_offset_fit)
      
      pval = x$`Pr(>F)`[2]
      d0_ret_b = coef(d0_ret_fit_weighted)[1]
      d2_ret_b = coef(d2_ret_fit_weighted)[1]
      d16_ret_b = coef(d16_ret_fit_weighted)[1]
      d0_spl_b = coef(d0_spl_fit_weighted)[1]
      d2_spl_b = coef(d2_spl_fit_weighted)[1]
      d16_spl_b = coef(d16_spl_fit_weighted)[1]
      hl_d0_ret = log(2)/-d0_ret_b
      hl_d2_ret = log(2)/-d2_ret_b
      hl_d16_ret = log(2)/-d16_ret_b
      hl_d0_spl = log(2)/-d0_spl_b
      hl_d2_spl = log(2)/-d2_spl_b
      hl_d16_spl = log(2)/-d16_spl_b
      output = c(pval, d0_ret_b, d0_spl_b, d2_ret_b, d2_spl_b, d16_ret_b, d16_spl_b, hl_d0_ret, hl_d0_spl, hl_d2_ret, hl_d2_spl, hl_d16_ret, hl_d16_spl)
      names(output) = c("pval", "d0_ret_b", "d0_spl_b", "d2_ret_b", "d2_spl_b", "d16_ret_b", "d16_spl_b", "hl_d0_ret", "hl_d0_spl", "hl_d2_ret", "hl_d2_spl", "hl_d16_ret", "hl_d16_spl")
      return(output)
    }, error = function(e){
      return(NA)
    }
  )
}

##############################
#### RUN IT ON THE EVENTS ####
##############################

pvals = data.frame(event=c(), pval=c(), 
                   d0_ret_b=c(), d0_spl_b=c(),
                   d2_ret_b=c(), d2_spl_b=c(), 
                   d16_ret_b=c(), d16_spl_b=c(),
                   hl_d0_ret=c(), hl_d0_spl=c(), 
                   hl_d2_ret=c(), hl_d2_spl=c(),
                   hl_d16_ret=c(), hl_d16_spl=c())

in_all = intersect(intersect(intersect(intersect(intersect(rownames(d0_ret_pct_matrix), rownames(d2_ret_pct_matrix)), rownames(d16_ret_pct_matrix)),
                   rownames(d0_spl_pct_matrix)), rownames(d2_spl_pct_matrix)), rownames(d16_spl_pct_matrix))

for(event in in_all){
  data = get_pval_from_combo_nls(event)
  pval = data[1]
  d0_ret_b = data[2]
  d0_spl_b = data[3]
  d2_ret_b = data[4]
  d2_spl_b = data[5]
  d16_ret_b = data[6]
  d16_spl_b = data[7]
  hl_d0_ret = data[8]
  hl_d0_spl = data[9]
  hl_d2_ret = data[10]
  hl_d2_spl = data[11]
  hl_d16_ret = data[12]
  hl_d16_spl = data[13]
  pvals = rbind(pvals, data.frame(event, pval, d0_ret_b, d0_spl_b, d2_ret_b, d2_spl_b, d16_ret_b, d16_spl_b, hl_d0_ret, hl_d0_spl, hl_d2_ret, hl_d2_spl, hl_d16_ret, hl_d16_spl))
}

# Filter out things where the half life estimate is >24 hours, or negative
pvals = na.omit(pvals) %>% filter(hl_d0_ret<24, hl_d0_spl<24,
                                  hl_d2_ret<24, hl_d2_spl<24,
                                  hl_d16_ret<24, hl_d16_spl<24,
                                  hl_d0_ret>0, hl_d0_spl>0,
                                  hl_d2_ret>0, hl_d2_spl>0,
                                  hl_d16_ret>0, hl_d16_spl>0)

# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, outfile, sep="\t")
