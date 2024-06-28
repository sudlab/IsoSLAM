library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

input_file =  args[1]

outfile = args[2]

#########################################
### Enter manually for script testing ###
#########################################
#input_file = "/shared/sudlab1/General/projects/stem_utrons/wet-lab/SLAMseq/MAIN/slam_3UIs/all_introns_counts_and_info/summarized/d2_summarized.tsv"

#outfile = "test_output.tsv"
#########################################

combined_df = fread(input_file)

combined_df = combined_df  %>% 
  separate(File, into=c("Day", "Misc"), sep="_") %>%
  dplyr::select(-Day) %>%
  separate(Misc, into=c("Hour", "Rep"), sep="hr") %>%
  group_by(Transcript_id, 
           Chr, 
           Strand,
           Start, 
           End, 
           Assignment,
           Hour, 
           Rep) %>% 
  summarise(Converted = `Conversions >= 1`, 
            counts = counts, 
            total_counts = sum(counts)) %>% 
  mutate(pct = counts/total_counts) %>%
  na.omit() # na.omit here removes the no4su samples

non_captured = combined_df %>% 
  group_by(Transcript_id, Chr, Strand, Start, End, Assignment, Hour, Rep) %>% 
  summarise(n=n(),
            Converted=dplyr::first(Converted)) %>% 
  filter(n==1) %>%
  filter(Converted==FALSE) %>%
  inner_join(combined_df) 

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
                Hour,
                Rep,
                Converted_counts = counts,
                Total_counts = total_counts,
                Pct_converted = pct)

########################################################
#### DECIDE WHETHER WE FILTER OR GIVE IT ALL EVENTS ####
#### EITHER GIVE NAMED_OVER10 OR NAMED_COMBINED_DF  ####
########################################################

named_combined_df = combined_df %>% 
  rowwise() %>% 
  mutate(new_col_name = paste(Hour, Rep, sep="_")) %>%
  mutate(eventName = paste(Transcript_id, Chr, Strand, Start, End, sep=":"))

ret_pct_matrix = named_combined_df %>%
  filter(Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

ret_count_matrix = named_combined_df %>%
  filter(Assignment=="Ret") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

spl_pct_matrix = named_combined_df %>%
  filter(Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

spl_count_matrix = named_combined_df%>%
  filter(Assignment=="Spl") %>%
  pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) %>%
  column_to_rownames("eventName")

rm(named_combined_df)

########################################################
#### FUNCTION TO GET PVALUES AND ESTIMATED B VALUES ####
########################################################

get_pval_from_combo_nls <- function(event_to_interogate){
  set.seed(123)
  
  spliced_counts =spl_count_matrix[event_to_interogate,] %>% unlist() 
  spliced_pct = spl_pct_matrix[event_to_interogate,] %>% unlist()
  spliced_times = as.numeric(gsub("_.*", "", names(spliced_counts)))
  retained_counts = ret_count_matrix[event_to_interogate,] %>% unlist()
  retained_pct = ret_pct_matrix[event_to_interogate,] %>% unlist()
  retained_times = as.numeric(gsub("_.*", "", names(retained_counts)))
  
  spliced_df = data.frame(total_counts = spliced_counts,
                          value = spliced_pct,
                          time = spliced_times,
                          isoform = "Spliced")
  retained_df = data.frame(total_counts = retained_counts,
                           value = retained_pct,
                           time = retained_times,
                           isoform = "Retained")
  actual_data = rbind(spliced_df, retained_df)
  
  exponential_func_forced_1 <- function(b, time) {
    #as we are going to force that y intercept to be 1 we do not need an 'a' param
    predicted <- exp(b * time)
    return(predicted)
  }
  
  average0_spliced = actual_data %>% filter(time==0, isoform=="Spliced")
  average0_spliced_weighted_mean = weighted.mean(average0_spliced$value, average0_spliced$total_counts)
  average0_retained = actual_data %>% filter(time==0, isoform=="Retained")
  average0_retained_weighted_mean = weighted.mean(average0_retained$value, average0_retained$total_counts)
  
  actual_data_forced_1 = actual_data %>% 
    mutate(value0 = ifelse(isoform=="Spliced", average0_spliced_weighted_mean, average0_retained_weighted_mean),
           value = value/value0)
  
  # we should also normalize the weighting factor so that it reflects the within-isoform differences
  # this will be important for examples where there are loads more spliced reads than retained, or vv
  
  total_counts_spliced = actual_data_forced_1 %>% filter(isoform=="Spliced") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_retained = actual_data_forced_1 %>% filter(isoform=="Retained") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  actual_data_forced_1 = actual_data_forced_1 %>%
    mutate(sum_total_counts = ifelse(isoform=="Spliced", total_counts_spliced, total_counts_retained),
           normalized_total_counts = total_counts/sum_total_counts)
  actual_data_forced_1$isoform = factor(actual_data_forced_1$isoform, levels=c("Spliced", "Retained"))
  
  tryCatch(
    {  
      spliced_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                  start = c(b = 0),
                                  data = actual_data_forced_1[actual_data_forced_1$isoform == "Spliced", ],
                                  weights = normalized_total_counts)
      
      retained_fit_weighted <- nls(value ~ exponential_func_forced_1(b, time),
                                   start = c(b = 0),
                                   data = actual_data_forced_1[actual_data_forced_1$isoform == "Retained", ],
                                   weights = normalized_total_counts)
      
      shared_params_fit = nls(value ~ exponential_func_forced_1(b, time),
                              start = c(b=0),
                              data = actual_data_forced_1,
                              weights= normalized_total_counts)
      
      combined_params_fit <- nls(value ~ exponential_func_forced_1(b[isoform], time),
                                 start = list(b = c(Spliced = 0, Retained = 0)),
                                 data = actual_data_forced_1,
                                 weights = normalized_total_counts)
      x =anova(combined_params_fit, shared_params_fit)
      pval = x$`Pr(>F)`[2]
      spl_b = coef(spliced_fit_weighted)[1]
      ret_b = coef(retained_fit_weighted)[1]
      hl_spl = log(2)/-spl_b
      hl_ret = log(2)/-ret_b
      output = c(pval, spl_b, ret_b, hl_spl, hl_ret)
      names(output) = c("pval", "spl_b", "ret_b", "hl_spl", "hl_ret")
      return(output)
    }, error = function(e){
      return(NA)
    }
  )
}

##############################
#### RUN IT ON THE EVENTS ####
##############################

pvals = data.frame(event=c(), pval=c(), spl_b=c(), ret_b=c(), hl_spl=c(), hl_ret=c())

in_both = intersect(rownames(ret_pct_matrix), rownames(spl_pct_matrix))

for(event in in_both){
  data = get_pval_from_combo_nls(event)
  pval = data[1]
  spl_b = data[2]
  ret_b = data[3]
  hl_spl = data[4]
  hl_ret = data[5]
  pvals = rbind(pvals, data.frame(event, pval, spl_b, ret_b, hl_spl, hl_ret))
}

# Filter out things where the half life estimate is >24 hours, or negative
pvals = na.omit(pvals) %>% filter(hl_spl<24, hl_ret<24, hl_spl>0, hl_ret>0)

# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, outfile, sep="\t")
