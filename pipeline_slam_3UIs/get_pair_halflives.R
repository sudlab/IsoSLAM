library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

input_folder =  args[1]

input_regex = args[2]

num_bootstraps = as.numeric(args[3])

outfile = args[4]

#########################################
### Enter manually for script testing ###
#########################################
#input_folder = "/shared/sudlab1/General/projects/stem_utrons/wet-lab/SLAMseq/MAIN/slam_3UIs/3UI_counts_and_info/"

#input_regex = "d2_"

#num_bootstraps = 1000

#outfile = "test_output.tsv"
#########################################

files = list.files(input_folder, pattern=input_regex, full.names = TRUE)

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

### Filter for events we have pvals for (0hr > half-life > 24hr)
pvals = fread(paste0(input_folder,"/pairs/",input_regex,"_pvalues.tsv"))
named_combined_df = named_combined_df %>% filter(eventName %in% pvals$event)

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

###########################
#### FITTING FUNCTIONS ####
###########################

# function to optimize
exponential_func <- function(params, time) {
  a <- params[1] # start point
  b <- params[2] # decay rate (technically rate of change, but we expect it to be negative) 
  predicted <- a * exp(b * time) # predict values at each point based on the time point
  return(predicted)
}

# function to fit
fit_with_nls = function(values,times){
  tryCatch(
    {
      fit = nls(values ~ exponential_func(c(a,b), times),
                start=c(a=1, b=0))
      return(fit)
    }, error = function(e){
      return(NA)
    }
  )
}

#######################
#### BOOTSTRAPPING ####
#######################

set.seed(123)
i = 1

output_df = data.frame(matrix(nrow=0, ncol=9))
colnames(output_df) = c("event","bs_spliced_a", "bs_spliced_a_sd", "bs_retained_a", "bs_retained_a_sd",
                        "bs_spliced_b", "bs_spliced_b_sd", "bs_retained_b","bs_retained_b_sd")

in_both = intersect(rownames(ret_pct_matrix), rownames(spl_pct_matrix))

for(event in in_both){
  i = i+1
  #print(i)
  #if(i>10){
  #  break
  #}
  
  spliced_values = spl_pct_matrix[event,] %>% unlist()
  spliced_times = as.numeric(gsub("_.*", "", names(spliced_values)))
  retained_values = ret_pct_matrix[event,] %>% unlist()
  retained_times = as.numeric(gsub("_.*", "", names(retained_values)))
  
  spliced_counts = spl_count_matrix[event,] %>% unlist()
  retained_counts = ret_count_matrix[event,] %>% unlist()
  
  spliced_fit = fit_with_nls(spliced_values, spliced_times)
  retained_fit = fit_with_nls(retained_values, retained_times)
  
  spliced_sampling_a = c()
  spliced_sampling_b = c()
  retained_sampling_a = c()
  retained_sampling_b = c()
  
  for(bootstrap_i in 1:num_bootstraps){
    #################################
    ###         SPLICED           ###
    #################################
    sampled_spliced_values = c()
    
    for(timepoint in unique(combined_df$Hour)){
      values =  spliced_values[grep(paste0("^",timepoint,"_"), names(spliced_values))]
      weights = spliced_counts[grep(paste0("^",timepoint,"_"), names(spliced_counts))]
      
      #both matrices were made in the same manner therefore the order should be the same
      #but just to be safe, if this is not the case throw an error
      stopifnot(all(names(values)==names(weights)))
      
      weighted_mean = sum(values * weights) / sum(weights)
      weighted_sd = sqrt(sum(weights * (values-weighted_mean)^2) / sum(weights))
      sampled_values = rnorm(length(values), weighted_mean, weighted_sd)
      names(sampled_values) = paste0(timepoint, "_", 1:length(values))
      sampled_spliced_values = c(sampled_spliced_values, sampled_values)
    }
    
    times = as.numeric(gsub("_.*", "", names(sampled_spliced_values)))
    bootstrap_spliced_fit = fit_with_nls(sampled_spliced_values, times)
    
    if(!any(is.na(bootstrap_spliced_fit))){
      spliced_sampling_a = c(spliced_sampling_a, coefficients(bootstrap_spliced_fit)[1])
      spliced_sampling_b = c(spliced_sampling_b, coefficients(bootstrap_spliced_fit)[2])
    }  
    
    ##################################
    ###         RETAINED           ###
    ##################################
    sampled_retained_values = c()
    
    for(timepoint in unique(combined_df$Hour)){
      values =  retained_values[grep(paste0("^",timepoint,"_"), names(retained_values))]
      weights = retained_counts[grep(paste0("^",timepoint,"_"), names(retained_counts))]
      
      #both matrices were made in the same manner therefore the order should be the same
      #but just to be safe, if this is not the case throw an error
      stopifnot(all(names(values)==names(weights)))
      
      weighted_mean = sum(values * weights) / sum(weights)
      weighted_sd = sqrt(sum(weights * (values-weighted_mean)^2) / sum(weights))
      sampled_values = rnorm(length(values), weighted_mean, weighted_sd)
      names(sampled_values) = paste0(timepoint, "_", 1:length(values))
      sampled_retained_values = c(sampled_retained_values, sampled_values)
    }
    
    times = as.numeric(gsub("_.*", "", names(sampled_retained_values)))
    bootstrap_retained_fit = fit_with_nls(sampled_retained_values, times)
    
    if(!any(is.na(bootstrap_retained_fit))){
      retained_sampling_a = c(retained_sampling_a, coefficients(bootstrap_retained_fit)[1])
      retained_sampling_b = c(retained_sampling_b, coefficients(bootstrap_retained_fit)[2])
    }  
    
  } 
  
  
  bs_spliced_a = mean(spliced_sampling_a)
  bs_spliced_a_sd = sd(spliced_sampling_a)
  bs_spliced_b = mean(spliced_sampling_b)
  bs_spliced_b_sd = sd(spliced_sampling_b)
  
  bs_retained_a = mean(retained_sampling_a)
  bs_retained_a_sd = sd(retained_sampling_a)
  bs_retained_b = mean(retained_sampling_b)
  bs_retained_b_sd = sd(retained_sampling_b)
  
  to_add = data.frame(matrix(c(event,
                               bs_spliced_a,
                               bs_spliced_a_sd,
                               bs_retained_a,
                               bs_retained_a_sd,
                               bs_spliced_b,
                               bs_spliced_b_sd,
                               bs_retained_b,
                               bs_retained_b_sd),nrow=1, ncol=9))
  colnames(to_add) = c("event","bs_spliced_a", "bs_spliced_a_sd", "bs_retained_a", "bs_retained_a_sd",
                       "bs_spliced_b", "bs_spliced_b_sd", "bs_retained_b","bs_retained_b_sd")
  output_df = rbind(output_df, to_add)
  
}

output_df = output_df %>%
  mutate(event = as.character(event),
         bs_spliced_a = as.numeric(as.character(bs_spliced_a)),
         bs_spliced_a_sd = as.numeric(as.character(bs_spliced_a_sd)),
         bs_retained_a = as.numeric(as.character(bs_retained_a)),
         bs_retained_a_sd = as.numeric(as.character(bs_retained_a_sd)),
         bs_spliced_b = as.numeric(as.character(bs_spliced_b)),
         bs_spliced_b_sd = as.numeric(as.character(bs_spliced_b_sd)),
         bs_retained_b = as.numeric(as.character(bs_retained_b)),
         bs_retained_b_sd = as.numeric(as.character(bs_retained_b_sd)))

fwrite(output_df, outfile, sep="\t")
