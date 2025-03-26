
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

input_folder =  args[1]
input_folder <- "tests/resources/tsv/output"
input_regex = args[2]
input_regex <- "*.tsv$"
outfile = args[3]

#########################################
### Enter manually for script testing ###
#########################################
#input_folder = "/shared/sudlab1/General/projects/stem_utrons/wet-lab/SLAMseq/MAIN/slam_3UIs/3UI_counts_and_info/"

#input_regex = "d2_"

#outfile = "test_output.tsv"
#########################################

files = list.files(input_folder, pattern=input_regex, full.names = TRUE)

data_list = list()
for(file in files){
# for(file in files[1:3]){
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

summary_combined_df = combined_df  %>%
  separate(File, into=c("Day", "Misc"), sep="_") %>%
  # dplyr::select(-Day) %>%
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
    mutate(pct = counts/total_counts,
           Rep = as.numeric(stringr::str_replace(Rep, "\\.tsv", "")),
           Day = as.numeric(stringr::str_replace(Day, "d", "")),
           Hour = as.numeric(Hour),
           ) %>%
  na.omit() # na.omit here removes the no4su samples
readr::write_csv(summary_combined_df, "summary_combined.csv")

non_captured = summary_combined_df %>%
    group_by(
        Transcript_id,
        Chr,
        Strand, Start, End, Assignment, Day, Hour, Rep
    ) %>%
    summarise(
        n = n(),
        Converted = dplyr::first(Converted)
    ) %>%
    filter(n == 1) %>%
    filter(Converted == FALSE) |>
    ungroup() |>
    select(Transcript_id, Strand, Start, End, Assignment, Day, Hour, Rep, n, Converted) # %>%
readr::write_csv(non_captured, "filtered_no_conversions.csv")
print("non-captured")
non_captured |> nrow()
non_captured |> head()
print("summary_combined_df")
summary_combined_df |> nrow()
summary_combined_df |> head()
non_captured_joined <- non_captured |>
  inner_join(summary_combined_df)
print("non_captured_joined (post inner_join())")
non_captured_joined |> nrow()
non_captured_joined |> head()
non_captured_joined$Converted |> table()

# We have the events where Converted==FALSE account for 100% of the counts
# Use this to create the corresponding Converted==TRUE where counts=0. Can keep total_counts the same
non_captured_joined$counts = 0
non_captured_joined$Converted = TRUE
non_captured_joined$pct = 0
non_captured_joined = non_captured_joined %>% dplyr::select(-n)
readr::write_csv(non_captured_joined, "inner_join.csv")
nrow(non_captured_joined)
readr::write_csv(non_captured_joined, "non_captured_joined_r.csv")

## Alternatively we can...
##   1. Select the variables of interest
##   2. Create a set where there is only 1 instance False

# Add these back onto the combo df
summary_combined_df2 = rbind(summary_combined_df, non_captured_joined)
nrow(summary_combined_df2)

# summary_combined_df = summary_combined_df %>%
clean = summary_combined_df2 %>%
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
        Pct_converted = pct
        ) %>%
    dplyr::arrange(Transcript_id,
        Chr,
        Strand,
        Start,
        End,
        Assignment,
        Day,
        Hour,
        Rep)
nrow(clean)

readr::write_csv(clean, "final.csv")

########################################################
#### DECIDE WHETHER WE FILTER OR GIVE IT ALL EVENTS ####
#### EITHER GIVE NAMED_OVER10 OR NAMED_SUMMARY_COMBINED_DF  ####
########################################################

named_summary_combined_df = clean %>%
  rowwise() %>%
  mutate(new_col_name = paste(Hour, Rep, sep="_")) %>%
  mutate(eventName = paste(Transcript_id, Chr, Strand, Start, End, sep=":"))
named_summary_combined_df |> dplyr::ungroup() |> dplyr::select(new_col_name, eventName) |> head()

ret_pct_matrix = named_summary_combined_df %>%
  dplyr::filter(Assignment=="Ret") %>%
  tidyr::pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              #values_fill = 0
              ) %>%
  unnest(everything()) ## %>%
  ## column_to_rownames("eventName")

ret_count_matrix = named_summary_combined_df %>%
  dplyr::filter(Assignment=="Ret") %>%
  tidyr::pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) ## %>%
  ## column_to_rownames("eventName")

spl_pct_matrix = named_summary_combined_df %>%
  dplyr::filter(Assignment=="Spl") %>%
  tidyr::pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Pct_converted,
              # values_fill = 0
              ) %>%
  unnest(everything()) ## %>%
  ## column_to_rownames("eventName")

spl_count_matrix = named_summary_combined_df %>%
  dplyr::filter(Assignment=="Spl") %>%
  tidyr::pivot_wider(id_cols=eventName,
              names_from=new_col_name,
              values_from=Total_counts,
              values_fill = 0) %>%
  unnest(everything()) ## %>%
  ## column_to_rownames("eventName")

rm(named_summary_combined_df)

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

  average0_spliced = actual_data %>% dplyr::filter(time==0, isoform=="Spliced")
  average0_spliced_weighted_mean = weighted.mean(average0_spliced$value, average0_spliced$total_counts)
  average0_retained = actual_data %>% dplyr::filter(time==0, isoform=="Retained")
  average0_retained_weighted_mean = weighted.mean(average0_retained$value, average0_retained$total_counts)

  actual_data_forced_1 = actual_data %>%
    mutate(value0 = ifelse(isoform=="Spliced", average0_spliced_weighted_mean, average0_retained_weighted_mean),
           value = value/value0)

  # we should also normalize the weighting factor so that it reflects the within-isoform differences
  # this will be important for examples where there are loads more spliced reads than retained, or vv

  total_counts_spliced = actual_data_forced_1 %>% dplyr::filter(isoform=="Spliced") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
  total_counts_retained = actual_data_forced_1 %>% dplyr::filter(isoform=="Retained") %>% dplyr::select(total_counts) %>% unlist() %>% sum()
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
in_both = intersect(rownames(ret_pct_matrix$eventName), rownames(spl_pct_matrix$eventName))

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
pvals = na.omit(pvals) %>% dplyr::filter(hl_spl<24, hl_ret<24, hl_spl>0, hl_ret>0)

# Adjust p-values for multiple hypothesis testing
pvals$padj = p.adjust(pvals$pval, method="BH")

fwrite(pvals, outfile, sep="\t")
