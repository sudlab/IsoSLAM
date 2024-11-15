library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

input_folder =  args[1]

input_regex = args[2]

outfile = args[3]

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

fwrite(combined_df, outfile, sep="\t")