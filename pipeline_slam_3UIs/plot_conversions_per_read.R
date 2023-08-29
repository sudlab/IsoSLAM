library(dplyr)
library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
file = fread(args[1])

summarized = file %>% group_by(V1) %>% summarize(n=n())
five_plus = summarized %>% filter(V1>=5) %>% summarize(n=sum(n)) %>% mutate(V1 = "5+") %>% relocate(V1)
summarized = summarized %>% filter(V1<5) %>% mutate(V1= as.character(V1)) %>% add_row(five_plus)
summarized$V1 = factor(summarized$V1, levels=c("0", "1", "2", "3", "4", "5+"))
total = sum(summarized$n)
summarized = summarized %>% mutate(pct = (n/total)*100)

summarized %>% ggplot(aes(x=V1, y=pct)) + geom_col(color="black", fill="orange") +
  scale_y_continuous(breaks=seq(0,100,by=5)) +
  labs(x = "Number of conversions",
       y= "Percentage of reads") + 
  theme_bw(base_size=8) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=8)) -> saveme

png(args[2], width = 3.5, height=3.5, units="in", res=1200)
saveme
dev.off()