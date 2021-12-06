library(tidyverse)

num_days = 100
num_weeks = 5

args <- commandArgs(trailingOnly = TRUE)
swell <- read.table(file("stdin"), head=TRUE, sep='\t', na.strings=c("NaN", "-"))
metadata <- read.table(args[1], head=TRUE, sep='\t', fill=TRUE)
df <- merge(x=swell, y=metadata, by.x="header", by.y="fasta_header")
print(df)

# Filter df by the past x days
# filter_date = Sys.Date() - num_days
# df <- filter(df, sequencing_submission_date >= filter_date) 
# Calculate ISO week for each resulting entry
df$sequencing_submission_week = as.integer(strftime(df$sequencing_submission_date, format="%V"))
# Use last week as the final week in the filtering
this_week = as.integer(strftime(Sys.Date(), format="%V"))
# last_week = this_week - 1
# Filter df to be only the past x weeks up to and including last week
# df <- filter(df, sequencing_submission_week > last_week - num_weeks) 
df <- filter(df, sequencing_submission_week < this_week) 


# Percentage acgt scatterplot, sequencing org on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_point() + 
    theme_bw()

# Percentage ambiguous scatterplot, sequencing org on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_point() + 
    theme_bw()

# Percentage acgt histogram, faceted by sequencing org
ggplot(data = df, mapping = aes(x = pc_acgt)) +
    geom_histogram() + 
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()

# Percentage ambiguous histogram, faceted by sequencing org
ggplot(data = df, mapping = aes(x = pc_ambiguous)) +
    geom_histogram() +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()

# Average pc_acgt per ISO week, coloured by sequencing org
ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_acgt, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "avg")