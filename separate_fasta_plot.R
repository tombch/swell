library(tidyverse)

# Generating graphs similar to/the same as ones made by Nick Loman

alpha_val = 0.4

args <- commandArgs(trailingOnly = TRUE)
swell <- read.table(file("stdin"), head=TRUE, sep='\t', na.strings=c("NaN", "-"))
metadata <- read.table(args[1], head=TRUE, sep='\t', fill=TRUE)
df <- merge(x=swell, y=metadata, by.x="header", by.y="fasta_header")
# Calculate ISO week for each resulting entry
df$sequencing_submission_week = as.integer(strftime(df$sequencing_submission_date, format="%V"))

# Percentage acgt scatterplot, sequencing org on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_point(alpha = alpha_val) + 
    theme_bw()
# Percentage ambiguous scatterplot, sequencing org on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_point(alpha = alpha_val) + 
    theme_bw()

# Percentage acgt scatterplot, faceted by sequencing org
ggplot(data = df, mapping = aes(x = pc_acgt, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank()) 
# Percentage ambiguous scatterplot, faceted by sequencing org
ggplot(data = df, mapping = aes(x = pc_ambiguous, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())

# Percentage acgt scatterplot, faceted by sequencing org and run name
ggplot(data = df, mapping = aes(x = pc_acgt, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(sequencing_org_code ~ run_name, ncol=2, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank()) 
# Percentage ambiguous scatterplot, faceted by sequencing org and run name
ggplot(data = df, mapping = aes(x = pc_ambiguous, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(sequencing_org_code ~ run_name, ncol=2, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())

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

# Percentage acgt histogram, faceted by sequencing org and run name
ggplot(data = df, mapping = aes(x = pc_acgt)) +
    geom_histogram() + 
    facet_wrap(sequencing_org_code ~ run_name, ncol=2, scale="free_y") + 
    theme_bw()
# Percentage ambiguous histogram, faceted by sequencing org and run name
ggplot(data = df, mapping = aes(x = pc_ambiguous)) +
    geom_histogram() +
    facet_wrap(sequencing_org_code ~ run_name, ncol=2, scale="free_y") + 
    theme_bw()

# Average pc_acgt per ISO week, coloured by sequencing org
ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_acgt, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "avg")
# Average pc_ambiguous per ISO week, coloured by sequencing org
ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "avg")