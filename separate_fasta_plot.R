library(tidyverse)

swell <- read.table(file("stdin"), head=TRUE, sep='\t', na.strings=c("NaN", "-"))
metadata <- read.table("majora.metadata.matched.tsv", head=TRUE, sep='\t', fill=TRUE)
# Filter metadata by the past x days
filter_date = Sys.Date() - 7
past_week_metadata <- filter(metadata, collection_date >= filter_date)
# Merge ignores rows that don't have a matching fasta_header
df <- merge(x=swell, y=past_week_metadata, by.x="header", by.y="fasta_header")

# Percentage acgt scatterplot, sequencing org on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt, colour = instrument_make)) +
    geom_point() + 
    theme_bw()

# Percentage ambiguous scatterplot, sequencing org on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous, colour = instrument_make)) +
    geom_point() + 
    theme_bw()

# Percentage acgt histogram, faceted by sequencing org
ggplot(data = df, mapping = aes(x = pc_acgt, fill = instrument_make)) +
    geom_histogram() + 
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()

# Percentage ambiguous histogram, faceted by sequencing org
ggplot(data = df, mapping = aes(x = pc_ambiguous, fill = instrument_make)) +
    geom_histogram() +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()

# colour or facet by
# run_group or run_name