library(tidyverse)

swell <- read.table(file("stdin"), head=TRUE, sep='\t', na.strings=c("NaN", "-"))
metadata <- read.table("majora.metadata.matched.tsv", head=TRUE, sep='\t', fill=TRUE)

# Filter metadata by the past x days
filter_date = Sys.Date() - 300
past_week_metadata <- filter(metadata, collection_date >= filter_date)
# Merge ignores rows that don't have a matching fasta_header
df <- merge(x=swell, y=past_week_metadata, by.x="fasta_header", by.y="fasta_header")

# Percentage acgt scatterplot, organisation on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_point() + 
    theme_bw()

# Percentage ambiguous scatterplot, organisation on x axis
ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_point() + 
    theme_bw()

# Percentage acgt histogram, faceted by organisation
ggplot(data = df, mapping = aes(x = pc_acgt)) +
    geom_histogram() + 
    facet_wrap(. ~ sequencing_org_code) + 
    theme_bw()

# Percentage ambiguous histogram, faceted by organisation
ggplot(data = df, mapping = aes(x = pc_ambiguous)) +
    geom_histogram() +
    facet_wrap(. ~ sequencing_org_code) + 
    theme_bw()