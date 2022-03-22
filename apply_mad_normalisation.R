library(dplyr)
library(tidyr)

normalise_sample_to_constant_mad <- function(sample_values) {
    # Calculate the MAD
    sample_mad = mad(sample_values)

    normalised = (sample_values - median(sample_values)) / sample_mad

    return (normalised)
}

filename = snakemake@input[[1]]
output_file = snakemake@output[[1]]
expression = read.csv(filename, row.names=1, check.names=FALSE, sep='\t')

mad_expression = expression %>% mutate(across(everything(), normalise_sample_to_constant_mad))

mad_expression = cbind(probe_id=rownames(mad_expression),
                       mad_expression)
write.table(mad_expression, output_file, row.names=FALSE, col.names=TRUE, sep='\t')
