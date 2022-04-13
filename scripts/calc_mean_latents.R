library(data.table)

latentObs = do.call(rbind, lapply(snakemake@input[-c(1)], fread))
latentMeans = colMeans(latentObs)

obsData = read.csv(snakemake@input[[1]], sep='\t', row.names=1)

mat = matrix(latentMeans, ncol=ncol(obsData))
df = data.frame(mat)
colnames(df) = colnames(obsData)
rownames(df) = rownames(obsData)

write.csv(df, snakemake@output[["df"]])
