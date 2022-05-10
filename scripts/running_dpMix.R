library(DPMUnc)

read_matrix <- function(filename) {
    table <- read.table(filename, sep='\t', header=TRUE, row.names=1)
    mat <- as.matrix(table)
    colnames(mat) <- NULL
    rownames(mat) <- NULL
    return (mat)
}

# Snakemake variables
obsfile <- snakemake@input[['obs']]
varfile <- snakemake@input[['var']]
seed  <- snakemake@wildcards[['seed']]
iterations <- snakemake@params[['iterations']]
thinningFreq <- snakemake@params[['thinningFreq']]

obsData = read_matrix(obsfile)
obsVars = read_matrix(varfile)

alpha0 = 2; kappa0 = 0.5; beta0 = 0.2 * mean(apply(obsData, 2, var))

directory = dirname(snakemake@output[[1]])
DPMUnc(obsData, obsVars, saveFileDir=directory, seed=seed,
       kappa0=kappa0, alpha0=alpha0, beta0=beta0,
       nIts=iterations, thinningFreq=thinningFreq, saveLatentObs=TRUE, scaleData=FALSE)
