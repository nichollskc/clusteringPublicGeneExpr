library(data.table)
library(ggplot2)

##### Coulson

data=fread("grep -v 'log2 intensity' data/datasets/coulson/ifn_signature.data_matrix.tab", nrows=1000)
dat=as.matrix(data)[,-1]
rownames(dat)=data[[1]]

pca=prcomp(t(dat))
str(pca)

batch=as.numeric(as.factor(sub("_.*","",rownames(pca$x)))) # batch id, could be dropped if not known
disease=as.factor(sub("PBMC.*", "PBMC", sub("_.*", "", sub("^\\d_","",rownames(pca$x)))))
table(batch)
table(disease)
# Each batch has different mix of diseases. Only batch 1 has SLE
table(batch, disease)

df = data.frame(pca$x)
df$disease = disease
df$batch = as.factor(batch)

# Batch is very dominant in PCA
ggplot(df, aes(x=PC1, y=PC2, colour=batch)) + geom_point()
ggsave("plots/pca_coulson.png", width=6, height=6)

# Grouping is now clearly disease based
batch1_pca = prcomp(t(dat[, batch == 1]))
df = data.frame(batch1_pca$x)
df$disease = disease[batch == 1]
ggplot(df, aes(x=PC1, y=PC2, colour=disease)) + geom_point()
ggsave("plots/pca_coulson1.png", width=6, height=6)

##### Smith

data=read.table("data/datasets/smith/expression.tsv", nrows=1000, sep="\t", row.names=1, header=TRUE)

pca=prcomp(t(data))
str(pca)

disease=as.factor(sub("_.*","",rownames(pca$x)))
table(disease)

df = data.frame(pca$x)
df$disease = disease

# PCA not dominated by disease
ggplot(df, aes(x=PC1, y=PC2, colour=disease)) + geom_point()
ggsave("plots/pca_smith.png", width=6, height=6)

###### Chaussabel A

data=read.table("data/datasets/chaussabelA/expression.tsv", nrows=1000, sep="\t", row.names=1, header=TRUE)

# Load Chaussabel sample_info with ftp info
source("scripts/standardise_chaussabel.R")

pca=prcomp(t(data))
str(pca)

disease=as.factor(sub("_.*","",rownames(pca$x)))
table(disease)

df = data.frame(pca$x)
df$disease = disease
df$neat_sample_id = rownames(df)
df = merge(df, sample_info, by="neat_sample_id")

# PCA not dominated by disease
ggplot(df, aes(x=PC1, y=PC2, colour=disease)) + geom_point()
ggsave("plots/pca_chaussabelA.png", width=6, height=6)

# PCA not dominated by ftp file either
ggplot(df, aes(x=PC1, y=PC2, colour=ftp)) + geom_point()
ggsave("plots/pca_chaussabelA_ftp.png", width=6, height=6)

###### Chaussabel B

data=read.table("data/datasets/chaussabelB/expression.tsv", nrows=1000, sep="\t", row.names=1, header=TRUE)

pca=prcomp(t(data))
str(pca)

disease=as.factor(sub("_.*","",rownames(pca$x)))
table(disease)

df = data.frame(pca$x)
df$disease = disease
df$neat_sample_id = rownames(df)
df = merge(df, sample_info, by="neat_sample_id")

# PCA not dominated by disease
ggplot(df, aes(x=PC1, y=PC2, colour=disease)) + geom_point()
ggsave("plots/pca_chaussabelB.png", width=6, height=6)

# PCA not dominated by ftp file either
ggplot(df, aes(x=PC1, y=PC2, colour=ftp)) + geom_point()
ggsave("plots/pca_chaussabelB_ftp.png", width=6, height=6)

