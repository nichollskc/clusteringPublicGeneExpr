Copied data/pathways from ~/rds/rds-cew54-basis/People/KATH/clustering_geneExpr/data/pathways/

Saved this from affymetrix https://www.affymetrix.com/products_services/arrays/specific/hugene-1_1-st-v1_strip.affx#1_4
HuGene-1_1-st-v1.na36.hg19.probeset.csv.zip

snakemake  -k -j 1000 --cluster-config slurm_config.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -c {cluster.cpus-per-task}   -t {cluster.time} --output {cluster.error} -J {cluster.job}" all_signatures
