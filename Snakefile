localrules: find_signature_intersection

wildcard_constraints:
    dataset = "smith|chaussabelA|chaussabelB|coulson",
    normalisation = "\w*",

###############################################################################
# Smith data                                                                  #
###############################################################################

rule fetch_smith:
    output:
        "data/datasets/smith/processed.zip",
        "data/datasets/smith/PBMC_norm.txt",
    shell:
        "wget -O {output[0]} https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-145/E-MTAB-145.processed.1.zip && unzip {output[0]} -d data/datasets/smith"

rule fetch_smith_sample_info:
    output:
        "data/datasets/smith/sample_info.txt",
    shell:
        "wget -O {output} https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-145/E-MTAB-145.sdrf.txt"

rule fetch_smith_array:
    output:
        "data/datasets/smith/A-MEXP-1152.adf.txt",
    shell:
        "wget -O {output} https://www.ebi.ac.uk/arrayexpress/files/A-MEXP-1152/A-MEXP-1152.adf.txt"

rule download_biomart:
    output:
        "data/datasets/biomart.Rdata",
    script:
        "scripts/get_biomart.R"

rule standardise_smith:
    input:
        "data/datasets/smith/A-MEXP-1152.adf.txt",
        "data/datasets/smith/sample_info.txt",
        "data/datasets/smith/PBMC_norm.txt",
    output:
        "data/datasets/smith/sample_info.tsv",
        "data/datasets/smith/probe_to_gene.tsv",
        "data/datasets/smith/expression.tsv",
    script:
        "scripts/standardise_smith.R"

###############################################################################
# Chaussabel data                                                             #
###############################################################################

rule fetch_chaussabel:
    output:
        "data/datasets/chaussabel/processed.zip",
        "data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830054.txt",
        "data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830055.txt",
    shell:
        "wget -O {output[0]} https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-11907/E-GEOD-11907.processed.1.zip && unzip {output[0]} -d data/datasets/chaussabel"

rule fetch_chaussabel_sample_info:
    output:
        "data/datasets/chaussabel/sample_info.txt",
    shell:
        "wget -O {output} https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-11907/E-GEOD-11907.sdrf.txt"

#http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus
rule fetch_affymetrix_translation:
    output:
        "data/datasets/HG-U133_Plus_2.na36.annot.csv",
    shell:
        "wget -O tmp.zip http://www.affymetrix.com/Auth/analysis/downloads/na36/ivt/HG-U133_Plus_2.na36.annot.csv.zip && unzip tmp.zip -d data/datasets"

rule standardise_chaussabel:
    input:
        "data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830054.txt",
        "data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830055.txt",
        "data/datasets/chaussabel/sample_info.txt",
        "data/datasets/chaussabel/HG-U133_Plus_2.na36.annot.csv",
    output:
        "data/datasets/chaussabelA/sample_info.tsv",
        "data/datasets/chaussabelA/probe_to_gene.tsv",
        "data/datasets/chaussabelA/expression.tsv",
        "data/datasets/chaussabelB/sample_info.tsv",
        "data/datasets/chaussabelB/probe_to_gene.tsv",
        "data/datasets/chaussabelB/expression.tsv",
    script:
        "scripts/standardise_chaussabel.R"

###############################################################################
# Coulson data                                                                #
###############################################################################

rule fetch_coulson:
    output:
        "data/datasets/coulson/ifn_signature.data_matrix.tab",
    shell:
        "wget -O tmp.zip https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1724/E-MTAB-1724.processed.1.zip && unzip tmp.zip -d data/datasets/coulson"

rule fetch_coulson_sample_info:
    output:
        "data/datasets/coulson/sample_info.txt",
    shell:
        "wget -O {output} https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1724/E-MTAB-1724.sdrf.txt"

rule fetch_coulson_array_comments:
    output:
        "data/datasets/coulson/A-AFFY-168_comments.txt",
    shell:
        "wget -O {output} https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-168/A-AFFY-168_comments.txt"

rule fetch_coulson_array:
    output:
        "data/datasets/coulson/A-AFFY-168.adf.txt",
    shell:
        "wget -O {output} https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-168/A-AFFY-168.adf.txt"

rule standardise_coulson:
    input:
        "data/datasets/coulson/HuGene-1_1-st-v1.na36.hg19.probeset.csv",
        "data/datasets/coulson/sample_info.txt",
        "data/datasets/coulson/ifn_signature.data_matrix.tab",
    output:
        "data/datasets/coulson/sample_info.tsv",
        "data/datasets/coulson/probe_to_gene.tsv",
        "data/datasets/coulson/expression.tsv",
    script:
        "scripts/standardise_coulson.R"

###############################################################################
# Normalisation                                                               #
###############################################################################

rule apply_vsn:
    input:
        "data/datasets/{dataset}/expression.tsv",
    output:
        "data/datasets/{dataset}/expression_vsn.tsv",
    script:
        "scripts/apply_vsn.R"

rule apply_mad:
    input:
        "data/datasets/{dataset}/expression{normalisation}.tsv",
    output:
        "data/datasets/{dataset}/expression{normalisation}_mad.tsv",
    script:
        "scripts/apply_mad_normalisation.R"

rule calc_logfc:
    input:
        genelist="data/pathways/processed/{genelist_name}.csv",
        expression="data/datasets/{dataset}/expression{normalisation}.tsv",
        sample_info="data/datasets/{dataset}/sample_info.tsv",
        probe_to_gene="data/datasets/{dataset}/probe_to_gene.tsv",
    output:
        gene_means="data/datasets/{dataset}/gene_means{normalisation}-{genelist_name}.tsv",
        logfc="data/datasets/{dataset}/logfc{normalisation}-{genelist_name}.tsv",
        average_logfc="data/datasets/{dataset}/average_logfc{normalisation}-{genelist_name}.tsv",
        median_logfc="data/datasets/{dataset}/median_logfc{normalisation}-{genelist_name}.tsv",
    script:
        "scripts/calculate_logfc_signature.R"

rule rank_all_genes:
    input:
        expression="data/datasets/{dataset}/expression{normalisation}.tsv",
        probe_to_gene="data/datasets/{dataset}/probe_to_gene.tsv",
    output:
        ranks="data/datasets/{dataset}/ranked_expression{normalisation}.tsv",
        gene_means="data/datasets/{dataset}/gene_means{normalisation}.tsv",
    script:
        "scripts/rank_genes.R"

rule calc_mean_ranks:
    input:
        genelist="data/pathways/processed/{genelist_name}.csv",
        sample_info="data/datasets/{dataset}/sample_info.tsv",
        ranks="data/datasets/{dataset}/ranked_expression{normalisation}.tsv",
    output:
        ranks="data/datasets/{dataset}/ranks{normalisation}-{genelist_name}.tsv",
        mean_ranks="data/datasets/{dataset}/mean_ranks{normalisation}-{genelist_name}.tsv",
    script:
        "scripts/calc_mean_ranks.R"

rule all_signatures:
    input:
        expand("data/datasets/{dataset}/{summary}{normalisation}-{genelist_name}.tsv",
               summary=["average_logfc", "mean_ranks"],
               dataset=["smith", "chaussabelA", "chaussabelB", "coulson"],
               normalisation=["", "_vsn", "_vsn_mad", "_mad"],
               genelist_name=["ifn", "nk_dipp", "exhaustion_down_wherry",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),

rule compare_datasets:
    input:
        expand("data/datasets/{dataset}/{{summary}}{{normalisation}}-{{genelist_name}}.tsv",
               dataset=["smith", "chaussabelA", "chaussabelB", "coulson"]),
    wildcard_constraints:
        summary="average_logfc|mean_ranks"
    output:
        "plots/compare_datasets_{summary}--{normalisation}-{genelist_name}.png",
        "plots/compare_datasets_{summary}--se_{normalisation}-{genelist_name}.png",
    script:
        "scripts/compare_datasets.R"

rule all_comparisons:
    input:
        expand("plots/compare_datasets_{summary}--{normalisation}-{genelist_name}.png",
               summary=["average_logfc", "mean_ranks"],
               normalisation=["", "_vsn", "_vsn_mad", "_mad"],
               genelist_name=["ifn", "nk_dipp", "exhaustion_down_wherry",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),

rule find_signature_intersection:
    input:
        expand("data/datasets/{dataset}/logfc-{{genelist_name}}.tsv",
               dataset=["smith", "chaussabelA", "chaussabelB", "coulson"]),
    output:
        "data/pathways/processed/intersection_{genelist_name}.csv",
    script:
        "scripts/find_signature_intersection.R"

rule all_comparisons_intersection:
    input:
        expand("plots/compare_datasets_{summary}--{normalisation}-intersection_{genelist_name}.png",
               summary=["average_logfc", "mean_ranks"],
               normalisation=["", "_vsn", "_vsn_mad", "_mad"],
               genelist_name=["ifn", "nk_dipp", "exhaustion_down_wherry",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),

rule combine_datasets:
    input:
        expand("data/datasets/{dataset}/average_logfc{{normalisation}}-{genelist_name}.tsv",
               dataset=["smith", "chaussabelA", "chaussabelB", "coulson"],
               genelist_name=["ifn", "nk_dipp", "exhaustion_down_wherry",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),
    output:
        "data/datasets/combined/average_logfc{normalisation}.tsv",
    script:
        "scripts/combine_datasets.R"

rule make_individual_dataset:
    input:
        expand("data/datasets/{{dataset}}/average_logfc{{normalisation}}-{genelist_name}.tsv",
               genelist_name=["exhaustion_down_wherry",
                   "ifn", "nk_dipp",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),
    output:
        expand("data/processed/{{dataset}}{{normalisation}}_logfc{var}/{signature}/{outputfile}",
               var = ["", "_novar"],
               outputfile = ["obsVars.tsv", "obsData.tsv"],
               signature = ["signatures_v2", "signatures_v3",
                   "ifn", "nk_dipp", "exhaustion_down_wherry",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),
    script:
        "scripts/make_individual_dataset.R"

rule make_individual_dataset_fix:
    input:
        expand("data/datasets/{{dataset}}/average_logfc{{normalisation}}-{genelist_name}.tsv",
               genelist_name=["exhaustion_down_wherry",
                   "ifn", "nk_dipp",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"]),
    output:
        expand("data/processed/{{dataset}}{{normalisation}}_logfc{var}/{signature}/{outputfile}",
               var = ["", "_novar"],
               outputfile = ["obsVars.tsv", "obsData.tsv"],
               signature = ["signatures_v4"]),
    script:
        "scripts/make_individual_dataset.R"

rule all_datasets:
    input:
        expand("data/processed/{dataset}{processing}/signatures_v4/obsData.tsv",
               processing = ["_vsn_mad_logfc", "_vsn_mad_logfc_novar"],
               dataset=["smith", "chaussabelA", "chaussabelB", "coulson"]),
