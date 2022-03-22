wildcard_constraints:
    dataset = "\w+",
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
        "get_biomart.R"

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
        "standardise_smith.R"

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
        "standardise_chaussabel.R"

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
        "standardise_coulson.R"

###############################################################################
# Normalisation                                                               #
###############################################################################

rule apply_vsn:
    input:
        "data/datasets/{dataset}/expression.tsv",
    output:
        "data/datasets/{dataset}/expression_vsn.tsv",
    script:
        "apply_vsn.R"

rule apply_mad:
    input:
        "data/datasets/{dataset}/expression{normalisation}.tsv",
    output:
        "data/datasets/{dataset}/expression{normalisation}_mad.tsv",
    script:
        "apply_mad_normalisation.R"

rule calc_logfc:
    input:
        expression="data/datasets/{dataset}/expression{normalisation}.tsv",
        sample_info="data/datasets/{dataset}/sample_info.tsv",
        probe_to_gene="data/datasets/{dataset}/probe_to_gene.tsv",
    output:
        gene_means="data/datasets/{dataset}/gene_means{normalisation}-{genelist_name}.tsv",
        logfc="data/datasets/{dataset}/logfc{normalisation}-{genelist_name}.tsv",
        average_logfc="data/datasets/{dataset}/average_logfc{normalisation}-{genelist_name}.tsv",
    script:
        "calculate_logfc_signature.R"

rule all_signatures:
    input:
        expand("data/datasets/{dataset}/average_logfc{normalisation}-{genelist_name}.tsv",
               dataset=["smith", "chaussabelA", "chaussabelB", "coulson"],
               normalisation=["", "_vsn", "_vsn_mad", "_mad"],
               genelist_name=["ifn", "nk_dipp", "exhaustion_down_wherry",
                   "cd4_activation_green_30", "cd4_activation_yellow_30",
                   "cd4_activation_black_30"])
