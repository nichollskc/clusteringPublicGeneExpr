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

rule fetch_chaussabel:
    output:
        "data/datasets/chaussabel/processed.zip",
#        "data/datasets/smith/PBMC_norm.txt",
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
