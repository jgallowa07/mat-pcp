"""Top-level ``snakemake`` file that runs pipeline."""


import os
configfile: "config.yaml"

# TODO ?
# if "--use_gpu" in config["chronumental_kwargs"]:
#     envvars:
#         "CUDA_VISIBLE_DEVICES"

results_files = [
    "ref/coding_sites.csv",
    "mat/filtered_metadata.tsv.gz",
    "mat/filtered_tree.pb"
]

rule all:
    """Target rule with desired output files."""
    input:
        expand(
            [os.path.join("results_{mat}", f) for f in results_files],
            mat=config["mat_trees"],
        )

# TODO make this compatible with either a local or remote mat tree
#        """
#        curl {params.mat_url} > {output.mat}
#        curl {params.meta_url} > {output.meta}
#        """
rule get_mat_tree:
    """Get the pre-built mutation-annotated tree, nwk, and metadata"""
    params:
        mat=lambda wc: config["mat_trees"][wc.mat]["mat"],
        meta=lambda wc: config["mat_trees"][wc.mat]["meta"],
    output:
        mat="results_{mat}/mat/mat_tree.pb.gz",
        meta="results_{mat}/mat/tree_metadata.tsv.gz"
    shell:
        """
        ln -sr {params.mat} {output.mat}
        ln -sr {params.meta} {output.meta}
        """

rule filter_metadata:
    """filter the mat tree to remove invalid sequences"""
    input:
        meta="results_{mat}/mat/tree_metadata.tsv.gz"
    output:
        # invalid_strains="results_{mat}/mat/invalid_strains.txt",
        filtered_metadata="results_{mat}/mat/filtered_metadata.tsv.gz"
    params:
        qc_filter_options = config["qc_filter_options"]
    log:
        notebook="results_{mat}/log/filter_metadata.ipynb",
        invalid_strains="results_{mat}/log/invalid_strains.txt",
    notebook:
        "notebooks/filter_metadata.ipynb"
    
rule filter_mat_tree:
    """
    extract valid sequences from the mat tree,
    and convert that subsetted tree to newick for chronumental.
    see https://github.com/matsengrp/s2trajectory/issues/9
    for why we double extract.
    """    
    input:
        mat_tree = rules.get_mat_tree.output.mat,
        filtered_metadata = rules.filter_metadata.output.filtered_metadata
    output:
        filtered_nwk="results_{mat}/mat/filtered_tree.nwk",
        filtered_pb="results_{mat}/mat/filtered_tree.pb"
    log:
        stderr="results_{mat}/log/filter_mat_tree.stderr",
        stdout="results_{mat}/log/filter_mat_tree.stdout",
        valid_strains="results_{mat}/log/valid_strains.txt",
    conda:
        "envs/usher.yml"
    params:
        mat_subset_kwargs=config["mat_subset_kwargs"]
    shell:
        """
        zcat {input.filtered_metadata} | cut -f1 | sed 1,1d > {log.valid_strains}
        (
            matUtils extract \
                -i {input.mat_tree} \
                -s {log.valid_strains} \
                -o {output.filtered_pb}

            matUtils extract \
                -i {output.filtered_pb} \
                {params.mat_subset_kwargs} \
                -o {output.filtered_pb}

            matUtils extract \
                -i {output.filtered_pb} \
                -t {output.filtered_nwk}
        ) > {log.stdout} 2> {log.stderr}
        """


rule get_ref_fasta:
    """Get the reference FASTA."""
    params:
        url=config["ref_fasta"],
    output:
        ref_fasta="results_{mat}/ref/ref.fa",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_fasta}"


rule get_ref_gtf:
    """Get the reference GTF."""
    params:
        url=config["ref_gtf"],
    output:
        ref_gtf="results_{mat}/ref/original_ref.gtf",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_gtf}"


rule edit_ref_gtf:
    """Edit the reference GTF with manual additions."""
    input:
        gtf=rules.get_ref_gtf.output.ref_gtf
    output:
        gtf="results_{mat}/ref/edited_ref.gtf",
    params:
        edits=config["add_to_ref_gtf"],
    script:
        "scripts/edit_ref_gtf.py"


rule ref_coding_sites:
    """Get all sites in reference that are part of a coding sequence."""
    input:
        gtf=rules.edit_ref_gtf.output.gtf,
        fasta=rules.get_ref_fasta.output.ref_fasta,
    output:
        csv="results_{mat}/ref/coding_sites.csv",
    script:
        "scripts/ref_coding_sites.py"
