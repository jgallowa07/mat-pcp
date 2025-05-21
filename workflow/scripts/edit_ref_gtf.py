input_gtf = snakemake.input.gtf
output_gtf = snakemake.output.gtf
edits = snakemake.params.edits

import csv

import pandas as pd

gtf = pd.read_csv(
    input_gtf,
    sep="\t",
    header=None,
    names=[
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ],
)

seqname = gtf["seqname"].unique()
assert len(seqname) == 1, seqname
seqname = seqname[0]

source = gtf["source"].unique()
assert len(source) == 1, source
source = source[0]

assert (
    gtf.query("feature == 'transcript'")["start"].tolist()
    == sorted(gtf.query("feature == 'transcript'")["start"])
), "transcript starts not sorted"

for gene, [start, end] in edits.items():
    print(f"Adding {gene=}, {start=}, {end=}")
    # get index of first entry in GTF that starts **after** this feature
    after_gtf = gtf.query("feature == 'transcript'").query("start > @start")
    if len(after_gtf):
        i = after_gtf.index[0]
    else:
        i = len(gtf)
    records = []
    for feature in ["transcript", "exon", "CDS", "start_codon", "stop_codon"]:
        records.append(
            (
                seqname,
                source,
                feature,
                start if feature != "stop_codon" else end - 2,
                start + 2 if feature == "start" else end - 3 if feature == "CDS" else end,
                ".",
                "+",
                "." if feature in ["transcript", "exon"] else 0,
                (
                    f'gene_id "{gene}"; transcript_id "{gene};"'
                    if feature == "transcript"
                    else
                    f'gene_id "{gene}"; transcript_id "{gene}"; exon_number "1"; exon_id "{gene}";'
                )
            )
        )
    add_gtf = pd.DataFrame.from_records(records, columns=gtf.columns)
    
    gtf = pd.concat([gtf.loc[: i - 1], add_gtf, gtf.loc[i: ]], ignore_index=True)

gtf.to_csv(output_gtf, index=False, header=None, sep="\t", quoting=csv.QUOTE_NONE)
