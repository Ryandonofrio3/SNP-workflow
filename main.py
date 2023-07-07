import pandas as pd
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from Bio.Blast import NCBIXML
import primer3

# Configure logging
logging.basicConfig(
    filename="blast.log",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s:%(message)s",
)

assoc_path = "assoc_res2.assoc"

df = pd.read_csv(assoc_path, sep="\s+")

logging.info("Successfully read assoc file.")

header_names = df.columns.tolist()

# Filter DataFrame
df_filtered = df[
    ((df["F_A"] == 1.0) & (df["F_U"] == 0.0))
    | ((df["F_A"] == 0.0) & (df["F_U"] == 1.0))
]

# keep only columns SNP, BP, A1, A2
df_filtered = df_filtered.loc[:, ["SNP", "BP", "A1", "A2"]]

probes_df = pd.read_excel("probes.xlsx")
filtered_probes_df = probes_df[probes_df["Probe Set ID"].isin(df_filtered["SNP"])]
filtered_probes_df = filtered_probes_df.loc[:, ["Probe Set ID", "Flank"]]

# Merge DataFrames based on Probe Set ID / SNP
merged_df = pd.merge(
    df_filtered, filtered_probes_df, left_on="SNP", right_on="Probe Set ID"
)


# Create a new column 'Flank with LVP', replace brackets with corresponding A1 values
merged_df["Flank with LVP"] = [
    re.sub(r"\[.*?\]", a1, flank)
    for a1, flank in zip(merged_df["A1"], merged_df["Flank"])
]


# sort the data frame by BP smallest to largest
merged_df.sort_values(by=["BP"], inplace=True)

# Write the 'Flank with LVP' sequences to a FASTA file
with open("input.fasta", "w") as f:
    for i, row in merged_df.iterrows():
        sequence = SeqRecord(Seq(row["Flank with LVP"]), id=row["SNP"])
        SeqIO.write(sequence, f, "fasta")


subprocess.run(
    [
        "blastn",
        "-query",
        "input.fasta",
        "-db",
        "AAE5.fna",
        "-out",
        "output.xml",
        "-outfmt",
        "5",
    ]
)


# Open a new dataframe to hold the final results
final_df = pd.DataFrame()

with open("output.xml") as result_handle:
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        snp = record.query.split()[0]
        row = merged_df.loc[merged_df["SNP"] == snp].iloc[0]

        for alignment in record.alignments:
            for hsp in alignment.hsps:
                query = hsp.query
                match = hsp.match
                sbjct = hsp.sbjct

                if query == sbjct:
                    final_df = final_df.append(row)

print(final_df.head())

primer_params = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 55.0,
    "PRIMER_MIN_TM": 48.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 20.0,
    "PRIMER_MAX_GC": 80.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_PRODUCT_SIZE_RANGE": [[200, 500]],
}


def design_primers(sequence, seq_id):
    primer_params["SEQUENCE_ID"] = seq_id
    primer_params


primer_results = {}


# Load the whole genome
genome_record = SeqIO.read("AAE5.fna", "fasta")

# Assume you have a DataFrame with SNPs and their positions
for idx, row in final_df.iterrows():
    sequence = row["Flank with LVP"].replace("/", "")  # remove slash if present
    seq_id = row["SNP"]
    bp = row["BP"]

    # Define the region around the SNP (+-250 bp)
    start = max(0, bp - 250)
    end = bp + 250 + len(sequence)  # adjust for the length of the sequence

    # Extract the region from the genome
    region = str(genome_record.seq[start:end])

    # Design primers for this sequence
    results = design_primers(region, seq_id)

    # Store the results in the dictionary
    primer_results[seq_id] = results

# Convert the results to a DataFrame for easier viewing
primer_results_df = pd.DataFrame(primer_results).T
print(primer_results_df.head())
