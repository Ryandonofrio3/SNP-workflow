import pandas as pd
import logging
import re
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

logging.basicConfig(level=logging.INFO)

assoc_path = "assoc_res2.assoc"

# Read the assoc file with pandas
# sep="\s+" tells pandas to use one or more spaces as the separator
df = pd.read_csv(assoc_path, sep="\s+")

logging.info("Successfully read assoc file.")

header_names = df.columns.tolist()
print(header_names)

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


print(merged_df)

# Always tell NCBI who you are
Entrez.email = "ryandonofrio@gmail.com"

# Create an empty DataFrame to store final results
final_df = pd.DataFrame()


def blast_sequence(sequence):
    result_handle = NCBIWWW.qblast(
        "blastn", "nt", sequence, entrez_query="txid7159[Organism:exp]"
    )
    blast_records = NCBIXML.read(result_handle)
    result_handle.close()
    return blast_records


# For each row in the dataframe, call blast_sequence
for index, row in merged_df.iterrows():
    blast_records = blast_sequence(row["Flank"])
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.query == row["Flank with LVP"]:
                final_df = final_df.append(row)

# Save the final DataFrame
final_df.to_csv("final_df.csv", index=False)
