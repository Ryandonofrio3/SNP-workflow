import pandas as pd
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from Bio.Blast import NCBIXML
import primer3


def process_data(assoc_path, probes_path, input_list):
    df_filtered = pd.read_csv(assoc_path, sep="\s+")

    # keep only columns SNP, BP, A1, A2
    df_filtered = df_filtered.loc[:, ["SNP", "BP", "A1", "A2"]]

    probes_df = pd.read_excel(probes_path)
    filtered_probes_df = probes_df[probes_df["Probe Set ID"].isin(df_filtered["SNP"])]
    filtered_probes_df = filtered_probes_df.loc[:, ["Probe Set ID", "Flank"]]

    # Merge DataFrames based on Probe Set ID / SNP
    merged_df = pd.merge(
        df_filtered, filtered_probes_df, left_on="SNP", right_on="Probe Set ID"
    )

    # Create a new column 'Flank with LVP', replace brackets with corresponding A1
    merged_df["Flank with LVP"] = [
        re.sub(r"\[.*?\]", a1, flank)
        for a1, flank in zip(merged_df["A1"], merged_df["Flank"])
    ]

    merged_df["Flank with Orlando"] = [
        re.sub(r"\[.*?\]", a2, flank)
        for a2, flank in zip(merged_df["A2"], merged_df["Flank"])
    ]

    # sort the data frame by BP smallest to largest
    merged_df.sort_values(by=["BP"], inplace=True)

    merged_df = merged_df[merged_df["SNP"].isin(input_list)]

    return merged_df
