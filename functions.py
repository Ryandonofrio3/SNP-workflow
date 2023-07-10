import pandas as pd
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from Bio.Blast import NCBIXML
import primer3
import time
import warnings
import csv


def process_data(assoc_path, probes_path):
    df = pd.read_csv(assoc_path, sep="\s+")

    # Filter DataFrame
    df_filtered = df[
        ((df["F_A"] == 1.0) & (df["F_U"] == 0.0))
        | ((df["F_A"] == 0.0) & (df["F_U"] == 1.0))
    ]

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

    # sort the data frame by BP smallest to largest
    merged_df.sort_values(by=["BP"], inplace=True)

    return merged_df


def run_blast(YourName, merged_df, AAE5_loc, today):
    input_file_name = f"{YourName}_{today}_input.fasta"
    output_file_name = f"{YourName}_{today}_output.xml"

    with open(input_file_name, "w") as f:
        for i, row in merged_df.iterrows():
            sequence = SeqRecord(Seq(row["Flank with LVP"]), id=row["SNP"])
            SeqIO.write(sequence, f, "fasta")

    subprocess.run(
        [
            "blastn",
            "-query",
            input_file_name,
            "-db",
            AAE5_loc,
            "-out",
            output_file_name,
            "-outfmt",
            "5",
        ]
    )
    return output_file_name


def process_blast_results(output_file_name, merged_df):
    final_df = pd.DataFrame()

    with open(output_file_name) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        for record in blast_records:
            snp = record.query.split()[0]
            row = merged_df.loc[merged_df["SNP"] == snp].iloc[0]

            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    query = hsp.query
                    match = hsp.match
                    sbjct = hsp.sbjct
                    hit_from = hsp.sbjct_start
                    hit_to = hsp.sbjct_end
                    hit_length = hit_to - hit_from
                    hit_def = alignment.hit_def

                    # Extract the chromosome number from the hit_def string
                    match = re.search("chromosome (\d+),", hit_def)
                    if match:
                        chromosome = match.group(1)
                    else:
                        chromosome = None

                    if query == sbjct and chromosome is not None:
                        row["Hit from"] = hit_from
                        row["Hit to"] = hit_to
                        row["Hit length"] = hit_length
                        row["Chromosome"] = chromosome
                        final_df = final_df.append(row)

    final_df = final_df[final_df["Hit length"] == 70]
    logging.info("Successfully filtered results.")
    logging.info(final_df.head(5))

    return final_df


def process_sequences_and_concatenate(final_df, AAE5_loc, chromosome_ids):
    # Write the 'Flank with LVP' sequences to a FASTA file
    with open("final_input.fasta", "w") as f:
        for i, row in final_df.iterrows():
            sequence = SeqRecord(Seq(row["Flank with LVP"]), id=row["SNP"])
            SeqIO.write(sequence, f, "fasta")

    # Create a dictionary to hold the entire genome
    genome = SeqIO.to_dict(SeqIO.parse(AAE5_loc, "fasta"))

    concatenated_genome = ""
    for chrom_id in chromosome_ids:
        concatenated_genome += str(genome[chrom_id].seq)

    final_df["Sequence"] = ""

    return final_df, concatenated_genome


def generate_final_output(final_df, concatenated_genome, YourName, today):
    final_output_file_name = f"{YourName}_{today}_final_output.fasta"

    for i, row in final_df.iterrows():
        start = row["BP"] - 1000
        end = row["BP"] + 1000
        sequence = concatenated_genome[start - 1 : end]

        lvp = row["Flank with LVP"]
        sequence = sequence[:1000] + lvp + sequence[1000:]

        final_df.at[i, "Sequence"] = sequence

    with open(final_output_file_name, "w") as f:
        for i, row in final_df.iterrows():
            sequence = SeqRecord(Seq(row["Sequence"]), id=row["SNP"])
            SeqIO.write(sequence, f, "fasta")

    return final_df


def design_primers(final_df, desired_tm=60.0, desired_product_size=500):
    # Create new columns to store the primer information
    final_df["Best Pair"] = ""
    final_df["Score"] = ""
    final_df["Left Primer"] = ""
    final_df["Right Primer"] = ""
    final_df["Left TM"] = ""
    final_df["Right TM"] = ""
    final_df["Product Size"] = ""

    for index, row in final_df.iterrows():
        result = primer3.bindings.design_primers(
            seq_args={
                "SEQUENCE_ID": row["SNP"],
                "SEQUENCE_TEMPLATE": row["Sequence"],
            },
            global_args={
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_PICK_INTERNAL_OLIGO": 1,
                "PRIMER_INTERNAL_MAX_SELF_END": 8,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 25,
                "PRIMER_OPT_TM": 60.0,
                "PRIMER_MIN_TM": 57.0,
                "PRIMER_MAX_TM": 63.0,
                "PRIMER_MIN_GC": 20.0,
                "PRIMER_MAX_GC": 80.0,
                "PRIMER_MAX_POLY_X": 100,
                "PRIMER_INTERNAL_MAX_POLY_X": 100,
                "PRIMER_SALT_MONOVALENT": 50.0,
                "PRIMER_DNA_CONC": 50.0,
                "PRIMER_MAX_NS_ACCEPTED": 0,
                "PRIMER_MAX_SELF_ANY": 12,
                "PRIMER_MAX_SELF_END": 8,
                "PRIMER_PAIR_MAX_COMPL_ANY": 12,
                "PRIMER_PAIR_MAX_COMPL_END": 8,
            },
        )

        scores = {}

        for i in range(result["PRIMER_PAIR_NUM_RETURNED"]):
            left_tm_diff = abs(result[f"PRIMER_LEFT_{i}_TM"] - desired_tm)
            right_tm_diff = abs(result[f"PRIMER_RIGHT_{i}_TM"] - desired_tm)
            product_size = abs(
                result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"] - desired_product_size
            )

            score = left_tm_diff + right_tm_diff + 3 * product_size

            scores[i] = score

        best_pair = min(scores, key=scores.get)

        final_df.at[index, "Best Pair"] = best_pair
        final_df.at[index, "Score"] = scores[best_pair]
        final_df.at[index, "Left Primer"] = result[f"PRIMER_LEFT_{best_pair}_SEQUENCE"]
        final_df.at[index, "Right Primer"] = result[
            f"PRIMER_RIGHT_{best_pair}_SEQUENCE"
        ]
        final_df.at[index, "Left TM"] = result[f"PRIMER_LEFT_{best_pair}_TM"]
        final_df.at[index, "Right TM"] = result[f"PRIMER_RIGHT_{best_pair}_TM"]
        final_df.at[index, "Product Size"] = result[
            f"PRIMER_PAIR_{best_pair}_PRODUCT_SIZE"
        ]

    return final_df


# def convert_assoc_to_csv(assoc_file_path, csv_file_path):
#     csv_data = []
#     with open(assoc_file_path, "r") as assoc_file:
#         assoc_reader = csv.reader(
#             assoc_file, delimiter="\t"
#         )  # assuming your original file is tab-separated

#         for row in assoc_reader:
#             csv_data.append(row)

#     # Split the first row into separate columns
#     headers = csv_data[0][0].split()

#     # Remove the header row from csv_data
#     csv_data = csv_data[1:]

#     with open(csv_file_path, "w", newline="") as csv_file:
#         csv_writer = csv.writer(
#             csv_file, delimiter=","
#         )  # write CSV file using comma as delimiter

#         # Write the header row
#         csv_writer.writerow(headers)

#         # Write the remaining data rows
#         for row in csv_data:
#             csv_writer.writerow(row)

#     return csv_data
