import logging
import time
import warnings
from new_functions import process_data
import os
from Bio import SeqIO
import pandas as pd
from snapgene_reader import snapgene_file_to_seqrecord


# #########USER INPUTS##########

# Path to the association results file
assoc_path = r"Inputs/associaton_results.assoc"

YourName = "Ryan_Test"

today = time.strftime("%Y-%m-%d")

####################################


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s:%(message)s",
)
warnings.filterwarnings("ignore", category=FutureWarning)
logging.info("Successfully read assoc file.")
# AAE5_loc = r"Support/AAE5.fna"
# input_list = [
#     "AX-93263480",
#     "AX-93258685",
#     "AX-93239662",
#     "AX-93224995",
#     "AX-93246299",
#     "AX-93257379",
#     "AX-93217744",
#     "AX-93252573",
#     "AX-93258936",
#     "AX-93257863",
#     "AX-93239566",
#     "AX-93250737",
#     "AX-93243723",
#     "AX-93251560",
#     "AX-93251321",
#     "AX-93253587",
#     "AX-93226405",
#     "AX-93214935",
# ]
# logging.info("length of input list: " + str(len(input_list)))

# merged_df = process_data(assoc_path, r"Support/probes.xlsx", input_list)
# logging.info("Successfully processed data.")

candidate_list = pd.read_csv("cand_list.csv")

# candidate_list = merged_df


# Function to convert .dna to SeqRecord
def convert_to_seqrecord(directory):
    sequences = []
    for filename in os.listdir(directory):
        if filename.endswith(".dna"):
            file_path = os.path.join(directory, filename)
            seqrecord = snapgene_file_to_seqrecord(file_path)
            sequences.append(seqrecord)
    return sequences


# Use the function on your directories
orlando_seqs = convert_to_seqrecord("Orlando/")
liverpool_seqs = convert_to_seqrecord("Liverpool/")

expected_orl = len(
    [
        entry
        for entry in os.listdir("Orlando/")
        if os.path.isfile(os.path.join("Orlando/", entry))
    ]
)
expected_lvp = len(
    [
        entry
        for entry in os.listdir("Liverpool/")
        if os.path.isfile(os.path.join("Liverpool/", entry))
    ]
)


# Create a function that checks for a match allowing for a certain tolerance
def match_with_tolerance(seq, target, tolerance):
    seq_len = len(seq)
    for i in range(len(target) - seq_len + 1):
        sub = target[i : i + seq_len]
        differences = sum(a != b for a, b in zip(seq, sub))
        if differences <= tolerance:
            return True
    return False


# Iterate over the sequences and check for matches
matches_orlando = []
matches_liverpool = []
tolerance = 0

mismatches_orlando = []
mismatches_liverpool = []

for i, row in candidate_list.iterrows():
    match_found_orlando = False
    match_found_liverpool = False

    for seq in orlando_seqs:
        if match_with_tolerance(seq, row["Flank with Orlando"], tolerance):
            matches_orlando.append(row["SNP"])
            match_found_orlando = True

    for seq in liverpool_seqs:
        if match_with_tolerance(seq, row["Flank with LVP"], tolerance):
            matches_liverpool.append(row["SNP"])
            match_found_liverpool = True

    if not match_found_orlando:
        mismatches_orlando.append(row["SNP"])

    if not match_found_liverpool:
        mismatches_liverpool.append(row["SNP"])

matches_orlando = set(matches_orlando)
matches_liverpool = set(matches_liverpool)

print(f"Number of matches for Liverpool: {len(matches_liverpool)}")
print(f"Number of expected matches for Liverpool: {expected_lvp}")
print(f"SNPS with Liverpool mismatches: {mismatches_liverpool}")
print("")
print(f"Number of matches for Orlando: {len(matches_orlando)}")
print(f"Number of expected matches for Orlando: {expected_orl}")
print(f"SNPS with Orlando mismatches: {mismatches_orlando}")
