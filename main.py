import logging
import time
import warnings
import os
import shutil
from functions import (
    process_data,
    run_blast,
    process_blast_results,
    process_sequences_and_concatenate,
    generate_final_output,
    design_primers,
)


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
AAE5_loc = r"Support/AAE5.fna"


merged_df = process_data(assoc_path, r"Support/probes.xlsx")
logging.info("Successfully processed data.")


output_file_name = run_blast(YourName, merged_df, AAE5_loc, today)
logging.info("Successfully ran BLAST.")


final_df = process_blast_results(output_file_name, merged_df)
logging.info("Successfully processed BLAST results.")

chromosome_ids = ["NC_035107.1", "NC_035108.1", "NC_035109.1"]

final_df, concatenated_genome = process_sequences_and_concatenate(
    final_df, AAE5_loc, chromosome_ids
)
logging.info("Successfully processed sequences and concatenated.")

final_df = generate_final_output(final_df, concatenated_genome, YourName, today)
logging.info("Successfully generated output file.")

final_df = design_primers(final_df)
logging.info("Successfully designed primers.")

final_df.to_csv(f"{YourName}_{today}_final_output.csv", index=False)
logging.info("Successfully wrote final output to file.")


output_folder = "Output"
# Move the files to the output folder
os.makedirs(output_folder, exist_ok=True)
shutil.move(
    f"{YourName}_{today}_final_output.csv",
    f"{output_folder}/{YourName}_{today}_final_output.csv",
)
shutil.move(
    f"{YourName}_{today}_input.fasta", f"{output_folder}/{YourName}_{today}_input.fasta"
)
shutil.move(output_file_name, f"{output_folder}/{YourName}_{today}_output.xml")

shutil.move(
    f"{YourName}_{today}_final_output.fasta",
    f"{output_folder}/{YourName}_{today}_final_output.fasta",
)

logging.info("Successfully moved files to output folder.")
