import csv


def convert_assoc_to_csv(assoc_file_path, csv_file_path):
    csv_data = []
    with open(assoc_file_path, "r") as assoc_file:
        assoc_reader = csv.reader(
            assoc_file, delimiter="\t"
        )  # assuming your original file is tab-separated

        for row in assoc_reader:
            csv_data.append(row)

    # Split the first row into separate columns
    headers = csv_data[0][0].split()

    # Remove the header row from csv_data
    csv_data = csv_data[1:]

    with open(csv_file_path, "w", newline="") as csv_file:
        csv_writer = csv.writer(
            csv_file, delimiter=","
        )  # write CSV file using comma as delimiter

        # Write the header row
        csv_writer.writerow(headers)

        # Write the remaining data rows
        for row in csv_data:
            csv_writer.writerow(row)

    return csv_data
