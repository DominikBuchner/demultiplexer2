import glob, datetime
import pandas as pd
from demultiplexer2 import find_file_pairs
from pathlib import Path


def collect_primerset_information(primerset_path: str) -> object:
    """Function to parse the primerset file and return as a dataframe to chose from.
    Args:
        primerset_path (str): Path to the primerset file to be used.

    Returns:
        object: Formatted primerset for user output.
    """
    # parse the specific primers first
    try:
        general_information = pd.read_excel(
            primerset_path, sheet_name="general_information"
        )
    except (FileNotFoundError, ValueError):
        print(
            "{}: Please select a valid input for the primerset path.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        return None

    forward_primer = general_information["Forward primer (5' - 3')"].item()
    reverse_primer = general_information["Reverse primer (5' - 3')"].item()


def main(tagging_scheme_name: str, data_dir: str, primerset_path: str) -> None:
    """Main function to create a tagging scheme.

    Args:
        tagging_scheme_name (str): The name of the tagging scheme that is used for saving it.
        data_dir (str): The directory where the files to demultiplex are stored.
        primerset_path (str): The path to the primerset to use for demultiplexing.

    Returns:
        _type_: _description_
    """
    # collect all file pairs from input folder, give warning if unpaired files are found
    all_files = glob.glob(str(Path(data_dir).joinpath("*.fastq.gz")))
    file_pairs, singles = find_file_pairs.find_file_pairs(all_files)

    # give error message file cannot be matched
    if singles:
        for file in singles:
            print("{}: {}.".format(datetime.datetime.now().strftime("%H:%M:%S"), file))

        print(
            "{}: Found files that could not be matched as pairs. Please check your input.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        return None

    if not file_pairs:
        print(
            "{}: No file pairs where found. Please check your input.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        return None

    # collect primers / tags / names from primer set
    primerset_information = collect_primerset_information(Path(primerset_path))
    # request input for combinations (form?)
    pass
