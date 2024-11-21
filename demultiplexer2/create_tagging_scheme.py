import glob, datetime
import pandas as pd
import numpy as np
from demultiplexer2 import find_file_pairs
from pathlib import Path


def collect_primerset_information(primerset_path: str) -> tuple:
    """Function to parse the primerset file and return as a dataframe to chose from.
    Args:
        primerset_path (str): Path to the primerset file to be used.

    Returns:
        tuple: Forward primer name, reverse primer name and formatted primerset for user output.
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

        return None, None, None

    # extract the primer information
    forward_primer = (
        general_information["Forward primer (5' - 3')"].replace(np.nan, "").values[0]
    )
    reverse_primer = (
        general_information["Reverse primer (5' - 3')"].replace(np.nan, "").values[0]
    )

    # handle empty primer field accordingly
    if not forward_primer or not reverse_primer:
        print(
            "{}: Please add primer information to your primerset.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        return None, None, None

    # extract the tagging information
    forward_tags = pd.read_excel(
        primerset_path, sheet_name="forward_tags", index_col=0
    ).rename(columns={"name": "name_forward_tag"})
    reverse_tags = pd.read_excel(
        primerset_path, sheet_name="reverse_tags", index_col=0
    ).rename(columns={"name": "name_reverse_tag"})

    # concat both tags for easy to read output
    tag_information = pd.concat(
        [forward_tags, reverse_tags],
        axis=1,
    )

    # return primer names and tags
    return forward_primer, reverse_primer, tag_information


def request_combinations(tag_information: object) -> list:
    """Function to request the used combinations from the user.

    Args:
        tag_information (object): Tag information generated from the primerset.

    Returns:
        list: List of all tag combinations to be used to generate the tagging scheme.
    """
    print(
        "{}: All primer combinations used in your dataset are required for generating the tagging scheme.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    print(tag_information.to_string)

    return "test"


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
    forward_primer, reverse_primer, tag_information = collect_primerset_information(
        Path(primerset_path)
    )

    # request input for combinations (form?)
    combinations_to_use = request_combinations(tag_information)
