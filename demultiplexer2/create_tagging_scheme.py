import glob, datetime, sys, os
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

        sys.exit()

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
        sys.exit()

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
    print(
        "{}: Press any key to cotinue.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # wait for any key
    input()

    # repeat the information from the primerset for easy input
    print(tag_information.to_string() + "\n")

    # more user output
    print(
        "{}: Please indicate all COMBINATIONS by providing the respective line number in the format fwd-rev,fwd-rev...".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    print(
        "{}: Example: 1-1,2-2,2-3,4-4,5-5".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # ask for input
    combinations = input("Combinations used: ")

    # parse the input
    combinations = combinations.split(",")
    combinations = [tuple(combination.split("-")) for combination in combinations]
    combinations = [(int(comb[0]), int(comb[1])) for comb in combinations]

    # translate to primer names for clean output
    fwd_tag_names = dict(
        zip(tag_information.index, tag_information["name_forward_tag"])
    )
    rev_tag_names = dict(
        zip(tag_information.index, tag_information["name_reverse_tag"])
    )

    combinations = [
        (fwd_tag_names[combination[0]], rev_tag_names[combination[1]])
        for combination in combinations
    ]

    # return the combinations
    return combinations


def create_tagging_scheme_file(
    tagging_scheme_name: str, file_pairs: list, combinations_to_use: list
):
    """Function that creates a tagging scheme file and saves it to excel.

    Args:
        tagging_scheme_name (str): Name of the tagging scheme.
        file_pairs (list): List of all file pairs to demultiplex.
        combinations_to_use (list): Primer combinations that were used to generate the dataset.
    """
    # extract file names from the file pairs
    file_pair_names = [(pair[0].name, pair[1].name) for pair in file_pairs]

    # add the names to the paths to have all file data
    file_data = [
        (file_pair + file_pair_name)
        for file_pair, file_pair_name in zip(file_pairs, file_pair_names)
    ]

    # add the file data to a dataframe
    tagging_scheme = pd.DataFrame(
        file_data,
        columns=[
            "forward file path",
            "reverse file path",
            "forward file name",
            "reverse file name",
        ],
    )

    # generate the primer columns
    primer_columns = [
        "{}-{}".format(tag_fwd, tag_rev) for tag_fwd, tag_rev in combinations_to_use
    ]

    primer_columns = pd.DataFrame([], columns=primer_columns)

    # concat the final tagging scheme
    tagging_scheme = pd.concat([tagging_scheme, primer_columns], axis=1)

    # build the save path
    savepath = Path(os.getcwd()).joinpath(
        "{}_tagging_scheme.xlsx".format(tagging_scheme_name)
    )

    # save the scheme
    tagging_scheme.to_excel(savepath, index=False)

    # give user output
    print(
        "{}: Tagging scheme saved at {}.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), savepath
        )
    )


def main(tagging_scheme_name: str, data_dir: str, primerset_path: str) -> None:
    """Main function to create a tagging scheme.

    Args:
        tagging_scheme_name (str): The name of the tagging scheme that is used for saving it.
        data_dir (str): The directory where the files to demultiplex are stored.
        primerset_path (str): The path to the primerset to use for demultiplexing.
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

    # request input for combinations
    combinations_to_use = request_combinations(tag_information)

    # create the input file for the tagging scheme
    create_tagging_scheme_file(tagging_scheme_name, file_pairs, combinations_to_use)
