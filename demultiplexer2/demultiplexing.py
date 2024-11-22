import pandas as pd
from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
from demultiplexer2.create_tagging_scheme import collect_primerset_information


def extend_ambiguous_dna(seq: str) -> list:
    """Returns a list of all possible sequences given DNA input with ambiguous bases.

    Args:
        seq (str): DNA sequence.

    Returns:
        list: List of all possible combinations of the sequence.
    """
    d = ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


def extend_tags(tags_used: list, forward_primer: str, reverse_primer: str) -> dict:
    """Function to extend a given list of tag tuples with the primer sequence so that
    a) all tags are of the same length
    b) tags are unique

    Args:
        tags_used (list): List of tag tuples from the tagging scheme
        forward_primer (str): Forward primer used in the dataset
        reverse_primer (str): Reverse primer used in the dataset

    Returns:
        dict: Dictionary with old tag pairs as keys and extended tag pairs as values. Can be multiple tag pairs if ambiguous DNA is detected
    """


def update_tagging_scheme(tag_information: object, tagging_scheme_path: str) -> object:
    """Function to read the tagging scheme update it with primer sequences instead of names

    Args:
        tag_information (object): Dataframe holding the primerset.
        tagging_scheme_path (str): Path to the tagging scheme.

    Returns:
        Object: Dataframe with the updated tagging scheme.
    """
    tagging_scheme = pd.read_excel(tagging_scheme_path)

    # extract primer names and sequences from tag information
    forward_tags = dict(
        zip(
            tag_information["name_forward_tag"], tag_information["sequence_forward_tag"]
        )
    )

    reverse_tags = dict(
        zip(
            tag_information["name_reverse_tag"], tag_information["sequence_reverse_tag"]
        )
    )

    # translate the dataframe columns to primer information
    sequence_header = [
        column_name.split("-") for column_name in tagging_scheme.columns[4:]
    ]

    sequence_header = [
        (forward_tags[column[0]], reverse_tags[column[1]]) for column in sequence_header
    ]

    # update the tagging scheme
    tagging_scheme = tagging_scheme.rename(
        columns=dict(zip(tagging_scheme.columns[4:], sequence_header))
    )

    return tagging_scheme


def main(primerset_path: str, tagging_scheme_path: str, output_dir: str):
    """Main function to run the demultiplexing.

    Args:
        primerset_path (str): Path to the primerset to be used.
        tagging_scheme_path (str): Path to the tagging scheme to be used.
        output_dir (str): Directory to write demultiplexed files to.
    """
    # read the primerset again to collect the sequence information
    forward_primer, reverse_primer, tag_information = collect_primerset_information(
        primerset_path
    )

    # extract the primers used in the tagging scheme, directly translate everything that is needed for demultiplexing
    # input paths, tagging information, output files
    updated_tagging_scheme = update_tagging_scheme(tag_information, tagging_scheme_path)

    print(updated_tagging_scheme)
