import gzip
import pandas as pd
import numpy as np
from pathlib import Path
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


def check_tag_distances(tag_list: list) -> bool:
    """Function to check if all tags are unique. If this is true for forward and reverse tag, distance = 2 is met.

    Args:
        tag_list (list): List of DNA tags.

    Returns:
        bool: True if all distances >= dist.
    """
    # extend the tag list by removing ambigouities
    extended_tag_list = []

    for tag in tag_list:
        extended_tag_list += extend_ambiguous_dna(tag)

    # check if all tags are unique
    if len(set(extended_tag_list)) == len(extended_tag_list):
        return True
    else:
        return False


def extend_by_one(tag_list: list, primer_list: list) -> tuple:
    """Function to extend all tags in a list of tags with on base from the beginning of the primer
    list of tags and primers have to have the same length

    Args:
        tag_list (list): List of tags
        primer_list (list): List of primers

    Returns:
        tuple: Tuple with the extended tag list and shortened primer list.
    """
    # extend the tags by one, shorten the primers by one
    for idx in range(len(tag_list)):
        tag_list[idx] += primer_list[idx][:1]
        primer_list[idx] = primer_list[idx][1:]

    return tag_list, primer_list


def extend_tags(
    updated_tagging_scheme: object, forward_primer: str, reverse_primer: str
) -> dict:
    """Function calculate unambiguous Tags from a given tagging scheme so that
    a) all tags are of the same length
    b) tags are unique

    Args:
        updated_tagging_scheme (object): Tagging scheme as dataframe.
        forward_primer (str): Forward primer used in the dataset
        reverse_primer (str): Reverse primer used in the dataset

    Returns:
        dict: Dictionary with old tag pairs as keys and extended tag pairs as values. Can be multiple tag pairs if ambiguous DNA is detected
    """
    # extract the tag pairs from the scheme header
    tag_pairs = updated_tagging_scheme.columns[4:]

    # generate a list of primer lists of the same length for extending the tags
    primer_pairs = [(forward_primer, reverse_primer) for _ in tag_pairs]

    # find the length of the longest tag
    max_tag_length = max(
        [len(tag_pair[0]) for tag_pair in tag_pairs]
        + [len(tag_pair[1]) for tag_pair in tag_pairs]
    )

    # extend the tags to the maximum length
    # split the tuples into individual lists
    forward_tags, reverse_tags = (
        [tag_pair[0] for tag_pair in tag_pairs],
        [tag_pair[1] for tag_pair in tag_pairs],
    )

    # split the primers into individual lists
    forward_primers, reverse_primers = (
        [primer[0] for primer in primer_pairs],
        [primer[1] for primer in primer_pairs],
    )

    # extend the tags, shorten the primer sequences until all have the same length
    for idx in range(len(forward_tags)):
        # calculate the difference
        length_difference = max_tag_length - len(forward_tags[idx])
        # if there is a difference
        if length_difference:
            forward_tags[idx] += forward_primers[idx][:length_difference]
            forward_primers[idx] = forward_primers[idx][length_difference:]

        # calculate the difference
        length_difference = max_tag_length - len(reverse_tags[idx])
        # if there is a difference
        if length_difference:
            reverse_tags[idx] += reverse_primers[idx][:length_difference]
            reverse_primers[idx] = reverse_primers[idx][length_difference:]

    # check distances within tags --> is this unambiguous?
    while not check_tag_distances(forward_tags):
        forward_tags, forward_primers = extend_by_one(forward_tags, forward_primers)

    while not check_tag_distances(reverse_tags):
        reverse_tags, reverse_primers = extend_by_one(reverse_tags, reverse_primers)

    extended_tags = [
        (forward_tag, reverse_tag)
        for forward_tag, reverse_tag in zip(forward_tags, reverse_tags)
    ]

    return extended_tags


def generate_demultiplexing_data(updated_tagging_scheme: object) -> dict:
    """Creates a dict that holds all data needed for demultiplexing the files in the form of
    {(input_fwd, input_rev): {(tag_fwd, tag_rev) : file_name1}, (tag_fwd, tag_rev): file_name2...
     (input_fwd2, input_rev2): {(tag_fwd, tag_rev) : file_name3}, (tag_fwd, tag_rev): file_name4...
    }

    Args:
        updated_tagging_scheme (object): Dataframe with the updated tagging scheme

    Returns:
        dict: dict with demultiplexing data
    """
    # remove empty data from updated_tagging_scheme
    updated_tagging_scheme = updated_tagging_scheme.replace(np.nan, "")

    # extract the tag combinations from the column names
    tag_combinations = updated_tagging_scheme.columns[4:]

    # extend ambiguoities if there are any
    extended_combinations = {}

    for tag_combination in tag_combinations:
        fwd_tag, rev_tag = tag_combination[0], tag_combination[1]
        extended_fwd_tags, extended_rev_tags = extend_ambiguous_dna(
            fwd_tag
        ), extend_ambiguous_dna(rev_tag)

        # compute all combinations
        extended_combinations[tag_combination] = [
            (fwd_tag, rev_tag)
            for fwd_tag in extended_fwd_tags
            for rev_tag in extended_rev_tags
        ]

    # store all data needed for demultiplexing here
    demultiplexing_data = {}

    # go over the individual rows of the updated tagging scheme
    for idx, row in updated_tagging_scheme.iterrows():
        input_fwd, input_rev = row["forward file path"], row["reverse file path"]

        output_per_input = {}
        # loop over the tag combinations to retrieve the extended combinations
        for tag_combination in updated_tagging_scheme.columns[4:]:
            for combination in extended_combinations[tag_combination]:
                # only add data that from combinations connected to samples
                if row[tag_combination] != "":
                    output_per_input[combination] = row[tag_combination]

        # update the demultiplexing data
        demultiplexing_data[(input_fwd, input_rev)] = output_per_input

    return demultiplexing_data


def demultiplexing(
    demultiplexing_data_key: str, demultiplexing_data_value: dict, output_dir: str
):
    """Function to run the demultiplexing.

    Args:
        demultiplexing_data_key (str): Key e.g. input file pair from the demultiplexing data
        demultiplexing_data_value (dict): Value corresponsing to the key from the demultiplexing data, e.g. tag pair --> output
        output_dir (str): Directory to write to.
    """
    # generate a dict mapping sample names to actual output paths
    output_handles = {}

    for sample in demultiplexing_data_value.values():
        fwd_path = Path(output_dir).joinpath("{}_r1.fastq.gz".format(sample))
        rev_path = Path(output_dir).joinpath("{}_r2.fastq.gz".format(sample))

        # add paths to the output handles
        output_handles[sample] = (
            gzip.open(fwd_path, "wt", compresslevel=6),
            gzip.open(rev_path, "wt", compresslevel=6),
        )

    # count some basic statistics
    total_reads, unmatched_reads = 0, 0

    # extract the length of the fwd tags and reverse tags, can be done from a single tag, since all tags have the same length
    length_forward_tag, length_reverse_tag = (
        len(list(demultiplexing_data_value.keys())[0][0]),
        len(list(demultiplexing_data_value.keys())[0][1]),
    )

    print(length_forward_tag, length_reverse_tag)


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

    # extend tagging information to remove ambiguoity
    extended_tags = extend_tags(updated_tagging_scheme, forward_primer, reverse_primer)

    # update the extended tags into the updated tagging scheme
    updated_tagging_scheme = updated_tagging_scheme.rename(
        columns=dict(zip(updated_tagging_scheme.columns[4:], extended_tags))
    )

    # generate the data needed for demultiplexing as dict from the updated tagging scheme
    demultiplexing_data = generate_demultiplexing_data(updated_tagging_scheme)

    # TODO parallelize at the end
    for key in demultiplexing_data.keys():
        demultiplexing(key, demultiplexing_data[key], output_dir)
        break
