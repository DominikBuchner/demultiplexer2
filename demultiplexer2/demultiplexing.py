import datetime, dnaio, duckdb
import pandas as pd
import numpy as np
from tqdm import tqdm
from pathlib import Path
from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
from demultiplexer2.create_tagging_scheme import collect_primerset_information
from joblib import Parallel, delayed
import pyarrow as pa
import pyarrow.parquet as pq


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


def convert_to_parquet(
    forward_path,
    reverse_path,
    forward_file,
    reverse_file,
    forward_tag_length,
    reverse_tag_length,
    file_index,
    output_dir,
):
    # create a filename for the parquet output
    parquet_path = Path(output_dir).joinpath(f"{file_index}.parquet.snappy")

    # define the schema for writing
    schema = pa.schema(
        [
            ("file_forward", pa.string()),
            ("file_reverse", pa.string()),
            ("name_forward", pa.string()),
            ("name_reverse", pa.string()),
            ("sequence_forward", pa.string()),
            ("sequence_reverse", pa.string()),
            ("quality_forward", pa.string()),
            ("quality_reverse", pa.string()),
            ("tag_forward", pa.string()),
            ("tag_reverse", pa.string()),
        ]
    )

    with pq.ParquetWriter(parquet_path, schema, compression="snappy") as writer:
        with dnaio.open(forward_path, reverse_path, mode="r") as reader:
            # define the columns for the table
            columns = {
                "file_forward": [],
                "file_reverse": [],
                "name_forward": [],
                "name_reverse": [],
                "sequence_forward": [],
                "sequence_reverse": [],
                "quality_forward": [],
                "quality_reverse": [],
                "tag_forward": [],
                "tag_reverse": [],
            }

            batch_size = 250_000
            batch = {k: [] for k in columns.keys()}

            for idx, (fwd, rev) in enumerate(reader):
                batch["file_forward"].append(forward_file)
                batch["file_reverse"].append(reverse_file)
                batch["name_forward"].append(fwd.name)
                batch["name_reverse"].append(rev.name)
                batch["sequence_forward"].append(fwd.sequence)
                batch["sequence_reverse"].append(rev.sequence)
                batch["quality_forward"].append(fwd.qualities)
                batch["quality_reverse"].append(rev.qualities)
                batch["tag_forward"].append(fwd.sequence[:forward_tag_length])
                batch["tag_reverse"].append(rev.sequence[:reverse_tag_length])

                if (idx + 1) % batch_size == 0:
                    writer.write_table(pa.table(batch, schema=schema))
                    batch = {k: [] for k in columns.keys()}

            # write remaining rows
            if any(len(v) > 0 for v in batch.values()):
                writer.write_table(pa.table(batch, schema=schema))


def build_database(output_dir):
    # gather all parquet files
    parquet_path = Path(output_dir).joinpath("*.parquet.snappy")

    # build the database name
    read_database = Path(output_dir).joinpath("read_database.duckdb")

    # open the connection
    read_database_connection = duckdb.connect(read_database)

    # ingest the parquet
    read_database_connection.execute(
        f"""
        CREATE OR REPLACE TABLE read_database AS
        SELECT * FROM read_parquet('{parquet_path}')
        """
    )

    read_database_connection.close()


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

    # user output
    print(
        "{}: Primerset successfully loaded.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
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

    # user output
    print(
        "{}: Tagging scheme successfully loaded.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # user output
    print(
        "{}: Building read database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # extract the tag length of the forward / reverse tags
    fwd_tag_length, rev_tag_length = len(extended_tags[0][0]), len(extended_tags[0][1])

    # # prepare a list of delayed tasks
    # tasks = [
    #     delayed(convert_to_parquet)(
    #         row["forward file path"],
    #         row["reverse file path"],
    #         row["forward file name"],
    #         row["reverse file name"],
    #         fwd_tag_length,
    #         rev_tag_length,
    #         index,
    #         output_dir,
    #     )
    #     for index, row in updated_tagging_scheme.iterrows()
    # ]

    # # run in parallel
    # Parallel(n_jobs=-1)(tasks)

    # transform parquet to duckdb database
    build_database(output_dir)
