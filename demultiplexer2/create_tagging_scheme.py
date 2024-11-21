import glob, datetime
from demultiplexer2 import find_file_pairs
from pathlib import Path


def main(tagging_scheme_name: str, data_dir: str, primerset_path: str) -> None:
    # TODO
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
    # request input for combinations (form?)
    pass
