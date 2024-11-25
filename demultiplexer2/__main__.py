import argparse, datetime, sys
from importlib.metadata import version
from get_pypi_latest_version import GetPyPiLatestVersion
from demultiplexer2 import create_primerset, create_tagging_scheme, demultiplexing


def main() -> None:
    """Function to define the commandline interface"""
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=35)

    # define the parser
    parser = argparse.ArgumentParser(
        prog="demultiplexer2",
        description="A Python package to demultiplex Illumina reads by inline tags.",
        formatter_class=formatter,
    )

    # display help when no argument is called
    parser.set_defaults(func=lambda x: parser.print_help())

    # add the subparsers
    subparsers = parser.add_subparsers(dest="function")

    # add the create primerset parser, add arguments
    parser_create_primerset = subparsers.add_parser(
        "create_primerset", help="Create a primerset for demultiplexer2."
    )

    parser_create_primerset.add_argument(
        "--name",
        required=True,
        help="Define the name for the primerset to create.",
        type=str,
    )

    parser_create_primerset.add_argument(
        "--n_primers",
        required=True,
        help="Define the number of forward and reverse primers in the primerset.",
        type=int,
    )

    # add the create_tagging_scheme parser
    parser_create_tagging_scheme = subparsers.add_parser(
        "create_tagging_scheme", help="Create a tagging scheme for demultiplexer2."
    )

    parser_create_tagging_scheme.add_argument(
        "--name",
        required=True,
        help="Define the name for the tagging scheme to create",
        type=str,
    )

    parser_create_tagging_scheme.add_argument(
        "--data_dir",
        required=True,
        help="Path to the directory that contains the files to demultiplex.",
        type=str,
    )

    parser_create_tagging_scheme.add_argument(
        "--primerset_path",
        required=True,
        help="Path to the primerset to be used for demultiplexing.",
        type=str,
    )

    # add the demultiplex parser
    parser_demultiplex = subparsers.add_parser(
        "demultiplex", help="Start demultiplexing."
    )

    parser_demultiplex.add_argument(
        "--primerset_path",
        required=True,
        help="Path to the primerset to be used for demultiplexing.",
        type=str,
    )

    parser_demultiplex.add_argument(
        "--tagging_scheme_path",
        required=True,
        help="Path to the tagging scheme to be used for demultiplexing.",
        type=str,
    )

    parser_demultiplex.add_argument(
        "--output_dir",
        required=True,
        help="Directory to write the demultiplexed files to.",
        type=str,
    )

    add version control (TODO)
    current_version = version("demultiplexer2")
    obtainer = GetPyPiLatestVersion()
    latest_version = obtainer("demultiplexer2")

    # give a user warning if the latest version is not installed
    if current_version != latest_version:
        print(
            "{}: Your demultiplexer2 version is outdated. Consider updating to the latest version.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

    # add the version argument
    parser.add_argument(
        "--version", action="version", version=version("demultiplexer2")
    )

    # parse the arguments
    arguments = parser.parse_args()

    # print help if no argument is provided
    if len(sys.argv) == 1:
        arguments.func(arguments)
        sys.exit()

    # define behaviour for different arguments that might be called
    if arguments.function == "create_primerset":
        # create a new primerset
        create_primerset.create_primerset(arguments.name, arguments.n_primers)

    if arguments.function == "create_tagging_scheme":
        # create a new tagging scheme
        create_tagging_scheme.main(
            arguments.name, arguments.data_dir, arguments.primerset_path
        )

    if arguments.function == "demultiplex":
        # run the demultiplexing
        demultiplexing.main(
            arguments.primerset_path,
            arguments.tagging_scheme_path,
            arguments.output_dir,
        )


# run only if called as a top level script
if __name__ == "__main__":
    main()
