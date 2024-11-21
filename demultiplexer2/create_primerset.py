import datetime
import pandas as pd
import importlib.resources as resources
from pathlib import Path


def create_primerset(primerset_name: str, n_primers: int):
    """Function to create a primerset in an Excel sheet.

    Args:
        primerset_name (str): Name of the primerset.
        n_primers (int): Number of primers in the primerset.
    """
    # generate a path to save the primerset. This is wherever demultiplexer2 is installed.
    primerset_savepath = Path(resources.files(__package__)).joinpath(
        "data/primersets/{}_primerset.xlsx".format(primerset_name)
    )

    if primerset_savepath.is_file():
        print(
            "{}: A primerset with this name already exists. Please choose another name.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        return None

    # write output to excel
    with pd.ExcelWriter(primerset_savepath, mode="w", engine="openpyxl") as writer:
        # create the general information
        general_information = pd.DataFrame(
            [
                [
                    "Put the forward primer used for amplification here.",
                    "Put the reverse primer used for amplification here.",
                ]
            ],
            columns=["Forward primer (5' - 3')", "Reverse primer (5' - 3')"],
        )

        general_information.to_excel(
            writer, sheet_name="general_information", index=False
        )

        # create the forward sets
        forward_tags = pd.DataFrame(
            [["", ""]] * n_primers, columns=["name", "sequence"]
        )
        forward_tags.index = forward_tags.index + 1

        # since forward and reverse tag sheets are identical they can be written from the same dataset
        for sheet_name in ["forward_tags", "reverse_tags"]:
            forward_tags.to_excel(writer, sheet_name=sheet_name, index=True)

    print(
        "{}: Primerset was saved at:\n{}.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), primerset_savepath
        )
    )
