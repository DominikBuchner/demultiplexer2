import glob, os
from pathlib import Path


def common_prefix_length(word1: str, word2: str) -> int:
    """Function to return the length of a common prefix of the words
    e.g. FileA_r1, FileA_r2 --> len(FileA_r) --> 7

    Args:
        word1 (str): Input word one.
        word2 (str): Input word two.

    Returns:
        int: Length of common prefix.
    """
    return len(os.path.commonprefix([word1, word2]))


def find_file_pairs(path_list: list) -> tuple:
    """Function to find all file pairs with matching names in a list of paths.
    File pairs are matched by the longest prefix match, e.g. if the input is
    [nameA_r1, nameA_r2, nameB_r1, nameB_r2] the output is
    [(nameA_r1, nameA_r2), (nameB_r1, nameB_r2)], []. Input strings are grouped by length
    so this will match any common extension. Any name that does not get a match
    will be returned as a list of single names.

    Args:
        path_list (list): List of file paths.

    Returns:
        tuple: List of matching pairs, list of single names without match.
    """
    # convert all string to Path object
    path_list = [Path(path) for path in path_list]

    # collect pairs and singles here
    file_pairs = []
    singles = []

    # return if no input is given
    if not path_list:
        return file_pairs, singles

    # group the input by length before starting to increase performance
    path_lengths = sorted(set([len(path.stem) for path in path_list]))
    grouped_paths = [
        [path for path in path_list if len(path.stem) == length]
        for length in path_lengths
    ]

    # loop over the different lengths groups
    for path_list in grouped_paths:
        # run until all paths were checked
        while path_list:
            # pop one path from the list to start comparison
            current_path = path_list.pop(0)

            # compute longest match with all other pairs, use stems to speed up the search
            longest_matches = [
                common_prefix_length(current_path.stem, path.stem) for path in path_list
            ]

            # compute maximum match and count how often it appears
            maximum_match = max(longest_matches)
            maximum_match_count = longest_matches.count(maximum_match)

            # check if there is a single longest match
            # a file pair has been found
            if maximum_match_count == 1:
                match_index = longest_matches.index(maximum_match)
                # retrieve the match from the list
                match = path_list.pop(match_index)
                file_pairs.append((current_path, match))
            else:
                singles.append(current_path)

    return file_pairs, singles


def main(path_list: list) -> tuple:
    """Function to control the pair finding script.

    Args:
        path_list (list): List of paths e.g. created with glob.glob.

    Returns:
        tuple: All file pairs found and all single files found.
    """
    file_pairs, singles = find_file_pairs(path_list)

    # return
    return file_pairs, singles
