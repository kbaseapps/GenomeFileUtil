"""
This file contains miscellaneous utilities that are not tied to the concept of a genome.
"""


def validate_lists_have_same_elements(l1, l2):
    """
    Given two lists/sets of values (from different sources), verify that they match up.
    
    Useful for comparing ids in GFF and Fasta files, or across Genomes and Assemblies.
    """
    diff = set(l1) ^ (set(l2))  # get the symmetric difference of the sets
    # check if all ids are shared
    return len(diff) == 0


def get_int(putative_int, name, minimum=1):
    """
    Validates and returns an integer value.

    This function checks whether the provided value is an integer and if it meets the specified minimum value.
    If the checks are not passed, it raises a `ValueError` with a descriptive message.

    Args:
        putative_int (int or None): The value to be validated and returned. If `None`, it will be returned as is.
        name (str): A descriptive name for the value being checked. This is used in error messages.
        minimum (int, optional): The minimum acceptable value for `putative_int`. Defaults to 1.

    Returns:
        int: The validated integer if all checks are passed.

    Raises:
        ValueError: If `putative_int` is not an integer, or if it is less than `minimum`.
    """
    if putative_int is not None:
        if type(putative_int) is not int:
            raise ValueError(f"{name} must be an integer, got: {putative_int}")
        if putative_int < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    return putative_int
