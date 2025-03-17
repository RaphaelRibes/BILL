from Bio import GenBank as gb
import json

def get_orfs():
    """
    Extract ORFs (Open Reading Frames) from a GenBank file.

    Returns:
        dict: A dictionary where keys are ORF names and values are tuples containing the start and end positions of the ORFs.
    """
    # Open the GenBank file
    with open('data/Cyprinid herpesvirus 3 strain KHV-U, complete genome.gb') as f:
        record = gb.read(f)  # Read the record using the GenBank parser

    orfs = {}
    for feature in record.features:  # Iterate over the features
        complement = False

        if feature.key == 'gene':  # If the feature is a gene (ORF)
            orf = feature.qualifiers[0].value.replace('"', '').replace('CyHV3_', '')

            if 'complement' in feature.location:  # Check if the feature is in the complement strand
                loc = feature.location.replace('complement(', '').replace(')', '').split('..')
                complement = True
            else:
                loc = feature.location.split('..')

            orfs[orf] = (int(loc[0]), int(loc[1]), complement)  # Store the ORF information

    return orfs

def get_targeted_orfs(start: int, end: int):
    """
    Retrieve ORFs (Open Reading Frames) that overlap with a given range.

    Args:
        start (int): The start position of the target range.
        end (int): The end position of the target range.

    Returns:
        list: A list of tuples, each containing the ORF name, start position, and end position.
    """
    final = []
    for orf, (orf_start, orf_end, complement) in get_orfs().items():
        # Check if the ORF overlaps with the given range
        if (orf_start >= start or start <= orf_end) and end >= orf_start:
            final.append((orf, orf_start, orf_end, complement))
    return final
