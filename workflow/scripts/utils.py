import numpy as np


# TODO vectorize operations.


def make_is_paired_fn(file_path):
    """
    Return a function that outputs 0 for an unpaired site and 1 for a basepaired site.

    Parameters:
        file_path (str): The path to a file with secondary rna structure data.
    """
    with open(file_path) as f:
        lines = [line.rstrip().split() for line in f]
    paired = np.array([[int(x[0]), int(x[4])] for x in lines[1:]])
    paired_dict = dict(zip(paired[:, 0], paired[:, 1]))

    def is_paired(site):
        if site not in paired_dict or paired_dict[site] == 0:
            return 0
        else:
            return 1

    return is_paired


def mut_str_to_dict(mut_str):
    """
    Return a dictionary of mutations from a string mutations. Keys are site indices.
    Values are nucleotides after mutation.

    Parameters:
        mut_str (str): A whitespace delimited string of mutations. Individual mutations
        are like "A25G" for A->G at site 25.
    """

    return {int(mut[1:-1]): mut[-1] for mut in mut_str.split()}


def motif_for_site(site, founder_muts, ref_seq):
    """
    Return the 3mer motif centered at the given site of the reference sequence after
    applying the mutations.
    """
    muts = founder_muts
    which_nuc = lambda s: ref_seq[s - 1] if s not in muts else muts[s]
    return "".join(map(which_nuc, range(site - 1, site + 2)))


def motif_for_subtree_counts_row(row, ref_seq):
    """
    Return the 3mer motif centered at row.site of the reference sequence after
    applying row.mutations.
    """
    return motif_for_site(row.site, row.mutations, ref_seq)
