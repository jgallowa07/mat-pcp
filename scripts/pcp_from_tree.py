#!/usr/bin/env python
from pandas import DataFrame as DF
from collections import defaultdict
import bte
from time import time
import argparse

parser = argparse.ArgumentParser(description="Write out PCPs from a tree.")
parser.add_argument("tree_path", help="Path to the tree file")
# parser.add_argument("coding_site_path", help="Path to the coding site file")
parser.add_argument("fasta_path", help="Path to the fasta file")
parser.add_argument("gtf_path", help="Path to the gtf file")
parser.add_argument("all_pcps_path", help="Path to write out all pcps")

def load_reference_from_fasta(reference_file):
    """
    Provided a path to a fasta file containing a single sequence, load that sequence
    and return it as a string.
    """
    with open(reference_file, "r") as fh:
        next(fh)
        reference_seq = "".join(map(str.strip, fh))
    return reference_seq

def apply_muts(sequence, muts):        
    """                                                                       
    Return the sequence after applying mutations.
                                       
    Args:                                                                     
        sequence (str): The sequence.
        muts (dict): A dictionary of mutations. The keys are integers for the site of a
            mutation (site indices start at 1); the values are strings for the resulting
            nucleotide at the site. For example, muts[25] = 'G' indicates a mutation
            resulting in 'G' at site 25 (which is at index 24 of sequence).
    """
    if len(muts) == 0:
        return sequence
    sites = sorted(muts.keys())
    substrings = [sequence[: sites[0] - 1]]
    for start, stop in zip(sites[:-1], sites[1:]):
        tobase = muts[start]
        seq = sequence[start : stop - 1]
        substrings.append(tobase)
        substrings.append(seq)
    substrings.append(muts[sites[-1]])
    substrings.append(sequence[sites[-1] :])

    return "".join(substrings)

class PCPHelper:
    """Helper class to write out PCP's passing a filter from a tree."""

    def __init__(
        self,
        tree_path,
        fasta_path,
        gtf_path,
    ):

        self.tree = bte.MATree(tree_path)
        self.translations = defaultdict(list, self.tree.translate(gtf_path, fasta_path))
        self.nodes = self.tree.depth_first_expansion()
        self.int_nodes = [n for n in self.nodes if not n.is_leaf()]
        self.n_nodes = len(self.nodes)
        self.n_int_nodes = len(self.int_nodes)
        self.ref_seq = load_reference_from_fasta(fasta_path)


    def fails_filters(self, nt_muts, codon_muts):
        fail = len(nt_muts) > 4
        fail |= len(nt_muts) == 0
        fail |= sum(self.ref_seq[int(mut[1:-1]) - 1] == mut[-1] for mut in nt_muts) > 1
        fail |= len({(mut.gene, mut.aa_index) for mut in codon_muts}) < len(codon_muts)
        return fail

    def count_mutations_from_parent(self, parent):
        parent_node_haplotype = self.tree.get_haplotype(parent.id)
        p_muts = {int(mut[1:-1]): mut[-1] for mut in parent_node_haplotype}
        parent_seq = apply_muts(self.ref_seq, p_muts)
        pcps = (parent.id, parent_seq, [])

        n_filtered = 0
        for node in parent.children:
            nt_muts = node.mutations
            codon_muts = self.translations[node.id]
            if self.fails_filters(nt_muts, codon_muts):
                continue

            n_filtered += 1
            pcps[2].append((node.id, nt_muts))
        return n_filtered, pcps

    def count_mutations_on_tree(self):
        """
        # Cycle over nodes and record codon mutation counts
        ...returns int, Counter, Counter
        """
        # Use counters rather than a dataframe, since summing together many.
        n_filtered = 0
        all_pcps = []
        t0 = -time()
        for i, node in enumerate(self.int_nodes, 1):
            count, pcps = self.count_mutations_from_parent(node)
            n_filtered += count
            if count != 0:
                all_pcps.append(pcps)

            if i % 1000 == 0:
                print(f"{i} nodes processed in {t0 + time():0.1f} seconds")
        all_pcps_df = self.pcp_list_to_df(all_pcps)
        return n_filtered, all_pcps_df

    # TODO, currently, this slices out only spike.
    def pcp_list_to_df(self, pcps, gene_slice_boundaries=(21563, 25381)):
        def mut_dict(muts):
            return {int(mut[1:-1]): mut[-1] for mut in muts}
        def child_seq(seq, muts):
            return apply_muts(seq, mut_dict(muts))
        def slice_seq(seq):
            return (
                seq 
                if gene_slice_boundaries is None 
                else seq[gene_slice_boundaries[0] - 1 : gene_slice_boundaries[1]]
            )
        
        cols = ["parent_name", "child_name", "parent", "child", "branch_length"]
        rows = (
            (p_id, c_id, slice_seq(p_seq), slice_seq(child_seq(p_seq, c_muts)), len(c_muts))
            for (p_id, p_seq, child_entries) in pcps
            for c_id, c_muts in child_entries
        )
        # TODO, are we not wasting memory by converting to DF instead of writing out to the
        # csv directly? 
        the_df = DF(rows, columns=cols)
        return the_df

    def make_and_write_pcps(
        self,
        all_pcps_path,
    ):

        print("Getting pcp's from tree ...")
        n_filtered, all_pcps_df = self.count_mutations_on_tree()
        print(f"Number of nodes passing filter: {n_filtered}")
        all_pcps_df.to_csv(all_pcps_path)

        return None


def write_pcps(
    tree_path,
    # coding_site_path,
    fasta_path,
    gtf_path,
    all_pcps_path,
):
    p1 = tree_path, fasta_path, gtf_path
    PCPHelper(*p1).make_and_write_pcps(all_pcps_path)
    return None


if __name__ == "__main__":

    args = parser.parse_args()
    write_pcps(
        args.tree_path,
        args.fasta_path, 
        args.gtf_path, 
        args.all_pcps_path
    )
