import json
from typing import List, Dict

import numpy as np
import pandas as pd

from .taxon_manager import TaxonManager


def get_peptide_count_per_tax_id(proteins_per_taxon):
    """
    Convert tab-separated taxon-protein counts to a dictionary.

    Parameters
    ----------
    proteins_per_taxon: str,
        Path to file that contains the counts.

    """
    protein_counts_per_taxid = {}

    with open(proteins_per_taxon, "rb") as f:
        for line in f:
            # The file is encoded
            line_str = line.decode("utf-8").strip()
            # First column represents the TaxID, the second the count of peptides that are associated with that TaxID
            taxid, count = map(int, line_str.split("\t"))
            protein_counts_per_taxid[taxid] = int(count)

    return protein_counts_per_taxid


def get_lineage_at_specified_rank(tax_ids: List[int], taxa_rank: str) -> List[int]:
    """
    Returns the taxon ID of the specified rank in the lineage for all taxa ID given as argument.

    For example, given a taxon ID at strain level and "species" as value for the taxa_rank argument, this function
    will return the taxon ID at species level for the input taxon ID.

    Parameters
    -----
    tax_ids: [int]
         List of taxon_ids to get the lineage of
    taxa_rank:
         Rank at which you want to pin the taxa
    """

    # Get the full lineage for all given taxon IDs from Unipept
    lineages = TaxonManager.get_lineages_for_taxa(tax_ids)

    # Get the index of the NCBI rank that we're interested in. This index is required to extract the taxon IDs from the
    # correct place in the lineage.
    rank_idx = TaxonManager.NCBI_RANKS.index(taxa_rank)

    return [lineages[tax][rank_idx] for tax in tax_ids]


def perform_taxa_weighing(
    unipept_responses: List[any],
    pept_scores: Dict[str, Dict[str, float | int]],
    max_taxa,
    taxa_rank="species"
):
    """
    Weight inferred taxa based on their (1) degeneracy and (2) their proteome size.
    Parameters
    ----------
    unipept_responses: List[any]
        Peptide counts that have already been processed by Unipept before.
    pept_scores: Dict[str, Dict[str, float | int]]
        Dictionary that maps each peptide string onto an object containing the score associated to this peptide and the
        psm count.
    max_taxa: int
        Maximum number of taxa to include in the final graphical model.
    taxa_rank: str
        NCBI rank at which the Peptonizer analysis should be performed.

    Returns
    -------
    dataframe
        Top scoring taxa

    """
    print("Parsing Unipept responses from disk...")

    # Convert a JSON object into a Pandas DataFrame
    # record_path Parameter is used to specify the path to the nested list or dictionary that you want to normalize
    print("Normalizing peptides and converting to dataframe...")
    unipept_frame = pd.json_normalize(unipept_responses)
    # Merge psm_score and number of psms
    unipept_frame = pd.concat(
        [
            unipept_frame,
            pd.json_normalize(unipept_frame["sequence"].map(pept_scores)),
        ],
        axis=1,
    )

    # Score the degeneracy of a taxa, i.e.,
    # how conserved a peptide sequence is between taxa.
    # map all taxids in the list in the taxa column back to their taxid at species level (or the rank specified by the user)
    print("Started mapping all taxon ids to the specified rank...")
    unipept_frame["HigherTaxa"] = unipept_frame.apply(
        lambda row: get_lineage_at_specified_rank(row["taxa"], taxa_rank), axis=1
    )

    # Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with, exponentiated by 3
    print("Started dividing the number of PSMS of a peptide by the number the peptide is associated with...")
    unipept_frame["weight"] = unipept_frame["psms"].div(
        [len(element) ** 3 for element in unipept_frame["HigherTaxa"]]
    )
    mask = [len(element) == 1 for element in unipept_frame["HigherTaxa"]]
    unique_psm_taxa = set(i[0] for i in unipept_frame["HigherTaxa"][mask])
    unipept_frame = unipept_frame.explode("HigherTaxa", ignore_index=True)

    # Sum up the weights of a taxon and sort by weight
    print("Started summing the weights of a taxon and sorting them by weight...")
    unipept_frame["log_weight"] = np.log10(unipept_frame["weight"] + 1)
    tax_id_weights = unipept_frame.groupby("HigherTaxa")["log_weight"].sum().reset_index()

    # Since large proteomes tend to have more detectable peptides,
    # we adjust the weight by dividing by the size of the proteome i.e.,
    # the number of proteins that are associated with a taxon
    tax_id_weights["scaled_weight"] = tax_id_weights[
        "log_weight"
    ]  # / (TaxIDWeights["proteome_size"]) ** N

    # Retrieves the specified taxonomic rank taxid in the lineage of each of the species-level taxids returned by
    # Unipept for both the UnipeptFrame and the TaxIdWeightFrame
    higher_unique_psm_taxids = unique_psm_taxa  # set([GetLineageAtSpecifiedRank(i,TaxaRank) for i in UniquePSMTaxa])

    # group the duplicate entries of higher up taxa and sum their weights
    print("Started grouping duplicate entries of taxa situated higher up and sum their weights...")
    higher_taxid_weights = (
        tax_id_weights.groupby("HigherTaxa")["scaled_weight"]
        .sum()
        .reset_index()
        .sort_values(by=["scaled_weight"], ascending=False)
    )
    # HigherTaxidWeights = TaxIDWeights
    higher_taxid_weights["Unique"] = np.where(
        higher_taxid_weights["HigherTaxa"].isin(higher_unique_psm_taxids), True, False
    )

    try:
        higher_taxid_weights = higher_taxid_weights[
            higher_taxid_weights.HigherTaxa != 1869227
        ]
    except:
        pass

    if len(higher_taxid_weights.HigherTaxa) < 50:
        return unipept_frame, higher_taxid_weights
    else:
        taxa_to_include = set(higher_taxid_weights["HigherTaxa"][0:max_taxa])
        taxa_to_include.update(higher_unique_psm_taxids)
        return (
            unipept_frame[unipept_frame["HigherTaxa"].isin(taxa_to_include)],
            higher_taxid_weights,
        )
