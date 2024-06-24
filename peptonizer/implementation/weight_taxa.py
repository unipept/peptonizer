import argparse
import json
import numpy as np
import pandas as pd
from ete3 import NCBITaxa

ncbi = NCBITaxa()


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


def get_lineage_at_specified_rank(taxid, TaxaRank):
    """
    Returns the taxid of the specified rank in the lineage of any taxid given as argument
    -----
    taxid: int
         taxid to get the lineage of
    TaxaRank:
         rank you want to get
    """
    taxids = set()
    for tax in taxid:
        rank_taxid_dict = ncbi.get_rank(ncbi.get_lineage(tax))
        rank_taxid_dict = {rank: taxid for taxid, rank in rank_taxid_dict.items()}

        try:
            taxids.add(rank_taxid_dict[str(TaxaRank)])
        except:
            taxids.add(tax)

    return list(taxids)


def perform_taxa_weighing(
    unipept_response,
    pept_score_dict: str,
    max_tax,
    *peptides_per_taxon,
    chunks=True,
    n=0,
    select_rank=True,
    taxa_rank="species"
):
    """
    Weight inferred taxa based on their (1) degeneracy and (2) their proteome size.
    Parameters
    ----------
    unipept_response: str
        Path to Unipept response json file
    pept_score_dict: str
        Dictionary that contains peptide to score & number of PSMs map
    max_tax: int
        Maximum number of taxons to include in the graphical model
    peptides_per_taxon: str
        Path to the file that contains the size of the proteome per taxID (tab-separated)
    chunks: bool
        Allow memory-efficient reading of large json files
    n: int
        tbd

    Returns
    -------
    dataframe
        Top scoring taxa

    """
    print("Parsing Unipept responses from disk...")

    with open(pept_score_dict, "r") as file:
        pept_score_dict_loaded = json.load(file)

    if chunks:
        with open(unipept_response, 'r') as file:
            unipept_dict = {"peptides": []}
            for line in file:
                try:
                    # Get rid of the functional annotations from the response in order to speed up the JSON
                    # normalization further down in the script.
                    pept_data = json.loads(line)["peptides"]
                    for obj in pept_data:
                        del obj["fa"]
                    # Delete functional annotations (since we are not using these at the moment)
                    unipept_dict["peptides"].extend(pept_data)
                except:
                    # TODO: Pieter fixes internal server error
                    # in the meantime, we work with the incomplete mapping
                    # UnipeptDict["peptides"] = [json.loads(line)["peptides"]]
                    continue

    else:
        with open(unipept_response, "r") as file:
            unipept_dict = json.load(file)

    # Convert a JSON object into a Pandas DataFrame
    # record_path Parameter is used to specify the path to the nested list or dictionary that you want to normalize
    print("Normalizing peptides and converting to dataframe...")
    unipept_frame = pd.json_normalize(unipept_dict, record_path=["peptides"])
    # Merge psm_score and number of psms
    unipept_frame = pd.concat(
        [
            unipept_frame,
            pd.json_normalize(unipept_frame["sequence"].map(pept_score_dict_loaded)),
        ],
        axis=1,
    )

    # Score the degeneracy of a taxa, i.e.,
    # how conserved a peptide sequence is between taxa.
    # map all taxids in the list in the taxa column back to their taxid at species level
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
    # Retrieve the proteome size per taxid as a dictionary
    # This file was previously prepared by filtering a generic accession 2 taxid mapping file
    # to swissprot (i.e., reviewed) proteins only

    # Peptidome size: optional to include a weighting based on the size of the proteome, this didn't prove effective so
    # disabled for now.
    # PeptidomeSize: GetPeptideCountPerTaxID(PeptidesPerTaxon)
    # Map peptidome size and remove NAs
    # TaxIDWeights: TaxIDWeights[TaxIDWeights['taxa'].isin(PeptidomeSize.keys())].assign(proteome_size=lambda x: x['taxa'].map(PeptidomeSize))

    # Since large proteomes tend to have more detectable peptides,
    # we adjust the weight by dividing by the size of the proteome i.e.,
    # the number of proteins that are associated with a taxon
    tax_id_weights["scaled_weight"] = tax_id_weights[
        "log_weight"
    ]  # / (TaxIDWeights["proteome_size"]) ** N

    # Retrieves the specified taxonomic rank taxid in the lineage of each of the species-level taxids returned by
    # Unipept for both the UnipeptFrame and the TaxIdWeightFrame
    if select_rank:
        # TaxIDWeights['HigherTaxa'] = TaxIDWeights.apply(lambda row: GetLineageAtSpecifiedRank(row['taxa'],TaxaRank), axis = 1)
        # UnipeptFrame['HigherTaxa'] = UnipeptFrame.apply(lambda row: GetLineageAtSpecifiedRank(row['taxa'],TaxaRank), axis = 1)
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
        taxa_to_include = set(higher_taxid_weights["HigherTaxa"][0:max_tax])
        taxa_to_include.update(higher_unique_psm_taxids)
        return (
            unipept_frame[unipept_frame["HigherTaxa"].isin(taxa_to_include)],
            higher_taxid_weights,
        )
