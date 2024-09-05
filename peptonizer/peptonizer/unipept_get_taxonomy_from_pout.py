import requests
import json
import logging

from typing import Dict, Tuple, List, Any

from .taxon_manager import TaxonManager

UNIPEPT_URL = "http://rick.ugent.be"
UNIPEPT_PEPT2FILTERED_ENDPOINT = "/mpa/pept2filtered.json"

UNIPEPT_PEPTIDES_BATCH_SIZE = 500


# TODO check if we should still take into account the cutoff parameter? Maybe this is no longer an issue because
# of the new Unipept API...
def query_unipept_and_filter_taxa(peptides: List[str], taxa_filter: List[int]) -> List[Any]:
    url = UNIPEPT_URL + UNIPEPT_PEPT2FILTERED_ENDPOINT
    filtered_peptides = []

    # Convert to set in order to efficiently perform set intersection later on
    taxa_filter = set(taxa_filter)

    print("Querying unipept and filtering taxa...")

    batches = len(peptides) // UNIPEPT_PEPTIDES_BATCH_SIZE

    # Split the peptides into batches of 100
    for i in range(0, len(peptides), UNIPEPT_PEPTIDES_BATCH_SIZE):
        print(f"Now querying Unipept with batch {(i // UNIPEPT_PEPTIDES_BATCH_SIZE) + 1} out of {batches}.")
        batch = peptides[i:i+UNIPEPT_PEPTIDES_BATCH_SIZE]

        # Prepare the request payload
        payload = {
            "peptides": batch
        }

        # Perform the HTTP POST request
        response = requests.post(url, json=payload)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()

            # Process each peptide result in the response
            for peptide_data in data.get("peptides", []):
                original_taxa = peptide_data.get("taxa", [])

                # Find the intersection of original taxa and taxa_filter
                filtered_taxa = list(set(original_taxa) & taxa_filter)

                # Replace the taxa with the filtered taxa
                peptide_data["taxa"] = filtered_taxa

                # Append the modified peptide data to the result list
                filtered_peptides.append(peptide_data)
        else:
            logging.error(f"Failed to retrieve peptide data for batch {batch}. Status code: {response.status_code}")

    return filtered_peptides


def fetch_unipept_taxon_information(
    pep_score_psm: Dict[str, Tuple[float, int]],
    unipept_peptide_counts_file: str,
    unipept_response_file: str,
    taxonomy_query: str,
    rank: str,
    log_file: str
):
    unipept_peptides = dict()

    for pep in pep_score_psm.keys():
        unipept_peptides[pep] = {
            "score": pep_score_psm[pep][0],
            "psms": pep_score_psm[pep][1],
        }

    with open(unipept_peptide_counts_file, "a") as f_out:
        f_out.write(json.dumps(unipept_peptides))

    # Set up the logger such that we can write errors to the provided file
    logging.basicConfig(filename=log_file, level=logging.INFO)

    taxa_filter = TaxonManager.get_descendants_for_taxa([int(item) for item in taxonomy_query.split(",")], rank)
    unipept_responses = query_unipept_and_filter_taxa(list(unipept_peptides.keys()), taxa_filter)

    with open(unipept_response_file, "w") as f_out:
        json.dump(unipept_responses, f_out)

