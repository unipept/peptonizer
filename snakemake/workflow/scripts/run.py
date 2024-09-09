import json
from sys import getsizeof

from memory_profiler import profile

from peptonizer.peptonizer import fetch_unipept_taxon_information, parse_ms2rescore_output, perform_taxa_weighing, generate_pepgm_graph, run_belief_propagation, extract_taxon_scores

@profile
def run():
    # The PSM input should be provided to the parser as a list of strings
    psm_file = "/IdeaProjects/peptonizer/peptonizer_js/public/data/rescored_medium.psms.tsv"

    with open(psm_file, "r") as f:
        psms = f.read()

    print("Started parsing pout file from MS2Rescore...")
    parsed_input = parse_ms2rescore_output(psms, 0.01)
    print(f"Input has been parsed successfully...")

    print("Started fetching Unipept taxon information...")
    unipept_responses = fetch_unipept_taxon_information(
        parsed_input,
        "2",
        "species",
        "file_unipept_taxon_information_log"
    )
    print("Successfully fetched Unipept taxon information...")

    print("Started weighing taxa...")
    taxa_weights_df, _ = perform_taxa_weighing(
        unipept_responses,
        parsed_input,
        10,
        "species"
    )
    print("Successfully weighed taxa...")

    print(taxa_weights_df)

    print("Start creation of PepGM graph...")
    pepgm_graph = generate_pepgm_graph(taxa_weights_df)
    print("Successfully created PepGM graph...")

    print("Started running PepGM...")
    pepgm_results = run_belief_propagation(
        pepgm_graph,
        0.9,
        0.6,
        True,
        0.5
    )
    print("Successfully executed PepGM...")

    # Now convert the results from PepGM into a list of taxon IDs and the corresponding score values.
    final_scores = extract_taxon_scores(pepgm_results)

    print(json.dumps(final_scores))

run()
