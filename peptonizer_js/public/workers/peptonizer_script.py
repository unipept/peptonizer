import micropip
import json
import gc

from sys import getsizeof

await micropip.install("/lib/peptonizer-0.1-py3-none-any.whl")

import peptonizer

# The PSM input should be provided to the parser as a list of strings
psms = globals().get('input')

print("Started parsing pout file from MS2Rescore...")
parsed_input = peptonizer.parse_ms2rescore_output(psms, 0.01)
print(f"Input has been parsed successfully... --> size: {getsizeof(parsed_input)}")

print("Started fetching Unipept taxon information...")
unipept_responses = peptonizer.fetch_unipept_taxon_information(
    parsed_input,
    "2",
    "species",
    "file_unipept_taxon_information_log"
)
print("Successfully fetched Unipept taxon information...")

print("Started weighing taxa...")
taxa_weights_df, _ = peptonizer.perform_taxa_weighing(
    unipept_responses,
    parsed_input,
    10,
    "species"
)
print("Successfully weighed taxa...")

print("Start creation of PepGM graph...")
pepgm_graph = peptonizer.generate_pepgm_graph(taxa_weights_df)
print("Successfully created PepGM graph...")


def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


for name, size in sorted(((name, getsizeof(value)) for name, value in list(
        locals().items())), key=lambda x: -x[1])[:10]:
    print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))

print("Freeing up memory...")
del taxa_weights_df
del unipept_responses
del parsed_input
del input
del psms
del _
gc.collect()


def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


for name, size in sorted(((name, getsizeof(value)) for name, value in list(
        locals().items())), key=lambda x: -x[1])[:10]:
    print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))

print("Started running PepGM...")
pepgm_results = peptonizer.run_belief_propagation(
    pepgm_graph,
    0.9,
    0.6,
    True,
    0.5
)
print("Successfully executed PepGM...")

# Now convert the results from PepGM into a list of taxon IDs and the corresponding score values.
final_scores = peptonizer.extract_taxon_scores(pepgm_results)

json.dumps(final_scores)
