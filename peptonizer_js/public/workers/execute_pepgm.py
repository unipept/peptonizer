import json

print("Started running PepGM...")
pepgm_graph = globals().get('graph')
alpha = globals().get('alpha')
beta = globals().get('beta')
prior = globals().get('prior')

pepgm_results = peptonizer.run_belief_propagation(
    pepgm_graph,
    alpha,
    beta,
    True,
    prior
)
print("Successfully executed PepGM...")

# Now convert the results from PepGM into a list of taxon IDs and the corresponding score values.
final_scores = peptonizer.extract_taxon_scores(pepgm_results)

json.dumps(final_scores)