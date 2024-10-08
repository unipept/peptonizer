import { pepgmGraphGeneration } from "./py-worker-generation";
import { asyncPepgmExecution } from "./py-worker-execution";

/**
 * Start the peptonizer. This function takes in a PSM-file that has been read in earlier (so no file paths here). The
 * PSMS including their intensities are then used as input to the Peptonizer-pipeline. This pipeline will finally
 * return a Map in which NCBI taxon IDs are mapped onto their propabilities (as computed by the Peptonizer).
 *
 * @param psmContent Content of a PSM-file generated by upstream tools (such as MS2Rescore).
 * @return Mapping between NCBI taxon IDs (integer, > 0) and probabilities (float in [0, 1]).
 */
export async function peptonize(psmContent, pepGMParams) {
    return pepgmGraphGeneration({ psms: psmContent })
        .then(data => {
            if (data.error) return { error: data.error };
            const { alphas, betas, priors } = pepGMParams;
            return asyncPepgmExecution({ 
                graph: data.graph, 
                alphas: alphas, 
                betas: betas, 
                priors: priors});
        })
        .then(data => {
            if (data.error) return { error: data.error };
            return data.results;
        });
}

