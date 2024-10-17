/**
 * A worker pool that can be used to generate factor graphs for PepGM that can, in turn, be send as input to the belief
 * propagation algorithm.
 */
class GraphGenerationWorkerPool {
    static worker = new Worker("./workers/generate_pepgm_graph_worker.js");
    static callbacks = new Map();
    static currentId = 0;

    // Properly initialize the worker and all logic required to handle its asynchronous nature.
    static {
        // Called by the worker when it's sending a message back to the main thread (on completion or on error)
        this.worker.onmessage = (event) => {
            const { id, ...data } = event.data;
            const onSuccess = this.callbacks.get(id);
            this.callbacks.delete(id);
            // Return the generated graphs to the caller
            onSuccess(data);
        }
    }

    /**
     * Generate the PepGM factor graph for a given list of PSM's. The computation of this factor graph will be queued
     * until the worker is available and the results will be returned as soon as they are available.
     *
     * @param psms A string containing the contents of a PSM-file that have already been parsed before.
     * @returns {Promise}
     */
    static generatePepGmGraph(psms) {
        // the id could be generated more carefully
        this.currentId = (this.currentId + 1) % Number.MAX_SAFE_INTEGER;
        return new Promise((onSuccess) => {
            // This promise will be resolved when the worker returns the results to the main thread asynchronously.
            this.callbacks.set(this.currentId, onSuccess);
            // Start a worker to generate the graph.
            this.worker.postMessage({
                psms,
                id: this.currentId,
            });
        });
    }
}


export { GraphGenerationWorkerPool };
