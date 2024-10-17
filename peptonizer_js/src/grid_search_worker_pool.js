import async from "async";

class GridSearchWorkerPool {
    static callbacks = new Map();
    static currentId = 0;

    // Maximum amount of workers that are allowed to be run in parallel (using multiple CPU cores).
    static numberOfWorkers = 1;
    static workers = [];

    static queue;

    // Properly initialize the array of workers and their logic required to handle their asynchronous nature.
    static {
        for (let i = 0; i < this.numberOfWorkers; i++) {
            let worker = new Worker("./workers/execute_pepgm_worker.js");

            // Define what to do if a pepGM execution is done.
            worker.onmessage = (event) => {
                const { id, ...data } = event.data;
                const onSuccess = this.callbacks.get(id);
                this.callbacks.delete(id);
                onSuccess(data);
            };

            this.workers.push(worker);
        }

        this.queue = async.queue((task, callback) => {
            (async () => {
                // Retrieve worker from the pool.
                const worker = this.workers.pop();

                if (worker) {
                    try {
                        const result = await new Promise((onSuccess) => {
                            this.currentId = (this.currentId + 1) % Number.MAX_SAFE_INTEGER;
                            this.callbacks.set(this.currentId, onSuccess);
                            worker.postMessage({...task, id: this.currentId});
                        });

                        // Add worker back to the pool.
                        this.workers.push(worker);

                        // We're done, return the results to the callback
                        callback(result);
                    } catch (err) {
                        callback({ error: err })
                    }
                } else {
                    callback({ error: new Error("No workers available in the queue!") })
                }
            })();
        }, this.numberOfWorkers);
    }

    /**
     * Run the belief propagation algorithm on the provided graph for all possible combinations of the parameter
     * values.
     *
     * @param graph
     * @param alphas
     * @param betas
     * @param priors
     * @returns {Promise}
     */
    static performGridSearch(graph, alphas, betas, priors) {
        return new Promise((resolve, reject) => {
            const results = [];

            // Loop through all parameter combinations.
            alphas.forEach((alpha, i, as) => {
                betas.forEach((beta, j, bs) => {
                    priors.forEach((prior, k, ps) => {
                        // Push the tasks to the queue.
                        this.queue.push({
                            graph,
                            alpha,
                            beta,
                            prior
                        }, (result) => {
                            // Check if an error occurred. If so, stop running the grid search entirely.
                            if (result.error) {
                                reject(result.error);
                                return;
                            }

                            // No error, continue and add the results to the output
                            results.push(JSON.parse(result.results))
                        });
                    });
                });
            });

            // send pepGM onSucces if all tasks are processed
            this.queue.drain(() => {
                resolve(results);
            });
        });
    }
}

export { GridSearchWorkerPool };
