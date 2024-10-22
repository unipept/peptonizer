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
                const {id, type, ...data} = JSON.parse(event.data);
                const [onSuccess, onError, progressListener] = this.callbacks.get(id);

                if (type === "result") {
                    onSuccess(JSON.parse(data.result));
                    this.callbacks.delete(id);
                } else if (type === "progress") {
                    if (data["progress_type"] === "graph") {
                        progressListener.graphsUpdated(
                            data["current_graph"],
                            data["total_graphs"]
                        );
                    } else if (data["progress_type"] === "max_residual") {
                        progressListener.maxResidualUpdated(
                            data["max_residual"],
                            data["tolerance"]
                        );
                    } else if (data["progress_type"] === "iterations") {
                        progressListener.iterationsUpdated(
                            data["current_iterations"],
                            data["total_iterations"]
                        )
                    }
                } else if (type === "error") {
                    onError(data.error);
                    this.callbacks.delete(id);
                }
            };

            this.workers.push(worker);
        }

        this.queue = async.queue((task, callback) => {
            (async () => {
                // Retrieve worker from the pool.
                const worker = this.workers.pop();

                if (worker) {
                    try {
                        const result = await new Promise((onSuccess, onError) => {
                            this.currentId = (this.currentId + 1) % Number.MAX_SAFE_INTEGER;
                            this.callbacks.set(
                                this.currentId,
                                [onSuccess, onError, task.progressListener]
                            );
                            task.progressListener.gridUpdated(task.alpha, task.beta, task.prior);
                            // The progressListener itself cannot be sent to the worker!
                            delete task.progressListener
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
     * @param progressListener
     * @returns {Promise}
     */
    static performGridSearch(graph, alphas, betas, priors, progressListener) {
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
                            prior,
                            progressListener
                        }, (result) => {
                            // Check if an error occurred. If so, stop running the grid search entirely.
                            if (result.error) {
                                reject(result.error);
                                return;
                            }

                            // No error, continue and add the results to the output
                            results.push(result)
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
