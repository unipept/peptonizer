import async from "async";

const callbacks = {};
let id = 0; // identify a Promise

// Create numWorkers workers, which will execute pepGM for different parameters in parallel.
const numWorkers = 2
const executePepgmWorkers = []
for (let i = 0; i < numWorkers; i++) {
    let worker = new Worker("./workers/execute_pepgm_worker.js");
    // Define what to do if a pepGM execution is done.
    worker.onmessage = (event) => {
        const { id, ...data } = event.data;
        const onSuccess = callbacks[id];
        delete callbacks[id];
        onSuccess(data);
    };
    executePepgmWorkers.push(worker);
}

// Keep one pepgm result at each point of execution.
// Once a new result is ready, compare it to the current result and choose the best solution.
let result;
function aggregateResult(taskResult) {
    result = taskResult;
};

const queue = async.queue(async(task) => {
    // Retrieve worker from the pool.
    const worker = executePepgmWorkers.pop();

    if (worker) {
        const result = await new Promise((onSuccess) => {
            id = (id + 1) % Number.MAX_SAFE_INTEGER;
            callbacks[id] = onSuccess;
            worker.postMessage({ ...task, id });
        });

        // Add worker back to the pool.
        executePepgmWorkers.push(worker);

        // Aggregate the results
        aggregateResult(result);

        return result;
    } else {
        throw new Error("No workers available in the queue!");
    }
}, numWorkers);

// Execute pepGM with all given parameter combinations.
const asyncPepgmExecution = (() => {
    return (context) => {
        const { graph, alphas, betas, priors } = context;

        return new Promise((onSuccess) => {
            // Loop through all parameter combinations.
            alphas.forEach((alpha, i, as) => {
                betas.forEach((beta, j, bs) => {
                    priors.forEach((prior, k, ps) => {
                        // Push the tasks to the queue.
                        queue.push({
                            graph,
                            alpha,
                            beta,
                            prior
                        });
                    });
                });
            });

            // send pepGM onSucces if all tasks are processed
            queue.drain(function() {
                console.log(result);
                onSuccess(result);
            });
        });
    };
})();

export { asyncPepgmExecution };

