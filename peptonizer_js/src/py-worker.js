const generatePepgmGraphWorker = new Worker("./workers/generate_pepgm_graph_worker.js");

const callbacks = {};
let id = 0; // identify a Promise

// Called if graph generation is done.
generatePepgmGraphWorker.onmessage = (event) => {
    const { id, ...data } = event.data;
    const onSuccess = callbacks[id];
    delete callbacks[id];
    onSuccess(data);
};

// Start generating the graph.
const pepgmGraphGeneration = (() => {
    return (context) => {
        // the id could be generated more carefully
        id = (id + 1) % Number.MAX_SAFE_INTEGER;
        return new Promise((onSuccess) => {
            callbacks[id] = onSuccess;
            // Start a worker to generate the graph.
            generatePepgmGraphWorker.postMessage({
                ...context,
                id,
            });
        });
    };
})();

// Keep one pepgm result at each point of execution.
// Once a new result is ready, compare it to the current result and choose the best solution.
let result;
const aggregateResult = (taskResult) => {
    result = taskResult;
};

// Create numWorkers workers, which will execute pepGM for different parameters in parallel.
const numWorkers = 2
const executePepgmWorkers = []
let tasksDone = 0;
for (let i = 0; i < numWorkers; i++) {
    let worker = new Worker("./workers/execute_pepgm_worker.js");
    // Define what to do if a pepGM execution is done.
    worker.onmessage = (event) => {
        const { id, ...data } = event.data;
        const onSuccess = callbacks[id];
        delete callbacks[id];
        
        // Aggregate the results
        aggregateResult(data);

        // call onSuccess if all parameter combinations are done.
        tasksDone ++;
        if (tasksDone == tasksAmount) {
            onSuccess(result);
        }
    };
    executePepgmWorkers.push(worker);
}

// Execute pepGM with all given parameter combinations.
let tasksAmount = 0;
const asyncPepgmExecution = (() => {
    return (context) => {
        const { graph, alphas, betas, priors } = context;

        return new Promise((onSuccess) => {
            // Loop through all parameter combinations.
            alphas.forEach((alpha, i, as) => {
                betas.forEach((beta, j, bs) => {
                    priors.forEach((prior, k, ps) => {
                        id = (id + 1) % Number.MAX_SAFE_INTEGER;
                        callbacks[id] = onSuccess;
                        // Distribute the tasks among the workers.
                        executePepgmWorkers[tasksAmount % numWorkers].postMessage({
                            graph,
                            id,
                            alpha,
                            beta,
                            prior
                        });
                        // Count the amount of tasks started.
                        tasksAmount ++;
                    });
                });
            });
        });
    };
})();

export { pepgmGraphGeneration, asyncPepgmExecution };