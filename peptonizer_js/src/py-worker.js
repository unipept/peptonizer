const generatePepgmGraphWorker = new Worker("./workers/generate_pepgm_graph_worker.js");

const callbacks = {};
let id = 0; // identify a Promise

generatePepgmGraphWorker.onmessage = (event) => {
    const { id, ...data } = event.data;
    const onSuccess = callbacks[id];
    delete callbacks[id];
    onSuccess(data);
};

const pepgmGraphGeneration = (() => {
    return (context) => {
        // the id could be generated more carefully
        id = (id + 1) % Number.MAX_SAFE_INTEGER;
        return new Promise((onSuccess) => {
            callbacks[id] = onSuccess;
            generatePepgmGraphWorker.postMessage({
                ...context,
                id,
            });
        });
    };
})();

let result;
const aggregateResult = (taskResult) => {
    result = taskResult;
};

const numWorkers = 2
const executePepgmWorkers = []
let tasksDone = 0;
for (let i = 0; i < numWorkers; i++) {
    let worker = new Worker("./workers/execute_pepgm_worker.js");
    worker.onmessage = (event) => {
        const { id, ...data } = event.data;
        const onSuccess = callbacks[id];
        delete callbacks[id];
        
        aggregateResult(data);
        tasksDone ++;
        if (tasksDone == tasksAmount) {
            onSuccess(result);
        }
    };
    executePepgmWorkers.push(worker);
}

let tasksAmount = 0;
const asyncPepgmExecution = (() => {
    return (context) => {
        const { graph, alphas, betas, priors } = context;

        return new Promise((onSuccess) => {
            alphas.forEach((alpha, i, as) => {
                betas.forEach((beta, j, bs) => {
                    priors.forEach((prior, k, ps) => {
                        id = (id + 1) % Number.MAX_SAFE_INTEGER;
                        callbacks[id] = onSuccess;
                        executePepgmWorkers[tasksAmount % numWorkers].postMessage({
                            graph,
                            id,
                            alpha,
                            beta,
                            prior
                        });
                        tasksAmount ++;
                    });
                });
            });
        });
    };
})();

export { pepgmGraphGeneration, asyncPepgmExecution };