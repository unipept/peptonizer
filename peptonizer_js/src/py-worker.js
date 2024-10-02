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

const num_workers = 4
const executePepgmWorkers = []
for (let i = 0; i < num_workers; i++) {
    let worker = new Worker("./workers/execute_pepgm_worker.js");
    worker.onmessage = (event) => {
        const { id, ...data } = event.data;
        console.log(event.data);
        const onSuccess = callbacks[id];
        delete callbacks[id];
        onSuccess(data);
    };
    executePepgmWorkers.push(worker);
}

const alphas = [0.2, 0.5, 0.8];
const betas = [0.2, 0.5, 0.8];
const priors = [0.2, 0.5];
const asyncPepgmExecution = (() => {
    return (context) => {
        return new Promise((onSuccess) => {
            let workerId = 0;
            id = (id + 1) % Number.MAX_SAFE_INTEGER;
            callbacks[id] = onSuccess;
            executePepgmWorkers[workerId % num_workers].postMessage({
                ...context,
                id,
                alpha: 0.2,
                beta: 0.2,
                prior: 0.2
            });
            /*alphas.forEach((alpha, i, as) => {
                betas.forEach((beta, j, bs) => {
                    priors.forEach((prior, k, ps) => {
                        id = (id + 1) % Number.MAX_SAFE_INTEGER;
                        callbacks[id] = onSuccess;
                        executePepgmWorkers[workerId % num_workers].postMessage({
                            ...context,
                            id,
                            alpha,
                            beta,
                            prior
                        });
                        workerId ++;
                    });
                });
            });*/
        });
    };
})();

export { pepgmGraphGeneration, asyncPepgmExecution };