
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

export { pepgmGraphGeneration };