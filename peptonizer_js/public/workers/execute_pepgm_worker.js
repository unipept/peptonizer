importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.2/full/pyodide.js");

async function loadPyodideAndPackages() {
    self.pyodide = await loadPyodide();
    // Load all packages into the Pyodide runtime environment that are required by the Peptonizer
    await self.pyodide.loadPackage([
        'micropip',
        'requests',
    ]);
}
let pyodideReadyPromise = loadPyodideAndPackages();
    

self.onmessage = async (event) => {
    // make sure loading is done
    await pyodideReadyPromise;
    // Don't bother yet with this line, suppose our API is built in such a way:
    const { graph, id, alpha, beta, prior } = event.data;

    console.log(event.data);

    // Now is the easy part, the one that is similar to working in the main thread:
    try {
        // Set input for python code.
        pyodide.globals.set('graph', graph);
        pyodide.globals.set('alpha', alpha);
        pyodide.globals.set('beta', beta);
        pyodide.globals.set('prior', prior);
        pyodide.globals.set('execution_id', id);

        // Fetch the python code and execute it with Pyodide.
        let results = await fetch('./execute_pepgm.py')
                            .then(x => x.text())
                            .then(code => self.pyodide.runPythonAsync(code));

        self.postMessage(JSON.stringify({ type: "result", result: results, id }));
    } catch (error) {
        self.postMessage(JSON.stringify({ type: "error", error: error.message, id }));
    }
};