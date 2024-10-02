importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.2/full/pyodide.js");

async function loadPyodideAndPackages() {
    self.pyodide = await loadPyodide();
    // Load all packages into the Pyodide runtime environment that are required by the Peptonizer
    await self.pyodide.loadPackage([
        'numpy',
        'scipy',
        'networkx',
        'pandas',
        'micropip',
        'requests',
        'openssl'
    ]);
}
let pyodideReadyPromise = loadPyodideAndPackages();

self.onmessage = async (event) => {
    // make sure loading is done
    await pyodideReadyPromise;
    // Don't bother yet with this line, suppose our API is built in such a way:
    const { psms, id } = event.data;

    console.log(event.data);

    // Now is the easy part, the one that is similar to working in the main thread:
    try {
        pyodide.globals.set('input', psms);

        let results = await fetch('./generate_pepgm_graph.py')
                            .then(x => x.text())
                            .then(code => self.pyodide.runPythonAsync(code))
                            
        self.postMessage({ id: id, graph: results });
    } catch (error) {
        self.postMessage({ error: error.message, id: id });
    }
};