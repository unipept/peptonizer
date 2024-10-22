// Import Pyodide
import { loadPyodide, PyodideInterface } from 'pyodide';
import pythonCode from "./generate_pepgm_graph.py?raw";
// @ts-ignore
import peptonizerWhl from "./lib/peptonizer-0.1-py3-none-any.whl";

// Extend the DedicatedWorkerGlobalScope to include the pyodide property
interface DedicatedWorkerGlobalScope {
    pyodide: PyodideInterface;
}

interface EventData {
    psms: any;
    id: string;
}

declare const self: DedicatedWorkerGlobalScope & typeof globalThis;

async function loadPyodideAndPackages(): Promise<void> {
    self.pyodide = await loadPyodide({
        indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.2/full/'
    });
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
    // Use the imported .whl file URL directly with micropip
    await self.pyodide.runPythonAsync(`
        import micropip
        await micropip.install("${peptonizerWhl}")
    `);
}

let pyodideReadyPromise: Promise<void> = loadPyodideAndPackages();

self.onmessage = async (event: MessageEvent<EventData>): Promise<void> => {
    // Make sure loading is done
    await pyodideReadyPromise;

    // Destructure the data from the event
    const { psms, id } = event.data;

    try {
        // Set inputs for the Python code
        self.pyodide.globals.set('input', psms);

        // Fetch the Python code and execute it with Pyodide
        const results = await self.pyodide.runPythonAsync(pythonCode);

        // Post the result back to the main thread
        self.postMessage({ id: id, graph: results });
    } catch (error: any) {
        // Post the error back to the main thread
        self.postMessage({ error: error.message, id: id });
    }
};