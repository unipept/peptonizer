// Import Pyodide
import { loadPyodide, PyodideInterface } from 'pyodide';
import pythonCode from "./execute_pepgm.py?raw";
// @ts-ignore
import peptonizerWhl from "./lib/peptonizer-0.1-py3-none-any.whl";


// Extend the DedicatedWorkerGlobalScope to include the pyodide property
interface DedicatedWorkerGlobalScope {
    pyodide: PyodideInterface;
}

interface EventData {
    graph: any;
    id: string;
    alpha: number;
    beta: number;
    prior: any;
}

declare const self: DedicatedWorkerGlobalScope & typeof globalThis;

async function loadPyodideAndPackages(): Promise<void> {
    self.pyodide = await loadPyodide({
        indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.2/full/'
    });
    // Load all packages into the Pyodide runtime environment that are required by the Peptonizer
    await self.pyodide.loadPackage([
        'micropip',
        'requests',
    ]);

    // Use the imported .whl file URL directly with micropip
    await self.pyodide.runPythonAsync(`
        import micropip
        await micropip.install("${peptonizerWhl}")
    `);
}

let pyodideReadyPromise: Promise<void> = loadPyodideAndPackages();

self.onmessage = async (event: MessageEvent<EventData>): Promise<void> => {
    // Ensure that loading is complete
    await pyodideReadyPromise;

    // Destructure the data from the event
    const { graph, id, alpha, beta, prior } = event.data;

    try {
        // Set input for Python code
        self.pyodide.globals.set('graph', graph);
        self.pyodide.globals.set('alpha', alpha);
        self.pyodide.globals.set('beta', beta);
        self.pyodide.globals.set('prior', prior);
        self.pyodide.globals.set('execution_id', id);

        // Fetch the Python code and execute it with Pyodide
        const results = await self.pyodide.runPythonAsync(pythonCode);

        // Post the result back to the main thread
        self.postMessage(JSON.stringify({ type: "result", result: results, id }));
    } catch (error: any) {
        // Post the error back to the main thread
        self.postMessage(JSON.stringify({ type: "error", error: error.message, id }));
    }
};
