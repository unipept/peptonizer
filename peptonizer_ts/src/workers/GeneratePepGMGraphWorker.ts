// Import Pyodide
import { loadPyodide, PyodideInterface } from 'pyodide';
import pythonCode from "./generate_pepgm_graph.py?raw";
import peptonizerWhlBase64 from "./lib/peptonizer_base64.whl?raw";

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
        indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.3/full/'
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
        import base64
        from pathlib import Path
        
        import micropip

        # Decode base64 string to binary and write to a temporary file
        wheel_data = "${peptonizerWhlBase64}"
        wheel_binary = base64.b64decode(wheel_data)
        
        # Define a temporary path for the .whl file
        wheel_path = Path("/tmp/peptonizer-0.1-py3-none-any.whl")
        wheel_path.write_bytes(wheel_binary)

        # Install the wheel package
        await micropip.install("emfs:///tmp/peptonizer-0.1-py3-none-any.whl")

        # Clean up by deleting the temporary file
        wheel_path.unlink()
    `);
}

let pyodideReadyPromise: Promise<void> = loadPyodideAndPackages();

self.onmessage = async (event: MessageEvent<EventData>): Promise<void> => {
    console.log("Received request in generate graph worker...");

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