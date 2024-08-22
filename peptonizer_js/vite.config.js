import { defineConfig } from "vite";
import { copyFile, mkdir, readFile } from "fs/promises";
import { join } from "path";
import { createRequire } from 'module';

export default defineConfig({
    optimizeDeps: { exclude: ["pyodide"] },
    plugins: [
        {
            name: "vite-plugin-pyodide",
            config(config) {
                console.log('Vite is initializing with the following config:', config);
            },
            /**
             * This function is only called when building the library (not as part of the dev server!). This function
             * here handles adding the correct assets to the
             */
            generateBundle: async () => {
                console.log("Copying files...");
                const assetsDir = "dist/assets";
                await mkdir(assetsDir, { recursive: true });
                const files = [
                    "pyodide-lock.json",
                    "pyodide.asm.js",
                    "pyodide.asm.wasm",
                    "python_stdlib.zip",
                ];
                for (const file of files) {
                    await copyFile(
                        join("node_modules/pyodide", file),
                        join(assetsDir, file),
                    );
                }
            },
        },
    ],
});
