import { defineConfig } from "vite";
import { copyFile, mkdir } from "fs/promises";
import { join } from "path";

export default defineConfig({
    optimizeDeps: { exclude: ["pyodide"] },
    plugins: [
        {
            name: "vite-plugin-pyodide",
            /**
             * This function is only called when building the library (not as part of the dev server!). This function
             * here handles adding the correct assets to the final output.
             */
            generateBundle: async () => {
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
