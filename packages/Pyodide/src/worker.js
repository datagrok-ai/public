// worker.js

importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.2/full/pyodide.js");

let pyodide;

async function init() {
    pyodide = await loadPyodide();
    pyodide.setDebug(false);
    await pyodide.loadPackage(["numpy", "pandas", "matplotlib", /* "pyarrow" */]);
}

let initCompleted = init();

onmessage = async (event) => {
    await initCompleted;
    const { id, script, namespace, outputs } = event.data;
    let scriptNamespace;
    pyodide.setStdout({ batched: (msg) => postMessage({ id, message: msg }) });
    try {
        await pyodide.loadPackagesFromImports(script);
        scriptNamespace = pyodide.toPy(namespace);
        await pyodide.runPythonAsync(script, { globals: scriptNamespace });
        const result = {};
        for (const output of outputs) {
            let value = scriptNamespace.get(output);
            if (value instanceof pyodide.ffi.PyProxy)
                value = value.toJs({dict_converter : Object.fromEntries, create_proxies : false});
            result[output] = value;
        }
        postMessage({ id, result });
    } catch (error) {
        postMessage({ error: error.message, id });
    } finally {
        scriptNamespace?.destroy();
        pyodide.setStdout();
    }
}
