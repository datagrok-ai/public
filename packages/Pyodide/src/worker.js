// worker.js

importScripts("https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.js");
importScripts("https://cdnjs.cloudflare.com/ajax/libs/uuid/8.2.0/uuidv4.min.js");

let pyodide;
let currentCallsExecutions = {};

function valueToJS(value) {
    if (value instanceof pyodide.ffi.PyProxy)
        return value.toJs({dict_converter : Object.fromEntries, create_proxies : false});
    return value;
}

async function init() {
    pyodide = await loadPyodide();
    pyodide.setDebug(false);
    pyodide.registerJsModule("dg", {
        async execFuncCall(nqName, args) {
	    const jsArgs = {};
	    const typings = {};
	    const argsObject = valueToJS(args);
	    for (const [key, val] of Object.entries(argsObject)) {
		let jsVal = valueToJS(val);
		if (jsVal?.to_csv) {
		    jsVal = jsVal.to_csv.callKwargs({index: false});
		    jsVal = jsVal.trim();
		    typings[key] = "dataframe";
		} else if (jsVal?.isoformat) {
		    jsVal = jsVal.isoformat();
		    typings[key] = "datetime";
		}
		jsArgs[key] = jsVal;
	    }
	    const { promise, resolve, reject } = Promise.withResolvers();
	    const id = uuidv4();
	    currentCallsExecutions[id] = { resolve, reject };
	    postMessage({ id, type: 'startFuncCall', nqName, args: jsArgs, typings });
	    return promise;
        },
    });
}

let initCompleted = init();

async function handleWorkerScript(event) {
    const { id, script, namespace, outputs } = event.data;
    let scriptNamespace;
    let formatError;
    let reformatExceptionFunc;
    pyodide.setStdout({ batched: (msg) => postMessage({ id, type: 'stdoutMessage', message: msg }) });
    try {
        await pyodide.loadPackagesFromImports(script);
        scriptNamespace = pyodide.toPy(namespace);
        await pyodide.runPythonAsync(script, { globals: scriptNamespace });
        const result = {};
        for (const output of outputs) {
	    let value = scriptNamespace.get(output);
	    result[output] = valueToJS(value);
        }
        postMessage({ id, type: 'scriptResults', result });
    } catch (error) {
        reformatExceptionFunc = scriptNamespace.get("reformat_exception");
        formatError = reformatExceptionFunc ? reformatExceptionFunc() : [error?.message, error?.stack];
        error.message = formatError[0]?.trim();
        error.stack = formatError[1]?.trim();
        postMessage({ id, type: 'scriptResults', error });
    } finally {
        scriptNamespace?.destroy();
        formatError?.destroy();
        reformatExceptionFunc?.destroy();
        pyodide.setStdout();
    }
}

const convertDFScriptText = `
import pandas as pd
from io import StringIO
result = pd.read_csv(StringIO(arg), sep=",")
`;

const convertDateScriptText = `
from dateutil import parser
result = parser.parse(arg)
`;

async function handleFuncCallResults(event) {
    const data = event.data;
    const cbs = currentCallsExecutions[data.id];
    if (!cbs) {
	console.warn(`Pyodide worker unknown FuncCallResults id ${data.id}`);
	return;
    }
    delete currentCallsExecutions[data.id];
    if (data.error) {
        console.error(data.error);
	const pyError = pyodide.toPy(data.error);
	cbs.reject(pyError);
    } else {
	const {results} = data;
	// TODO: investigate do we need to destroy PyProxy here, since
	// it goes back pyhton
	for (const [name, typing] of Object.entries(data.typings ?? {})) {
	    if (typing === 'datetime') {
		const arg = data.results?.[name];
		results[name] = await runPyConvert(arg, convertDateScriptText);
	    } else if (typing === 'dataframe') {
	        const arg = data.results?.[name];
		results[name] = await runPyConvert(arg, convertDFScriptText);
	    }
	}
	const pyResults = pyodide.toPy(results);
	cbs.resolve(pyResults);
    }
}

async function runPyConvert(arg, script) {
    const globals = pyodide.toPy({arg});
    await pyodide.loadPackagesFromImports(script);
    await pyodide.runPythonAsync(script, { globals });
    return globals.get('result');
}

onmessage = async (event) => {
    await initCompleted;
    // console.log(event);
    switch (event.data.type) {
    case 'startWorkerScript':
	return await handleWorkerScript(event);
    case 'funcCallResults':
	return await handleFuncCallResults(event);
    }
}
