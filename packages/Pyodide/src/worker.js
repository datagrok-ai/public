// worker.js

importScripts("https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.js");
importScripts("https://cdnjs.cloudflare.com/ajax/libs/uuid/8.2.0/uuidv4.min.js");

let pyodide;
let currentCallsExecutions = {};

function eager_converter(item, fn) {
    const isInt64 = item?.type === 'numpy.int64';
    const res = fn(item);
    if (isInt64 && res.length === 1) {
        return Number(res[0]);
    }
    return res;
}

function valueToJS(value) {
    if (value instanceof pyodide.ffi.PyProxy)
        return value.toJs({dict_converter : Object.fromEntries, create_proxies : false, eager_converter});
    return value;
}

async function init() {
    pyodide = await loadPyodide();
    pyodide.setDebug(false);
    pyodide.registerJsModule('dg', {
        async execFuncCall(nqName, args) {
            const jsArgs = {};
            const typings = {};
            const argsObject = valueToJS(args);
            for (const [key, val] of Object.entries(argsObject)) {
                let jsVal = val;
                if (jsVal?.to_csv) {
                    jsVal = jsVal.to_csv.callKwargs({index: false});
                    jsVal = jsVal.trim();
                    typings[key] = 'dataframe';
                } else if (jsVal?.isoformat) {
                    jsVal = jsVal.isoformat();
                    typings[key] = 'datetime';
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

const loadedItems = new Set();
const execScriptName = '<exec>';
const pipName = 'micropip';

async function handleWorkerScript(event) {
    const { id, script, scriptName, headerLinesCount, namespace, outputs, dependencies } = event.data;
    let scriptNamespace;
    let formatError;
    let reformatExceptionFunc;
    pyodide.setStdout({ batched: (msg) => postMessage({ id, type: 'stdoutMessage', message: msg }) });
    try {
        if (dependencies?.length) {
            if (!loadedItems.has(pipName)) {
                await pyodide.loadPackage(pipName);
                loadedItems.add(pipName);
            }
            const micropip = pyodide.pyimport("micropip");
            for (const dep of dependencies) {
                if (loadedItems.has(dep))
                    continue;
                await micropip.install(dep);
                loadedItems.add(dep);
            }
        }
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
        reformatExceptionFunc = scriptNamespace?.get('reformat_exception');
        formatError = reformatExceptionFunc ? reformatExceptionFunc() : [error?.message, error?.stack, null];
        error.message = formatError[0]?.trim();
        error.stack = formatError[1]?.trim();
        if (formatError[2] != null) {
            const trace = valueToJS(formatError[2]).map(frameSummaryProxy => {
                const item = { filename: frameSummaryProxy.filename, lineno: frameSummaryProxy.lineno };
                frameSummaryProxy.destroy();
                return item;
            });
            const traceStart = trace.findIndex(item => item.filename === execScriptName);
            if (traceStart >= 0) {
                const actualTrace = trace.slice(traceStart);
                const strNum = actualTrace[0].lineno;
                const traceTail = actualTrace.slice(1);
                const strIdx = strNum - 1;
                const sourceStrNum = strNum - headerLinesCount;
                const linesArray = script.split('\n');
                linesArray[strIdx] = '> ' + linesArray[strIdx];
                const shownLines = linesArray.slice(Math.max(strIdx - 2, 0), strIdx + 3);
                const text = [
                    `Pyodide Error: ${error.message}`,
                    `Script ${scriptName}, Line ${sourceStrNum}`,
                    ...traceTail.map(({ filename, lineno }) => execScriptName === filename ?
                       `Script ${scriptName}, Line ${lineno - headerLinesCount}` :
                       `File ${filename}, Line ${lineno - (filename === '')}`),
                    ...shownLines
                ].join('\n');
                console.error(text);
            }
        }
        postMessage({ id, type: 'scriptResults', error });
    } finally {
        scriptNamespace?.destroy?.();
        formatError?.destroy?.();
        reformatExceptionFunc?.destroy?.();
        pyodide.setStdout();
    }
}

const convertScripts = {
    convertDFScriptText: `
import pandas as pd
from io import StringIO
result = pd.read_csv(StringIO(arg), sep=",")
`,
    convertDateScriptText: `
from dateutil import parser
result = parser.parse(arg)
`
}

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
        for (const [name, typing] of Object.entries(data.typings ?? {})) {
            if (typing === 'datetime') {
                const arg = data.results?.[name];
                results[name] = await runPyConvert(arg, 'convertDateScriptText');
            } else if (typing === 'dataframe') {
                const arg = data.results?.[name];
                results[name] = await runPyConvert(arg, 'convertDFScriptText');
            }
        }
        const pyResults = pyodide.toPy(results);
        cbs.resolve(pyResults);
    }
}

const loadedSystemScripts = new Set();

async function runPyConvert(arg, id) {
    const globals = pyodide.toPy({arg});
    const script = convertScripts[id];
    if (!loadedSystemScripts.has(id)) {
        await pyodide.loadPackagesFromImports(script);
        loadedSystemScripts.add(id);
    }
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
