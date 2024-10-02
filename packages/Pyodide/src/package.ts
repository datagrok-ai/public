/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import {TYPE} from 'datagrok-api/dg';

export const _package = new DG.Package();

let pyodide: any;


//tags: init
export async function initPyodide() {
  // @ts-ignore
  pyodide = await loadPyodide();
}


function inputParameter(scriptCall: DG.FuncCall, inputName: string) {
  const value = scriptCall.inputs[inputName];
  const p = scriptCall.inputParams[inputName].property;
  switch (p.propertyType) {
    case TYPE.FLOAT:
    case TYPE.INT: return `${p.name} = ${value}`;
    case TYPE.STRING: return `${p.name} = "${value}"`;
    case TYPE.BOOL: return `${p.name} = ${value ? 'True' : 'False'}`;
  }
}


//input: funccall scriptCall
//meta.languageHandler: pyodide
export async function pyodideLanguageHandler(scriptCall: DG.FuncCall): Promise<void> {
  // for (const paramName of Object.keys(scriptCall.outputParams))
  //   scriptCall.setParamValue(paramName, 42);
  // return;

  let s = '';
  for (const paramName of Object.keys(scriptCall.inputParams)) {
    s += inputParameter(scriptCall, paramName) + '\n';
  }

  s += (scriptCall.func as DG.Script).script + '\n\n';

  console.log(s);
  await pyodide.runPythonAsync(s);

  for (const paramName of Object.keys(scriptCall.outputParams)) {
    console.log(`${paramName}: ${pyodide.globals.get(paramName)}`);
    scriptCall.setParamValue(paramName, pyodide.globals.get(paramName));
  }
}


//meta.languageHandler: pyodide
//input: string pythonScript
//output: dynamic result
export async function runPython(pythonScript: string): Promise<any> {
  return pyodide.runPythonAsync(pythonScript);
}