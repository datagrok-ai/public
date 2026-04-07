import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function initPyodide() : Promise<void> {
  await PackageFunctions.initPyodide();
}

//input: script script 
//output: string result
export function makeVectorCode(script: any) : string {
  return PackageFunctions.makeVectorCode(script);
}

//input: funccall scriptCall 
//meta.scriptHandler.language: pyodide
//meta.scriptHandler.extensions: py
//meta.scriptHandler.commentStart: #
//meta.scriptHandler.templateScript: #name: Template\n#description: Calculates number of cells in the table\n#language: pyodide\n#sample: cars.csv\n#input: dataframe table [Data table]\n#output: int count [Number of cells in table]\n\ncount = table.shape[0] * table.shape[1]
//meta.scriptHandler.codeEditorMode: python
//meta.scriptHandler.vectorizationFunction: Pyodide:makeVectorCode
//meta.icon: files/pyodide.png
//meta.role: scriptHandler
export async function pyodideLanguageHandler(scriptCall: DG.FuncCall) : Promise<void> {
  await PackageFunctions.pyodideLanguageHandler(scriptCall);
}
