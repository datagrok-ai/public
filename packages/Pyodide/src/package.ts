/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { v4 as uuidv4 } from 'uuid';

export const _package = new DG.Package();

let pyodideWorker: Worker;
let convertDataFrame: (df: DG.DataFrame) => Promise<string | Uint8Array>;
let supportsArrow: boolean = false;
const currentExecutions: {[key: string]: Function} = {};

interface WorkerRequest {
  id: string;
  script: string;
  namespace: {[key: string]: any};
  outputs: string[];
}

interface WorkerResponse {
  id: string;
  error?: string;
  result?: {[key: string]: any};
}

//tags: init
export async function initPyodide() {
  pyodideWorker = new Worker(new URL('worker.js', import.meta.url));

  pyodideWorker.onmessage = (event) => {
    const result: WorkerResponse = event.data;
    const onSuccess = currentExecutions[result.id];
    delete currentExecutions[result.id];
    onSuccess(result);
  };
  const f: DG.Func[] = DG.Func.find({package: 'Arrow', name: 'toFeather'});
  supportsArrow = false /* f.length > 0 */;
  if (supportsArrow) {
    const toArrow: DG.Func = f[0];
    convertDataFrame = async (df: DG.DataFrame) => (await toArrow.apply({'table': df, 'asStream': true})) as Uint8Array;
  }
  else
    convertDataFrame = async (df: DG.DataFrame) => df.toCsv();
}


async function prepareRequest(scriptCall: DG.FuncCall): Promise<WorkerRequest> {
  let code = '';
  // refactor
  // @ts-ignore
  if (Object.values(scriptCall.inputParams).some((p: DG.FuncCallParam) => p.property.propertyType === DG.TYPE.DATA_FRAME)
      // @ts-ignore
      || Object.values(scriptCall.outputParams).some((p: DG.FuncCallParam) => p.property.propertyType === DG.TYPE.DATA_FRAME)) {
    if (supportsArrow)
      code += 'import pyarrow as pa\n';
    else {
      code += 'import pandas as pd\n';
      code += 'import io\n';
    }
  }
  const namespace: {[key: string]: any} = {};
  for (const paramName of Object.keys(scriptCall.inputParams)) {
    const value = scriptCall.inputs[paramName];
    const type: string = scriptCall.inputParams[paramName].property.propertyType;
    if (type === DG.TYPE.DATA_FRAME) {
      namespace[paramName] = await convertDataFrame(value);
      if (supportsArrow) {
        code += `${paramName} = pa.ipc.open_stream(${paramName})\n`;
        code += `${paramName} = ${paramName}.read_pandas()\n`;
      }
      else
        code += `${paramName} = pd.read_csv(io.StringIO(${paramName}), sep=",")\n`;
    }
    else if (type === DG.TYPE.FILE || type === DG.TYPE.BLOB) {
      if (value instanceof DG.FileInfo)
        namespace[paramName] = value.data;
      else
        namespace[paramName] = value;
    }
    else
      namespace[paramName] = value;
  }
  code += (scriptCall.func as DG.Script).clientCode;
  for (const paramName of Object.keys(scriptCall.outputParams)) {
    const type: string = scriptCall.outputParams[paramName].property.propertyType;
    if (type === DG.TYPE.DATA_FRAME) {
      if (supportsArrow) {
        code += `batch__${paramName} = pa.record_batch(${paramName})\n`;
        code += `sink__${paramName} = pa.BufferOutputStream()\n`;
        code += `with pa.ipc.new_stream(sink__${paramName}, batch__${paramName}.schema) as writer:\n`;
        code += `\twriter.write_batch(batch__${paramName})\n`;
        code += `${paramName} = sink__${paramName}.getvalue()\n`;
      }
      else
        code += `${paramName} = ${paramName}.to_csv(index=False)\n`;
    }
  }
  return {id: uuidv4(), script: code, namespace: namespace, outputs: Object.keys(scriptCall.outputParams)};
}

async function sendRequest(req: WorkerRequest): Promise<WorkerResponse> {
  return new Promise<WorkerResponse>((onSuccess) => {
    currentExecutions[req.id] = onSuccess;
    pyodideWorker.postMessage(req);
  });
}

//tags: scriptHandler
//meta.scriptHandler.language: pyodide
//meta.scriptHandler.extensions: py
//meta.scriptHandler.commentStart: #
//meta.scriptHandler.templateScript: #name: Template\n#description: Calculates number of cells in the table\n#language: pyodide\n#input: dataframe table [Data table]\n#output: int count [Number of cells in table]\n\ncount = table.shape[0] * table.shape[1]
//meta.scriptHandler.codeEditorMode: python
//meta.icon: files/pyodide.png
//input: funccall scriptCall
export async function pyodideLanguageHandler(scriptCall: DG.FuncCall): Promise<void> {
  const req: WorkerRequest = await prepareRequest(scriptCall);
  console.log(req.script);
  const response: WorkerResponse = await sendRequest(req); // spawn new worker if the current is busy?
  if (response.error)
    throw new Error(response.error);

  if (response.result) {
    for (const paramName of Object.keys(scriptCall.outputParams)) {

      const type: string = scriptCall.outputParams[paramName].property.propertyType;
      const value = response.result[paramName]!;
      if (value) {
        if (type == DG.TYPE.DATA_FRAME) {
          if (supportsArrow)
            scriptCall.setParamValue(paramName, await grok.functions.call('Arrow:fromFeather', {'bytes': value}));
          else
            scriptCall.setParamValue(paramName, DG.DataFrame.fromCsv(value));
        }
        else
          scriptCall.setParamValue(paramName, value);
      }
    }
  }
}