/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { v4 as uuidv4 } from 'uuid';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';

export const _package = new DG.Package();

let pyodideWorker: Worker;
let supportsArrow: boolean = false;
const currentExecutions: {[key: string]: Function} = {};
const currentCalls: {[key: string]: DG.FuncCall} = {};

interface StdoutMessage {
  id: string;
  message: string;
}

interface WorkerRequest {
  id: string;
  script: string;
  namespace: {[key: string]: any};
  outputs: string[];
}

interface WorkerResponse {
  id: string;
  error?: Error;
  result?: {[key: string]: any};
}

//tags: init
export async function initPyodide() {
  pyodideWorker = new Worker(new URL('worker.js', import.meta.url));

  pyodideWorker.onmessage = (event) => {
    if (Object.hasOwn(event.data, 'message')) {
      const result: StdoutMessage = event.data;
      currentCalls[result.id]?.debugLogger?.debug(result.message, {'IS_SERVICE_LOG': true});
    }
    else {
      const result: WorkerResponse = event.data;
      const onSuccess = currentExecutions[result.id];
      delete currentExecutions[result.id];
      delete currentCalls[result.id];
      onSuccess(result);
    }
  };
  supportsArrow = DG.Func.find({package: 'Arrow', name: 'toFeather'}).length > 0;
}

function makeCodeHeader(scriptCall: DG.FuncCall): string {
  let code: string = `def reformat_exception():
\tfrom traceback import format_exception_only, format_exception
\timport sys
\treturn ["".join(format_exception_only(sys.last_value)), "".join(format_exception(sys.last_value))]\n`;
  const inputParamsTypes: string[] = (Object.values(scriptCall.inputParams) as DG.FuncCallParam[]).map((p) =>  p.property.propertyType);
  const outputParamsTypes: string[] = (Object.values(scriptCall.outputParams) as DG.FuncCallParam[]).map((p) =>  p.property.propertyType);
  if (inputParamsTypes.some((type) => type === DG.TYPE.DATA_FRAME)
      || outputParamsTypes.some((type) => type === DG.TYPE.DATA_FRAME)) {
    if (supportsArrow)
      code += 'import pyarrow as pa\n';
    else {
      code += 'import pandas as pd\n';
      code += 'import io\n';
    }
  }
  if (outputParamsTypes.some((type) => type === DG.TYPE.GRAPHICS)) {
    code += 'import base64\n';
    code += 'import js\n';
    code += 'import io\n';
    code += 'import matplotlib\n';
    code += 'matplotlib.rcParams[\'figure.figsize\'] = (8.0, 8.0)\n';
    code += `class DomPatch__PlotLib:
\t\tdef __init__(self, *args, **kwargs):
\t\t\t\treturn
\t\tdef __getattr__(self, __name: str):
\t\t\t\treturn DomPatch__PlotLib
\t\tdef removeChild(self, *args, **kwargs):
\t\t\t\treturn\n`;
    code += 'js.document = DomPatch__PlotLib()\n';
  }
  if (inputParamsTypes.some((type) => type === DG.TYPE.DATE_TIME) || outputParamsTypes.some((type) => type === DG.TYPE.DATE_TIME)) {
    code += 'from datetime import datetime\n';
    code += 'from js import Date\n';
  }
  return code;
}

async function prepareInputs(scriptCall: DG.FuncCall, namespace: { [key: string]: any }): Promise<string> {
  let code: string = '';
  for (const paramName of Object.keys(scriptCall.inputParams)) {
    const value = scriptCall.inputs[paramName];
    const type: string = scriptCall.inputParams[paramName].property.propertyType;
    switch (type) {
      case DG.TYPE.DATA_FRAME:
        if (supportsArrow) {
          namespace[paramName] = await grok.functions.call('Arrow:toFeather', {'table': value, 'asStream': true});
          code += `${paramName} = pa.ipc.open_stream(${paramName})\n`;
          code += `${paramName} = ${paramName}.read_pandas()\n`;
        }
        else {
          namespace[paramName] = (value as DG.DataFrame)?.toCsv();
          code += `${paramName} = pd.read_csv(io.StringIO(${paramName}), sep=",")\n`;
        }
        break;
      case DG.TYPE.FILE:
      case DG.TYPE.GRAPHICS:
        throw new Error(`Parameter of type '${type}' is not supported as input`);
      case DG.TYPE.BLOB:
        namespace[paramName] = value instanceof DG.FileInfo ? await value?.readAsBytes() : value;
        break;
      case DG.TYPE.DATE_TIME:
        if (value)
          code += `${paramName} = datetime.fromtimestamp(${value.unix()})\n`;
        else
          code += `${paramName} = None\n`;
        break;
      case DG.TYPE.COLUMN:
        namespace[paramName] = (value as DG.Column)?.name;
        break;
      case DG.TYPE.COLUMN_LIST:
        namespace[paramName] = (value as DG.ColumnList)?.toList()?.map((c) => c.name);
        break;
      default:
        namespace[paramName] = value;
    }
  }
  return code;
}

function prepareOutputs(scriptCall: DG.FuncCall, outputs: string[]): string {
  let code: string = '';
  for (const paramName of Object.keys(scriptCall.outputParams)) {
    const type: string = scriptCall.outputParams[paramName].property.propertyType;
    switch (type) {
      case DG.TYPE.DATA_FRAME:
        if (supportsArrow) {
          code += `arrow_table__${paramName} = pa.Table.from_pandas(${paramName})\n`;
          code += `sink__${paramName} = pa.BufferOutputStream()\n`;
          code += `with pa.ipc.new_stream(sink__${paramName}, arrow_table__${paramName}.schema) as writer__${paramName}:\n`;
          code += `\twriter__${paramName}.write_table(arrow_table__${paramName})\n`;
          code += `${paramName} = sink__${paramName}.getvalue()\n`;
        }
        else
          code += `${paramName} = ${paramName}.to_csv(index=False)\n`;
        outputs.push(paramName);
        break;
      case DG.TYPE.DATE_TIME:
        code += `${paramName} = None if ${paramName} is None else Date.new(${paramName}.timestamp() * 1000)\n`;
        outputs.push(paramName);
        break;
      case DG.TYPE.GRAPHICS:
        code += `try:
\timport matplotlib.pyplot as plt__${paramName}__
\t${paramName}__io = io.BytesIO()    
\tplt__${paramName}__.savefig(${paramName}__io, format='png')
\t${paramName}__io.seek(0)
\t${paramName} = base64.b64encode(${paramName}__io.read()).decode()
except:
\t${paramName} = None
finally:
\tplt__${paramName}__.close()
\tif ${paramName}__io is not None:
\t\t${paramName}__io.close()\n`;
        outputs.push(paramName);
        break;
      case DG.TYPE.FILE:
      case DG.TYPE.COLUMN:
      case DG.TYPE.COLUMN_LIST:
        throw new Error(`Parameter of type '${type}' is not supported as output`);
      default:
        outputs.push(paramName);
    }
  }
  return code;
}

async function prepareRequest(scriptCall: DG.FuncCall): Promise<WorkerRequest> {
  let code = makeCodeHeader(scriptCall);
  const namespace: { [key: string]: any } = {};
  code += await prepareInputs(scriptCall, namespace);
  code += (scriptCall.func as DG.Script).clientCode;
  const outputs: string[] = [];
  code += prepareOutputs(scriptCall, outputs);
  return {id: uuidv4(), script: code, namespace: namespace, outputs: outputs};
}

async function sendRequest(req: WorkerRequest): Promise<WorkerResponse> {
  return new Promise<WorkerResponse>((onSuccess) => {
    currentExecutions[req.id] = onSuccess;
    pyodideWorker.postMessage(req);
  });
}

async function setOutputs(scriptCall: DG.FuncCall, response: WorkerResponse): Promise<void> {
  if (response.result) {
    for (const paramName of Object.keys(scriptCall.outputParams)) {
      const type: string = scriptCall.outputParams[paramName].property.propertyType;
      const value = response.result[paramName];
      if (value !== undefined && value !== null) {
        switch (type) {
          case DG.TYPE.DATA_FRAME:
            if (supportsArrow) {
              const df: DG.DataFrame = await grok.functions.call('Arrow:fromFeather', {'bytes': new Uint8Array(value.buffer)});
              df.name = paramName;
              scriptCall.setParamValue(paramName, df);
            }
            else
              scriptCall.setParamValue(paramName, DG.DataFrame.fromCsv(value.trim()));
            break;
          case DG.TYPE.FILE:
          case DG.TYPE.BLOB:
            scriptCall.setParamValue(paramName, DG.FileInfo.fromBytes(paramName, value));
            break;
          case DG.TYPE.DATE_TIME:
            scriptCall.setParamValue(paramName, dayjs.utc(value));
            break;
          case DG.TYPE.GRAPHICS:
            scriptCall.aux.set(paramName, 'image/png');
            scriptCall.setParamValue(paramName, value);
            break;
          default:
            scriptCall.setParamValue(paramName, value);
        }
      }
    }
  }
}

//name: makeVectorCode
//input: script script
//output: string code
export function makeVectorCode(script: DG.Script): string {
  const inputTable = DG.Script.vecInputTableName;
  let code: string = `input_table_length = len(${inputTable}.index)\n`;
  for (const param of script.inputs.filter((p) => p.isVectorizable))
    code += `${param.vectorName} = ${inputTable}["${param.name}"]\n`;
  for (const param of script.outputs.filter((p) => p.isVectorizable))
    code += `${param.vectorName} = [None] * input_table_length\n`;
  code += 'for vec_loop_idx in range(0, input_table_length):\n';

  for (const param of script.inputs.filter((p) => p.isVectorizable))
    code += `\t${param.name} = ${param.vectorName}[vec_loop_idx]\n`;
  code += `${script.clientCode.split('\n').map((l) => l.trim()).filter((l) => l.length > 0).map((l) => `\t${l}`).join('\n')}\n`;
  for (const param of script.outputs.filter((p) => p.isVectorizable))
    code += `\t${param.vectorName}[vec_loop_idx] = ${param.name}\n`;
  code += `${DG.Script.vecOutputTableName} = pd.DataFrame({${script.outputs.filter((p) => p.isVectorizable).map((p) => `"${p.vectorName}": ${p.vectorName}`).join(', ')}})\n`;
  return code;
}

//tags: scriptHandler
//meta.scriptHandler.language: pyodide
//meta.scriptHandler.extensions: py
//meta.scriptHandler.commentStart: #
//meta.scriptHandler.templateScript: #name: Template\n#description: Calculates number of cells in the table\n#language: pyodide\n#sample: cars.csv\n#input: dataframe table [Data table]\n#output: int count [Number of cells in table]\n\ncount = table.shape[0] * table.shape[1]
//meta.scriptHandler.codeEditorMode: python
//meta.scriptHandler.vectorizationFunction: Pyodide:makeVectorCode
//meta.icon: files/pyodide.png
//input: funccall scriptCall
export async function pyodideLanguageHandler(scriptCall: DG.FuncCall): Promise<void> {
  const req: WorkerRequest = await prepareRequest(scriptCall);
  currentCalls[req.id] = scriptCall;
  scriptCall.debugLogger?.debug(`Final code:\n${req.script}`);
  const response: WorkerResponse = await sendRequest(req); // spawn new worker if current is busy?
  if (response.error)
    throw response.error;
  await setOutputs(scriptCall, response);
}
