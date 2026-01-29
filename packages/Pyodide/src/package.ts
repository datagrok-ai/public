/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { v4 as uuidv4 } from 'uuid';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
export * from './package.g';
export const _package = new DG.Package();

let pyodideWorker: Worker;
let supportsArrow: boolean = false;
const currentExecutions: {[key: string]: {resolve: Function, reject: Function}} = {};
const currentScriptCalls: {[key: string]: DG.FuncCall} = {};

interface StdoutMessage {
  id: string;
  type: 'stdoutMessage';
  message: string;
}

interface StartWorkerScript {
  id: string;
  type: 'startWorkerScript';
  script: string;
  scriptName: string;
  headerLinesCount: number;
  namespace: {[key: string]: any};
  dependencies?: string[];
  outputs: string[];
}

interface ScriptResults {
  id: string;
  type: 'scriptResults';
  error?: Error;
  result?: {[key: string]: any};
}

interface StartFuncCall {
  id: string;
  type: 'startFuncCall';
  nqName: string;
  args?: {[key: string]: any};
  typings?: {[key: string]: any};
}

interface FuncCallResults {
  id: string;
  type: 'funcCallResults';
  error?: Error;
  results?: {[key: string]: any};
  typings?: {[key: string]: any};
}

function handleStdoutMessage(result: StdoutMessage) {
  console.log(result.message);
  currentScriptCalls[result.id]?.debugLogger?.debug(result.message, {'IS_SERVICE_LOG': true});
}

async function handleScriptResults(result: ScriptResults) {
  const cbs = currentExecutions[result.id];
  delete currentExecutions[result.id];
  delete currentScriptCalls[result.id];
  if (!cbs) {
    console.warn(`Pyodide worker unknown script id ${result.id}`);
    return;
  }
  if (result.error)
    cbs.reject(result.error);
  else
    cbs.resolve(result.result);
}

async function handleStartFuncCall(data: StartFuncCall) {
  const func = DG.Func.byName(data.nqName);
  const fc = func.prepare();

  for (const paramName of Object.keys(fc.inputParams)) {
    let inputVal = data.args?.[paramName];
    if (data.typings?.[paramName]) {
      if (data.typings[paramName] === 'datetime')
        inputVal = dayjs(inputVal)
      if (data.typings[paramName] === 'dataframe')
        inputVal = DG.DataFrame.fromCsv(inputVal)
    }
    fc.inputs[paramName] = inputVal;
  }

  const jsRes: Record<string, any> = {};
  const typings: Record<string, any> = {};
  let error: any = undefined;

  try {
    await fc.call();

    for (const paramName of Object.keys(fc.outputParams)) {
      let outputVal = fc.outputs[paramName];
      const param  = fc.outputParams[paramName];
      const type = param.property.propertyType;

      if (type === DG.TYPE.DATA_FRAME) {
        outputVal = outputVal.toCsv();
        typings[paramName] = DG.TYPE.DATA_FRAME;
      } else if (type === DG.TYPE.DATE_TIME) {
        outputVal = outputVal.format();
        typings[paramName] = DG.TYPE.DATE_TIME;
      }
      jsRes[paramName] = outputVal;
    }
  } catch (e) {
    error = e;
  }

  const payload: FuncCallResults = {
    id: data.id,
    type: 'funcCallResults',
    error,
    results: jsRes,
    typings,
  };
  pyodideWorker.postMessage(payload);
}

function makeCodeHeader(scriptCall: DG.FuncCall): string {
  let code: string = `def reformat_exception():
\tfrom traceback import format_exception_only, format_exception, extract_tb
\timport sys
\treturn ["".join(format_exception_only(sys.last_value)), "".join(format_exception(sys.last_value)), extract_tb(sys.last_value.__traceback__)]\n`;
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
    code += '\n';
    const type: string = scriptCall.outputParams[paramName].property.propertyType;
    switch (type) {
      case DG.TYPE.DATA_FRAME:
        if (supportsArrow) {
          code += `batch__${paramName} = pa.record_batch(${paramName})\n`;
          code += `sink__${paramName} = pa.BufferOutputStream()\n`;
          code += `with pa.ipc.new_stream(sink__${paramName}, batch__${paramName}.schema) as writer:\n`;
          code += `\twriter.write_batch(batch__${paramName})\n`;
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

async function prepareRequest(scriptCall: DG.FuncCall): Promise<StartWorkerScript> {
  let code = makeCodeHeader(scriptCall);
  const func = (scriptCall.func as DG.Script);
  const scriptName = func.nqName;
  const namespace: { [key: string]: any } = {};
  code += await prepareInputs(scriptCall, namespace);
  const headerLinesCount = Math.max(0, code.split('\n').length - 1);
  code += func.clientCode;
  const outputs: string[] = [];
  code += prepareOutputs(scriptCall, outputs);
  let dependencies;
  if (scriptCall.func.options['dependencies'])
    dependencies = JSON.parse(scriptCall.func.options['dependencies']);
  return {id: uuidv4(), type: 'startWorkerScript', script: code, headerLinesCount, scriptName, namespace, dependencies, outputs};
}

async function sendRequest(req: StartWorkerScript): Promise<Record<string, any>> {
  // @ts-ignore-next-line
  const { promise, resolve, reject } = Promise.withResolvers();
  currentExecutions[req.id] = { resolve, reject };
  pyodideWorker.postMessage(req);
  return promise;
}

async function setOutputs(scriptCall: DG.FuncCall, results: Record<string, any>): Promise<void> {
  for (const paramName of Object.keys(scriptCall.outputParams)) {
    const type: string = scriptCall.outputParams[paramName].property.propertyType;
    const value = results[paramName];
    if (value !== undefined && value !== null) {
      switch (type) {
        case DG.TYPE.DATA_FRAME:
          if (supportsArrow)
            scriptCall.setParamValue(paramName, await grok.functions.call('Arrow:fromFeather', {'bytes': value}));
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


export class PackageFunctions{
  @grok.decorators.init()
  static async initPyodide() : Promise<void> {
  
    pyodideWorker = new Worker(new URL('worker.js', import.meta.url));

    pyodideWorker.onmessage = async (event: MessageEvent<StdoutMessage | ScriptResults | StartFuncCall>) => {
      // console.log(event);
      switch (event.data.type) {
        case 'stdoutMessage':
          return handleStdoutMessage(event.data);
        case 'scriptResults':
          return await handleScriptResults(event.data);
        case 'startFuncCall':
          return await handleStartFuncCall(event.data);
      }
    };
    /* supportsArrow = DG.Func.find({package: 'Arrow', name: 'toFeather'}).length > 0; */
  }

  @grok.decorators.func()
  static makeVectorCode(
    @grok.decorators.param({type: 'script'}) script: DG.Script): string {
  
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


  @grok.decorators.func({
    'meta': {
      'scriptHandler.language': 'pyodide',
      'scriptHandler.extensions': 'py',
      'scriptHandler.commentStart': '#',
      'scriptHandler.templateScript': '#name: Template\\n#description: Calculates number of cells in the table\\n#language: pyodide\\n#sample: cars.csv\\n#input: dataframe table [Data table]\\n#output: int count [Number of cells in table]\\n\\ncount = table.shape[0] * table.shape[1]',
      'scriptHandler.codeEditorMode': 'python',
      'scriptHandler.vectorizationFunction': 'Pyodide:makeVectorCode',
      'icon': 'files/pyodide.png',
      'role': 'scriptHandler',
    },
  })
  static async pyodideLanguageHandler(
    scriptCall: DG.FuncCall): Promise<void> {
  
    const req: StartWorkerScript = await prepareRequest(scriptCall);
    currentScriptCalls[req.id] = scriptCall;
    scriptCall.debugLogger?.debug(`Final code:\n${req.script}`);
    const response = await sendRequest(req); // spawn new worker if current is busy?
    await setOutputs(scriptCall, response);
  }
}
