import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {IKnimeClient} from './knime-client';
import {KnimeInputParam, KnimeOutputParam, KnimeParamType, KnimeDeployment, KnimeExecutionResult} from './types';
import {dataFrameToKnimeTable, knimeTableToDataFrame, knimeSpecDataToDataFrame} from './data-conversion';
import {pollJobUntilComplete} from './utils';

interface ParamMeta {
  sanitizedName: string;
  originalName: string;
  type: KnimeParamType;
  group?: string;
  groupDescription?: string;
  description?: string;
  defaultValue?: any;
}

const knimeTypeToDgType: Record<KnimeParamType, string> = {
  'table': 'dataframe',
  'file': 'file',
  'string': 'string',
  'int': 'int',
  'double': 'double',
  'boolean': 'bool',
  'json': 'string',
};

export function sanitizeFuncName(name: string, id: string): string {
  let result = name.replace(/[^a-zA-Z0-9_]/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '');
  if (!result)
    result = id.replace(/[^a-zA-Z0-9_]/g, '_').replace(/_+/g, '_').substring(0, 20);
  return result;
}

export function sanitizeParamName(name: string): string {
  let result = name.replace(/[^a-zA-Z0-9_]/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '');
  if (!result || /^\d/.test(result))
    result = 'p_' + result;
  return result;
}

function buildParamMeta(params: KnimeInputParam[]): ParamMeta[] {
  const meta: ParamMeta[] = [];
  const usedNames = new Set<string>();

  for (const param of params) {
    let baseName = param.group
      ? sanitizeParamName(param.group) + '_' + sanitizeParamName(param.name)
      : sanitizeParamName(param.name);
    let name = baseName;
    let n = 2;
    while (usedNames.has(name)) {
      name = `${baseName}_${n}`;
      n++;
    }
    usedNames.add(name);

    meta.push({
      sanitizedName: name,
      originalName: param.name,
      type: param.type,
      group: param.group,
      groupDescription: param.groupDescription,
      description: param.description,
      defaultValue: param.defaultValue,
    });
  }

  return meta;
}

function determineReturnType(outputs: KnimeOutputParam[]): string {
  if (outputs.length === 0)
    return 'object';
  if (outputs.length === 1)
    return knimeTypeToDgType[outputs[0].type] ?? 'object';
  const tableOutput = outputs.find((o) => o.type === 'table');
  if (tableOutput)
    return 'dataframe';
  return knimeTypeToDgType[outputs[0].type] ?? 'object';
}

function buildSignature(funcName: string, paramMeta: ParamMeta[], outputs: KnimeOutputParam[]): string {
  const parts: string[] = [];
  for (const pm of paramMeta) {
    const dgType = knimeTypeToDgType[pm.type] ?? 'string';
    parts.push(`${dgType} ${pm.sanitizedName}`);
  }
  const returnType = determineReturnType(outputs);
  return `${returnType} ${funcName}(${parts.join(', ')})`;
}

function createRunCallback(
  deployment: KnimeDeployment,
  paramMeta: ParamMeta[],
  outputs: KnimeOutputParam[],
  client: IKnimeClient,
): Function {
  const returnType = determineReturnType(outputs);

  return async (...args: any[]) => {
    const inputs: {[key: string]: any} = {};

    for (let i = 0; i < paramMeta.length; i++) {
      const pm = paramMeta[i];
      const val = args[i];
      if (val === null || val === undefined)
        continue;

      let converted: any;
      if (pm.type === 'table')
        converted = dataFrameToKnimeTable(val as DG.DataFrame);
      else if (pm.type === 'file') {
        const fileInfo = val as DG.FileInfo;
        const bytes = await fileInfo.readAsBytes();
        converted = new File([bytes.buffer as ArrayBuffer], fileInfo.name);
      }
      else if (pm.type === 'json') {
        try { converted = val ? JSON.parse(val as string) : null; }
        catch { converted = val; }
      }
      else
        converted = val;

      if (pm.group) {
        if (!inputs[pm.group])
          inputs[pm.group] = {};
        inputs[pm.group][pm.originalName] = converted;
      }
      else
        inputs[pm.originalName] = converted;
    }

    const isRestService = deployment.id.startsWith('rest:');
    let result: KnimeExecutionResult;

    if (isRestService)
      result = await client.executeSyncWorkflow(deployment.id, inputs);
    else {
      const jobId = await client.startAsyncJob(deployment.id, inputs);
      result = await pollJobUntilComplete(client, jobId);
    }

    return extractReturnValue(result, returnType);
  };
}

function extractReturnValue(result: KnimeExecutionResult, returnType: string): any {
  const tables: DG.DataFrame[] = [];
  const variables: {[key: string]: any} = {};

  if (result.outputTables)
    for (const t of result.outputTables) {
      const df = t?.['table-spec'] && t?.['table-data']
        ? knimeSpecDataToDataFrame(t['table-spec'], t['table-data'])
        : Array.isArray(t) ? knimeTableToDataFrame(t) : knimeTableToDataFrame([t]);
      tables.push(df);
    }

  if (result.fetchedResources)
    for (const {name: fileName, resource} of result.fetchedResources) {
      if (resource.json !== undefined && Array.isArray(resource.json)) {
        const df = knimeTableToDataFrame(resource.json, fileName);
        df.name = fileName;
        tables.push(df);
      }
      else if (resource.text !== undefined) {
        const isCsv = resource.contentType.includes('csv') || fileName.endsWith('.csv');
        const isTsv = resource.contentType.includes('tab-separated') || fileName.endsWith('.tsv');
        if (isCsv || isTsv) {
          try {
            const df = DG.DataFrame.fromCsv(resource.text);
            df.name = fileName;
            tables.push(df);
          }
          catch { /* not parseable */ }
        }
      }
    }

  if (returnType === 'dataframe')
    return tables.length > 0 ? tables[0] : DG.DataFrame.create();

  if (returnType === 'object') {
    if (tables.length > 0)
      return tables[0];
    if (Object.keys(variables).length > 0)
      return variables;
    return result.outputs;
  }

  // Scalar types — return first variable value
  const vals = Object.values(variables);
  if (vals.length > 0)
    return vals[0];
  return null;
}

export async function getOrRegisterFunc(
  deployment: KnimeDeployment,
  client: IKnimeClient,
): Promise<DG.Func> {
  const funcName = sanitizeFuncName(deployment.name, deployment.id);

  const existing = DG.Func.find({package: 'KnimeLink', name: funcName});
  if (existing.length > 0)
    return existing[0];

  const spec = await client.getWorkflowInputs(deployment.id);
  const paramMeta = buildParamMeta(spec.inputs);
  const signature = buildSignature(funcName, paramMeta, spec.outputs);
  const runCallback = createRunCallback(deployment, paramMeta, spec.outputs, client);

  const func = grok.functions.register({
    signature,
    run: runCallback,
    isAsync: true,
    namespace: 'KnimeLink',
  });

  console.log('^^^^^^^^^^^^^^ function registered');
  console.log(func.nqName);
  //console.log(func.package);

  // Set descriptions and default values on input properties
  for (const prop of func.inputs) {
    const pm = paramMeta.find((m) => m.sanitizedName === prop.name);
    if (pm) {
      if (pm.group) {
        const parts: string[] = [];
        if (pm.description)
          parts.push(pm.description);
        const groupLabel = pm.groupDescription
          ? `Part of "${pm.group}" parameter (${pm.groupDescription})`
          : `Part of "${pm.group}" parameter`;
        parts.push(groupLabel);
        prop.description = parts.join('. ');
      }
      else if (pm.description)
        prop.description = pm.description;
      if (pm.defaultValue !== undefined && pm.defaultValue !== null) {
        if (pm.type === 'string' || pm.type === 'json') {
          const str = pm.type === 'json' ? JSON.stringify(pm.defaultValue, null, 2) : String(pm.defaultValue);
          prop.initialValue = `"${str.replace(/\\/g, '\\\\').replace(/"/g, '\\"')}"`;
        }
        else
          prop.initialValue = String(pm.defaultValue);
      }
    }
  }

  return func;
}
