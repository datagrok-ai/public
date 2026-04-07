import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {IKnimeClient} from './knime-client';
import {KnimeInputParam, KnimeOutputParam, KnimeParamType, KnimeDeployment, KnimeExecutionResult, KnimeWorkflowSpec} from './types';
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
  tableSpec?: {name: string; type: string}[];
}

const knimeInputTypeToDgType: Record<KnimeParamType, string> = {
  'table': 'dataframe',
  'file': 'file',
  'string': 'string',
  'int': 'int',
  'double': 'double',
  'boolean': 'bool',
  'json': 'string',
};

const knimeOutputTypeToDgType: Record<KnimeParamType, string> = {
  'table': 'dataframe',
  'file': 'string',
  'string': 'string',
  'int': 'int',
  'double': 'double',
  'boolean': 'bool',
  'json': 'string',
};

const MOLECULE_KEYWORDS = ['smiles', 'molecule', 'mol', 'compound', 'structure', 'chemical'];

function containsMoleculeKeyword(name: string): boolean {
  const lower = name.toLowerCase();
  return MOLECULE_KEYWORDS.some((kw) => lower.includes(kw));
}

export function sanitizeFuncName(name: string, id: string): string {
  let result = name.replace(/[^a-zA-Z0-9_]/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '');
  if (!result)
    result = id.replace(/[^a-zA-Z0-9_]/g, '_').replace(/_+/g, '_').substring(0, 20);
  return 'Knime_' + result;
}

export function sanitizeParamName(name: string): string {
  let result = name.replace(/[^a-zA-Z0-9_]/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '');
  if (!result || /^\d/.test(result))
    result = 'p_' + result;
  return result;
}

export function buildParamMeta(params: KnimeInputParam[]): ParamMeta[] {
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
      tableSpec: param.tableSpec,
    });
  }

  return meta;
}

export interface OutputMeta {
  sanitizedName: string;
  originalName: string;
  dgType: string;
}

export function buildOutputMeta(outputs: KnimeOutputParam[]): OutputMeta[] {
  const meta: OutputMeta[] = [];
  const usedNames = new Set<string>();
  for (const out of outputs) {
    let name = sanitizeParamName(out.name);
    let n = 2;
    const baseName = name;
    while (usedNames.has(name)) {
      name = `${baseName}_${n}`;
      n++;
    }
    usedNames.add(name);
    meta.push({sanitizedName: name, originalName: out.name, dgType: knimeOutputTypeToDgType[out.type] ?? 'object'});
  }
  return meta;
}

export function buildSignature(funcName: string, paramMeta: ParamMeta[], outputs: KnimeOutputParam[]): string {
  const parts: string[] = [];
  for (const pm of paramMeta) {
    const dgType = knimeInputTypeToDgType[pm.type] ?? 'string';
    parts.push(`${dgType} ${pm.sanitizedName}`);
  }
  const returnType = outputs.length === 1
    ? (knimeOutputTypeToDgType[outputs[0].type] ?? 'object')
    : 'object';
  return `${returnType} ${funcName}(${parts.join(', ')})`;
}

function createRunCallback(
  deployment: KnimeDeployment,
  paramMeta: ParamMeta[],
  outputMeta: OutputMeta[],
  client: IKnimeClient,
): Function {
  return async (...args: any[]) => {
    // TODO: table schema validation disabled — the tableSpec comes from OpenAPI examples
    // and doesn't reliably represent the actual required schema.
    // for (let i = 0; i < paramMeta.length; i++) {
    //   const pm = paramMeta[i];
    //   const val = args[i];
    //   if (pm.type === 'table' && pm.tableSpec && pm.tableSpec.length > 0 && val instanceof DG.DataFrame) {
    //     const errors: string[] = [];
    //     for (const expected of pm.tableSpec) {
    //       const col = val.columns.byName(expected.name);
    //       if (!col) {
    //         errors.push(`Missing column "${expected.name}"`);
    //         continue;
    //       }
    //       const compatibleTypes = knimeTypeToDgColumnTypes(expected.type);
    //       if (!compatibleTypes.includes(col.type))
    //         errors.push(`Column "${expected.name}": expected ${expected.type}, got ${col.type}`);
    //     }
    //     if (errors.length > 0)
    //       throw new Error(`Input "${pm.originalName}" schema mismatch: ${errors.join('; ')}`);
    //   }
    // }

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

    return extractReturnValue(result, outputMeta);
  };
}

function parseTable(t: any): DG.DataFrame {
  if (t?.['table-spec'] && t?.['table-data'])
    return knimeSpecDataToDataFrame(t['table-spec'], t['table-data']);
  return Array.isArray(t) ? knimeTableToDataFrame(t) : knimeTableToDataFrame([t]);
}

function parseResource(fileName: string, resource: any): DG.DataFrame | null {
  if (resource.json !== undefined && Array.isArray(resource.json)) {
    const df = knimeTableToDataFrame(resource.json, fileName);
    df.name = fileName;
    return df;
  }
  if (resource.text !== undefined) {
    const isCsv = resource.contentType.includes('csv') || fileName.endsWith('.csv');
    const isTsv = resource.contentType.includes('tab-separated') || fileName.endsWith('.tsv');
    if (isCsv || isTsv) {
      try {
        const df = DG.DataFrame.fromCsv(resource.text);
        df.name = fileName;
        return df;
      }
      catch { /* not parseable */ }
    }
  }
  return null;
}

function extractReturnValue(result: KnimeExecutionResult, outputMeta: OutputMeta[]): any {
  // Build named maps of parsed outputs
  const parsedTables = new Map<string, DG.DataFrame>();
  const parsedScalars = new Map<string, any>();

  // Parse outputValues — inline key-value results
  if (result.outputs)
    for (const [key, val] of Object.entries(result.outputs)) {
      if (val && typeof val === 'object' && val['table-spec'] && val['table-data'])
        parsedTables.set(key, parseTable(val));
      else
        parsedScalars.set(key, val);
    }

  // Parse outputTables (unnamed, positional)
  const unnamedTables: DG.DataFrame[] = [];
  if (result.outputTables)
    for (const t of result.outputTables)
      unnamedTables.push(parseTable(t));

  // Parse fetchedResources (named)
  if (result.outputResources)
    for (const fileName of Object.keys(result.outputResources)) {
      const resource = result.outputResources[fileName];
      if (typeof resource === 'string') {
        parsedScalars.set(fileName, resource);
        continue;
      }
      const df = parseResource(fileName, resource);
      if (df)
        parsedTables.set(fileName, df);
    }

  // Single output — return the value directly for backward compatibility
  if (outputMeta.length === 1) {
    const om = outputMeta[0];
    if (om.dgType === 'dataframe')
      return parsedTables.get(om.originalName) ?? unnamedTables[0] ?? DG.DataFrame.create();
    return parsedScalars.get(om.originalName) ?? null;
  }

  // No known outputs — legacy fallback
  if (outputMeta.length === 0) {
    if (unnamedTables.length > 0)
      return unnamedTables[0];
    if (parsedTables.size > 0)
      return parsedTables.values().next().value;
    return result.outputs;
  }

  // Multiple outputs — return named object keyed by sanitized output names
  const ret: {[key: string]: any} = {};
  let tableIdx = 0;
  for (const om of outputMeta) {
    if (om.dgType === 'dataframe') {
      ret[om.sanitizedName] = parsedTables.get(om.originalName)
        ?? unnamedTables[tableIdx++]
        ?? DG.DataFrame.create();
    }
    else
      ret[om.sanitizedName] = parsedScalars.get(om.originalName) ?? null;
  }
  return ret;
}

export async function registerFuncFromSpec(
  deployment: KnimeDeployment,
  spec: KnimeWorkflowSpec,
  client: IKnimeClient,
): Promise<DG.Func> {
  const funcName = sanitizeFuncName(deployment.name, deployment.id);
  const paramMeta = buildParamMeta(spec.inputs);
  const outputMeta = buildOutputMeta(spec.outputs);
  const signature = buildSignature(funcName, paramMeta, spec.outputs);
  const runCallback = createRunCallback(deployment, paramMeta, outputMeta, client);

  const registrationOutputs = outputMeta.length > 1
    ? outputMeta.map((om) => ({name: om.sanitizedName, type: om.dgType}))
    : undefined;

  const func = grok.functions.register({
    signature,
    run: runCallback,
    isAsync: true,
    namespace: 'KnimeLink',
    options: {'engine': 'KNIME'},
    outputs: registrationOutputs,
  });

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
      if (pm.type === 'string' && containsMoleculeKeyword(pm.originalName)) {
        if (pm.defaultValue != null) {
          try {
            if (await grok.functions.call('Chem:isSmiles', {s: String(pm.defaultValue)}))
              prop.semType = 'Molecule';
          }
          catch { /* Chem package not available */ }
        }
        else
          prop.semType = 'Molecule';
      }
    }
  }

  return func;
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
  return registerFuncFromSpec(deployment, spec, client);
}
