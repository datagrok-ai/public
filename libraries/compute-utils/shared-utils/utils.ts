import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import wu from 'wu';
import {VIEWER_PATH, viewerTypesMapping} from './consts';
import {deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {OUTPUT_OUTDATED_PATH} from '../reactive-tree-driver/src/runtime/funccall-utils';

export const replaceForWindowsPath = (rawName: string, stringToInsert?: string) => {
  const regExpForWindowsPath = /(\/|\\|\:|\*|\?|\"|\<|\>|\|)/g;

  return rawName.replaceAll(regExpForWindowsPath, stringToInsert ?? '_');
};

const helpCache = new DG.LruCache<string, string>();

export const getContextHelp = async (func: DG.Func) => {
  const helpPath = func.options['help'] ?? func.options['readme'];

  if (!helpPath) return undefined;

  if (helpCache.get(func.id)) return helpCache.get(func.id);

  // some bug when saving package script from ui
  let packageName: string = '';
  try {
    packageName = func.package.name;
  } catch {}

  const packagePath = `System:AppData/${packageName}/${helpPath}`;
  if (packageName && await grok.dapi.files.exists(packagePath)) {
    const readme = await grok.dapi.files.readAsText(packagePath);
    helpCache.set(func.id, readme);
    return readme;
  }

  const homePath = `${grok.shell.user.name}.home/${helpPath}`;
  if (await grok.dapi.files.exists(homePath)) {
    const readme = await grok.dapi.files.readAsText(homePath);
    helpCache.set(func.id, readme);
    return readme;
  }

  return undefined;
};

export const hasContextHelp = (func?: DG.Func) => {
  return !!(func?.options?.['help']);
};

export const isRunningOnInput = (func: DG.Func) => {
  return func.options['runOnInput'] === 'true';
};

export const getFeatures = (func: DG.Func) => {
  return JSON.parse(func.options['features'] ?? '{}');
};

export const getDockSpawnConfig = (func: DG.Func) => {
  return JSON.parse(func.options['dockSpawnConfig'] ?? '{}');
};

export const getRunLabel = (func: DG.Func) => {
  return func.options['runLabel'];
};

export const getCustomExports = (func: DG.Func) => {
  return JSON.parse(func.options['customExports'] ?? '[]');
};

export const getFeature = (features: Record<string, boolean> | string[], featureName: string, defaultValue: boolean) => {
  if (features instanceof Array)
    return features.includes(featureName);

  if (features instanceof Object)
    return features[featureName] ?? defaultValue;

  return defaultValue;
};

export const getStartedOrNull = (run?: DG.FuncCall) => {
  try {
    return run?.started;
  } catch {
    return;
  }
};

export const isIncomplete = (run: DG.FuncCall) => {
  const isOutputOutdated = deserialize(run.options[OUTPUT_OUTDATED_PATH] ?? 'false');
  return isOutputOutdated || !getStartedOrNull(run) || !run.id;
};

export const deepCopy = (call: DG.FuncCall) => {
  // if there is no ':' in nqName, it might be ambiguous, so try to
  // reuse call.func, otherwise use nqName search
  const deepClone = call.func?.nqName?.includes(':') ? DG.Func.byName(call.func.nqName).prepare() : call.func.prepare();

  deepClone.id = call.id;

  call.options.forEach((key: string) => deepClone.options[key] = call.options[key]);

  const definedOutputs = wu(deepClone.outputParams.values())
    .filter((output) => call.outputs[output.name] != null);
  for (const output of definedOutputs) {
    if (output.property.propertyType === DG.TYPE.DATA_FRAME)
      deepClone.outputs[output.name] = call.outputs[output.name].clone();
    else
      deepClone.outputs[output.name] = call.outputs[output.name];
  }

  const definedInputs = wu(deepClone.inputParams.values())
    .filter((input) => call.inputs[input.name] != null);
  for (const input of definedInputs) {
    if (input.property.propertyType === DG.TYPE.DATA_FRAME)
      deepClone.inputs[input.name] = call.inputs[input.name].clone();
    else
      deepClone.inputs[input.name] = call.inputs[input.name];
  }

  if (getStartedOrNull(call)) deepClone.started = call.started;

  return deepClone;
};

export const getPropViewers = (prop: DG.Property): {name: string, config: Record<string, string | boolean>[]} => {
  const viewersRawConfig = prop.options[VIEWER_PATH];
  return viewersRawConfig ?
  // true and false values are retrieved as string, so we parse them separately
    {name: prop.name, config: JSON.parse(viewersRawConfig, (k, v) => {
      if (v === 'true') return true;
      if (v === 'false') return false;
      // Converting internal Dart labels to JS DG.VIEWER labels
      if (k === 'type') return viewerTypesMapping[v] || v;

      if (!k.toLowerCase().includes('color')) {
        const parsed = Number.parseFloat(v);

        if (!Number.isNaN(parsed))
          return parsed;
      }

      return v;
    })}:
    {name: prop.name, config: []};
};

export const getFuncCallDefaultFilename = (funcCall: DG.FuncCall) => {
  return `${funcCall.func.nqName} - ${getStartedOrNull(funcCall) ?? 'Not completed'}.xlsx`;
};

export const getStarted = (call: DG.FuncCall) => {
  try {
    return call.started.toDate()
      .toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'});
  } catch {
    return 'Not completed';
  }
};

export const delay = (delayInms: number) => {
  return new Promise((resolve) => setTimeout(resolve, delayInms));
};
