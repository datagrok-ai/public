import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {VIEWER_PATH, viewerTypesMapping} from './consts';
import wu from 'wu';

export const deepCopy = (call: DG.FuncCall) => {
  const deepClone = call.clone();

  const dfOutputs = wu(call.outputParams.values() as DG.FuncCallParam[])
    .filter((output) => output.property.propertyType === DG.TYPE.DATA_FRAME);
  for (const output of dfOutputs)
    deepClone.outputs[output.name] = call.outputs[output.name].clone();

  const dfInputs = wu(call.inputParams.values() as DG.FuncCallParam[])
    .filter((input) => input.property.propertyType === DG.TYPE.DATA_FRAME);
  for (const input of dfInputs)
    deepClone.inputs[input.name] = call.inputs[input.name].clone();

  return deepClone;
};

export const boundImportFunction = (func: DG.Func): string | undefined => {
  return func.options['getRealData'];
};

export const getPropViewers = (prop: DG.Property): {name: string, config: Record<string, string | boolean>[]} => {
  const viewersRawConfig = prop.options[VIEWER_PATH];
  return (viewersRawConfig !== undefined) ?
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

export const getFuncRunLabel = (func: DG.Func) => {
  return func.options['runLabel'];
};
