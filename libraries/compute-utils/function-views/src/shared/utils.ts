import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DIRECTION, VIEWER_PATH, viewerTypesMapping} from './consts';

// DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13015
export const getDataFrame = (call: DG.FuncCall, name: string, direction: DIRECTION): DG.DataFrame => {
  if (direction === DIRECTION.INPUT)
    return call.inputs[name] instanceof DG.DataFrame ? call.inputs[name] : call.inputs[name].dataFrame;
  else
    return call.outputs[name] instanceof DG.DataFrame ? call.outputs[name] : call.outputs[name].dataFrame;
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

      return v;
    })}:
    {name: prop.name, config: []};
};

export const getFuncFriendlyName = (func: DG.Func) => {
  return func.options['friendlyName'];
};

export const getFuncRunLabel = (func: DG.Func) => {
  return func.options['runLabel'];
};
