import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {historyUtils} from '../../../history-utils';
import {PipelineSerializedState} from '../config/PipelineInstance';
import {FuncCallAdapter} from './FuncCallAdapters';
import {serialize, deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {FuncCallInstancesBridge, RestrictionState} from './FuncCallInstancesBridge';
import {AdapterInitData} from './StateTreeNodes';
import {ItemMetadata} from '../view/ViewCommunication';
import {loadIsFavorite, saveIsFavorite} from '../../../shared-utils/history';
import dayjs from 'dayjs';

export const RESTRICTIONS_PATH = 'INPUT_RESTRICTIONS';
export const OUTPUT_OUTDATED_PATH = 'OUTPUT_OUTDATED';
export const RUN_ERROR_PATH = 'RUN_ERROR';

const CONFIG_PATH = 'PIPELINE_CONFIG';

// Marker key for a consistency DataFrame stored outside the restrictions blob as a table id.
const DF_REF_TOKEN = '_DG_CONSISTENCY_DF_REF_';

async function getAllUsersGroup() {
  return grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
}

// Stores restriction DataFrames as uploaded tables, replacing each with an id reference, so the
// serialized blob stays small. Scalar restrictions are serialized as before.
export async function serializeRestrictions(restrictions: Record<string, RestrictionState | undefined>) {
  const out: Record<string, RestrictionState | undefined> = {};
  let allGroup: DG.Group | undefined;
  for (const [name, restriction] of Object.entries(restrictions)) {
    if (restriction && restriction.assignedValue instanceof DG.DataFrame) {
      const df = restriction.assignedValue;
      const id = await grok.dapi.tables.uploadDataFrame(df);
      allGroup ??= await getAllUsersGroup();
      await grok.dapi.permissions.grant(df.getTableInfo(), allGroup, false);
      out[name] = {type: restriction.type, assignedValue: {[DF_REF_TOKEN]: id}};
    } else
      out[name] = restriction;
  }
  return serialize(out, {useJsonDF: false});
}

// Reads both the new id-reference format and the legacy embedded-DataFrame format. deserialize()
// reconstructs old embedded DataFrame tokens directly; id references are resolved via the tables API.
export async function deserializeRestrictions(str: string): Promise<Record<string, RestrictionState | undefined>> {
  const restrictions: Record<string, RestrictionState | undefined> = deserialize(str);
  for (const restriction of Object.values(restrictions)) {
    const ref = restriction?.assignedValue?.[DF_REF_TOKEN];
    if (ref)
      restriction!.assignedValue = await grok.dapi.tables.getTable(ref);
  }
  return restrictions;
}

// TODO: probably move to core
export function getFuncallDefaults(func: DG.Func) {
  const defaultValues: Record<string, any> = {};
  for (const prop of func.inputs) {
    const name = prop.name;
    if (prop.options.default)
    // sometimes it is parsable, sometimes not
    {
      try {
        defaultValues[name] = JSON.parse(prop.options.default);
      } catch {
        defaultValues[name] = prop.options.default;
      }
    }
  }
  return defaultValues;
}

export async function makeFuncCall(
  nqName: string, isReadonly: boolean,
): Promise<AdapterInitData> {
  const func = DG.Func.byName(nqName);
  const defaultValues = getFuncallDefaults(func);
  const fc = func.prepare(defaultValues);
  fc.newId();
  const adapter = new FuncCallAdapter(fc, isReadonly);
  return {adapter, restrictions: {}, runError: undefined, isOutputOutdated: true};
}

export async function saveFuncCall(bridge: FuncCallInstancesBridge) {
  const fc = bridge.getInstance()?.getFuncCall();
  if (!fc)
    throw new Error(`Attempting to save an empty FuncCall adapter`);
  if (bridge.isRunning$.value)
    throw new Error(`FuncCall node saving during being updated`);

  fc.options[RESTRICTIONS_PATH] = await serializeRestrictions(bridge.inputRestrictions$.value);
  fc.options[OUTPUT_OUTDATED_PATH] = serialize(bridge.isOutputOutdated$.value, {useJsonDF: false});
  fc.options[RUN_ERROR_PATH] = bridge.runError$.value;

  fc.newId();
  await historyUtils.saveRun(fc);
  return fc;
}

export async function loadFuncCall(
  id: string,
  isReadonly: boolean,
): Promise<AdapterInitData> {
  const fc = await historyUtils.loadRun(id);
  const restrictions = await deserializeRestrictions(fc.options[RESTRICTIONS_PATH] ?? '{}');
  const isOutputOutdated = deserialize(fc.options[OUTPUT_OUTDATED_PATH] ?? 'false');
  const runError = fc.options[RUN_ERROR_PATH];
  const adapter = new FuncCallAdapter(fc, isReadonly);
  return {adapter, restrictions, runError, isOutputOutdated};
}

export async function makeMetaCall(nqName: string) {
  const func = DG.Func.byName(nqName);
  const metaCall = func.prepare({});
  return metaCall;
}

export async function saveInstanceState(
  nqName: string,
  state: any,
  metaData?: ItemMetadata,
  version?: string,
) {
  const metaCall = await makeMetaCall(nqName);
  metaCall.options[CONFIG_PATH] = serialize(state, {useJsonDF: false});
  if (metaData?.title) metaCall.options['title'] = metaData.title;
  if (metaData?.description) metaCall.options['description'] = metaData.description;
  if (metaData?.tags) metaCall.options['tags'] = metaData.tags;
  if (version) metaCall.options['version'] = version;
  metaCall.newId();
  await metaCall.call(undefined, undefined, {processed: true, report: false});
  metaCall.started = dayjs();
  await historyUtils.saveRun(metaCall);
  await saveIsFavorite(metaCall, metaData?.isFavorite ?? false);
  return metaCall;
}

export async function loadInstanceState(id: string) {
  const metaCall = await historyUtils.loadRun(id);
  const isFavorite = await loadIsFavorite(metaCall);
  const config: PipelineSerializedState = deserialize(metaCall.options[CONFIG_PATH] ?? '{}');
  return [config, metaCall, isFavorite] as const;
}
