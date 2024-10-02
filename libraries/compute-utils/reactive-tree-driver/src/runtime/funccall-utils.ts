import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {historyUtils} from '../../../history-utils';
import {PipelineSerializedState} from '../config/PipelineInstance';
import {FuncCallAdapter, IFuncCallAdapter} from './FuncCallAdapters';
import {serialize, deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {FuncCallInstancesBridge, RestrictionState} from './FuncCallInstancesBridge';

const RESTRICTIONS_PATH = 'INPUT_RESTRICTIONS';
const OUTPUT_OUTDATED_PATH = 'OUTPUT_OUTDATED';
const CONFIG_PATH = 'PIPELINE_CONFIG';

export async function makeFuncCall(
  nqName: string, isReadonly: boolean,
): Promise<[IFuncCallAdapter, Record<string, RestrictionState>, boolean]> {
  const func: DG.Func = await grok.functions.eval(nqName);
  const fc = func.prepare({});
  fc.newId();
  const adapter = new FuncCallAdapter(fc, isReadonly);
  return [adapter, {}, true];
}

export async function saveFuncCall(bridge: FuncCallInstancesBridge) {
  const fc = bridge.getInstance()?.getFuncCall();
  if (!fc)
    throw new Error(`Attempting to save an empty FuncCall adapter`);
  if (bridge.isRunning$.value)
    throw new Error(`FuncCall node saving during being updated`);

  // TODO: DF restrictions better handling
  fc.options[RESTRICTIONS_PATH] = serialize(bridge.inputRestrictions$.value, {useJsonDF: true});
  fc.options[OUTPUT_OUTDATED_PATH] = serialize(bridge.isOutputOutdated$.value, {useJsonDF: true});
  fc.newId();
  await historyUtils.saveRun(fc);
  return fc;
}

export async function loadFuncCall(
  id: string,
  isReadonly: boolean,
): Promise<[IFuncCallAdapter, Record<string, RestrictionState | undefined>, boolean]> {
  const fc = await historyUtils.loadRun(id, false);
  const restrictions = deserialize(fc.options[RESTRICTIONS_PATH] ?? '{}');
  const outputState = deserialize(fc.options[OUTPUT_OUTDATED_PATH] ?? 'false');
  const adapter = new FuncCallAdapter(fc, isReadonly);
  return [adapter, restrictions, outputState];
}

export async function makeMetaCall(nqName: string) {
  const func: DG.Func = await grok.functions.eval(nqName);
  const metaCall = func.prepare({});
  return metaCall;
}

export async function saveInstanceState(nqName: string, state: any) {
  const metaCall = await makeMetaCall(nqName);
  metaCall.options[CONFIG_PATH] = serialize(state, {useJsonDF: true});
  metaCall.newId();
  await historyUtils.saveRun(metaCall);
  return metaCall;
}

export async function loadInstanceState(id: string) {
  const metaCall = await historyUtils.loadRun(id, false);
  const config: PipelineSerializedState = deserialize(metaCall.options[CONFIG_PATH] ?? '{}');
  return config;
}
