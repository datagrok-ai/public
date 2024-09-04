import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {historyUtils} from '../../../history-utils';
import {PipelineSerializedState} from '../config/PipelineInstance';
import {FuncCallAdapter, IFuncCallAdapter} from './FuncCallAdapters';
import {serialize, deserialize} from '@datagrok-libraries/utils/src/json-serialization';

const RESTRICTIONS_PATH = 'INPUT_RESTRICTIONS';
const OUTPUT_OUTDATED_PATH = 'OUTPUT_OUTDATED';
const CONFIG_PATH = 'PIPELINE_CONFIG';

export async function makeFuncCall(nqName: string, isReadonly: boolean) {
  const func: DG.Func = await grok.functions.eval(nqName);
  const fc = func.prepare({});
  fc.newId();
  const bridge = new FuncCallAdapter(fc, isReadonly);
  return bridge;
}

export async function saveFuncCall(adapter: IFuncCallAdapter) {
  const fc = adapter.getFuncCall();
  // TODO: DF restrictions better handling
  const restrictions = adapter.inputRestrictions$.value;
  const outputState = adapter.isOutputOutdated$.value;
  fc.options[RESTRICTIONS_PATH] = serialize(restrictions, {useJsonDF: true});
  fc.options[OUTPUT_OUTDATED_PATH] = serialize(outputState, {useJsonDF: true});
  fc.newId();
  return historyUtils.saveRun(fc);
}

export async function loadFuncCall(id: string, isReadonly: boolean): Promise<IFuncCallAdapter> {
  const fc = await historyUtils.loadRun(id, false);
  const restrictions = deserialize(fc.options[RESTRICTIONS_PATH] ?? '{}');
  const outputState = deserialize(fc.options[OUTPUT_OUTDATED_PATH] ?? 'false');
  const adapter = new FuncCallAdapter(fc, isReadonly);
  adapter.inputRestrictions$.next(restrictions);
  adapter.isOutputOutdated$.next(outputState);
  return adapter;
}

export async function makeMetaCall(nqName: string) {
  const func: DG.Func = await grok.functions.eval(nqName);
  const metaCall = func.prepare({});
  return metaCall;
}

export async function saveInstanceState(nqName: string, state: any, currentMetaCall?: DG.FuncCall) {
  const metaCall = currentMetaCall ?? await makeMetaCall(nqName);
  metaCall.options[CONFIG_PATH] = serialize(state, {useJsonDF: true});
  metaCall.newId();
  const fc = await historyUtils.saveRun(metaCall);
  return fc;
}

export async function loadInstanceState(id: string) {
  const metaCall = await historyUtils.loadRun(id, false);
  const config: PipelineSerializedState = deserialize(metaCall.options[CONFIG_PATH] ?? '{}');
  return [metaCall, config] as const;
}
