import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {historyUtils} from '../../../history-utils';
import {PipelineSerializedState} from '../config/PipelineInstance';
import {FuncCallAdapter, IFuncCallAdapter} from './FuncCallAdapters';

const RESTRICTIONS_PATH = 'INPUT_RESTRICTIONS';
const CONFIG_PATH = 'PIPELINE_CONFIG';

export async function makeFuncCall(nqName: string) {
  const func: DG.Func = await grok.functions.eval(nqName);
  const fc = func.prepare({});
  fc.newId();
  const bridge = new FuncCallAdapter(fc);
  return bridge;
}

export async function saveFuncCall(bridge: IFuncCallAdapter) {
  const fc = bridge.getFuncCall();
  // TODO: DF handling (?)
  const restrictions = bridge.inputRestrictions$.value;
  fc.options[RESTRICTIONS_PATH] = JSON.stringify(restrictions);
  return historyUtils.saveRun(fc);
}

export async function loadFuncCall(id: string): Promise<IFuncCallAdapter> {
  const fc = await historyUtils.loadRun(id, false);
  const restrictions = JSON.parse(fc.options[RESTRICTIONS_PATH] ?? {});
  const adapter = new FuncCallAdapter(fc);
  adapter.inputRestrictions$.next(restrictions);
  return adapter;
}

export async function saveInstanceState(nqName: string, stateJson: string) {
  const metaCall = await makeFuncCall(nqName);
  metaCall.instance.options[CONFIG_PATH] = stateJson;
  const fc = await historyUtils.saveRun(metaCall.instance);
  return fc.id;
}

export async function loadInstanceState(id: string) {
  const fc = await historyUtils.loadRun(id, false);
  const config = JSON.parse(fc.options[CONFIG_PATH] ?? {});
  return config as PipelineSerializedState;
}
