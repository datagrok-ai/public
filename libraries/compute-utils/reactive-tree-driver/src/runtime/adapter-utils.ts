import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {historyUtils} from '../../../history-utils';
import {PipelineSerializedState} from '../config/PipelineInstance';
import {FuncCallAdapter, IFuncCallAdapter} from './FuncCallAdapters';

const RESTRICTIONS_PATH = 'INPUT_RESTRICTIONS';
const OUTPUT_OUTDATED_PATH = 'OUTPUT_OUTDATED';
const CONFIG_PATH = 'PIPELINE_CONFIG';

export async function makeFuncCall(nqName: string) {
  const func: DG.Func = await grok.functions.eval(nqName);
  const fc = func.prepare({});
  fc.newId();
  const bridge = new FuncCallAdapter(fc);
  return bridge;
}

export async function saveFuncCall(adapter: IFuncCallAdapter) {
  const fc = adapter.getFuncCall();
  // TODO: DF restrictions handling (?)
  const restrictions = adapter.inputRestrictions$.value;
  const outputState = adapter.isOutputOutdated$.value;
  fc.options[RESTRICTIONS_PATH] = JSON.stringify(restrictions);
  fc.options[OUTPUT_OUTDATED_PATH] = JSON.stringify(outputState);
  fc.newId();
  return historyUtils.saveRun(fc);
}

export async function loadFuncCall(id: string): Promise<IFuncCallAdapter> {
  const fc = await historyUtils.loadRun(id, false);
  const restrictions = JSON.parse(fc.options[RESTRICTIONS_PATH] ?? {});
  const outputState = JSON.parse(fc.options[OUTPUT_OUTDATED_PATH] ?? false);
  const adapter = new FuncCallAdapter(fc);
  adapter.inputRestrictions$.next(restrictions);
  adapter.isOutputOutdated$.next(outputState);
  return adapter;
}

export async function saveInstanceState(nqName: string, stateJson: string) {
  const metaCall = await makeFuncCall(nqName);
  metaCall.instance.options[CONFIG_PATH] = stateJson;
  metaCall.instance.newId();
  const fc = await historyUtils.saveRun(metaCall.instance);
  return fc.id;
}

export async function loadInstanceState(id: string) {
  const fc = await historyUtils.loadRun(id, false);
  const config = JSON.parse(fc.options[CONFIG_PATH] ?? {});
  return config as PipelineSerializedState;
}
