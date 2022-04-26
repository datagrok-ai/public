/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HistoricalRunService} from '../common/service-interfaces';

export class DefaultHistoricalRunService implements HistoricalRunService {
  saveRun(call: DG.FuncCall): Promise<string> {
    throw new Error('Method is not implemented');
    /* await grok.dapi.functions.calls.save(call);*/
  }
  loadRun(runId: string): Promise<DG.FuncCall> {
    throw new Error('Method is not implemented');
    /* await grok.dapi.functions.calls.find(call.id);*/
  }
  pullRuns(): Promise<DG.FuncCall[]> {
    throw new Error('Method is not implemented');
    /* await grok.dapi.functions.calls.list();*/
  }
}
