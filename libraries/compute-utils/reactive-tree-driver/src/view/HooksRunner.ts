import * as DG from 'datagrok-api/dg';
import {PipelineState} from '../config/PipelineInstance';
import {HooksSpec, PipelineHooks} from '../config/PipelineConfiguration';
import {callHandler} from '../utils';
import {Observable, concat} from 'rxjs';
import {ItemPathArray} from '../data/common-types';

export class HooksRunner {
  constructor(private rt: any, private hooks: HooksSpec[]) {}

  onInit() {
    return this.execHooks('onInit');
  }

  beforeLoadFuncCall(id: string, editorFunc: DG.Func) {
    return this.execHooks('beforeLoadFuncCall', {id, editorFunc});
  }

  afterLoadFuncCall(id: string, funcCall: DG.FuncCall) {
    return this.execHooks('afterLoadFuncCall', {id, funcCall});
  }

  beforeInputFormRender(id: string, funcCall: DG.FuncCall) {
    return this.execHooks('beforeInputFormRender', {id, funcCall});
  }

  afterInputFormRender(id: string, inputForm: DG.InputForm, funcCall: DG.FuncCall) {
    return this.execHooks('afterInputFormRender', {id, inputForm, funcCall});
  }

  beforeViewerRender(id: string, name: string, funcCall: DG.FuncCall) {
    return this.execHooks('beforeViewerRender', {id, name, funcCall});
  }

  afterViewerRender(id: string, name: string, viewer: DG.Viewer, funcCall: DG.FuncCall) {
    return this.execHooks('afterViewerRender', {id, name, viewer, funcCall});
  }

  beforeLoadRun(id: string) {
    return this.execHooks('beforeLoadRun', {id});
  }

  afterLoadRun(id: string, funcCall: DG.FuncCall, config: PipelineState) {
    return this.execHooks('afterLoadRun', {id, funcCall, config});
  }

  beforeSaveRun(funcCall: DG.FuncCall, config: PipelineState) {
    return this.execHooks('beforeSaveRun', {funcCall, config});
  }

  afterSaveRun(funcCall: DG.FuncCall, config: PipelineState) {
    return this.execHooks('afterSaveRun', {funcCall, config});
  }

  onClose() {
    this.rt!.pipelineState.closed.next(true);
    return this.execHooks('onClose');
  }

  private execHooks(category: keyof PipelineHooks<ItemPathArray>, additionalParams: Record<string, any> = {}) {
    const observables: Observable<any>[] = [];
    for (const {hooks, path} of this.hooks!) {
      const items = hooks[category];
      for (const item of items ?? []) {
        const handler = item.handler!;
        // const ctrlConf = new ControllerConfig(path, item.from, item.to);
        // const controller = new RuntimeControllerImpl(item.id, ctrlConf, this.rt!);
        // const params = {...additionalParams, controller};
        // const obs$ = callHandler(handler, params);
        // observables.push(obs$);
      }
    }
    return concat(observables);
  }
}
