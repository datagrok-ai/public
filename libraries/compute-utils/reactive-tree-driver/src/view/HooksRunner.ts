import * as DG from 'datagrok-api/dg';
import {PipelineRuntime} from '../runtime/PipelineRuntime';
import {HookSpec} from './CompositionPipelineView';
import {PipelineHooks} from '../config/PipelineConfiguration';
import {callHandler} from '../utils';
import {ControllerConfig} from '../runtime/ControllerConfig';
import {RuntimeControllerImpl} from '../runtime/RuntimeControllerImpl';
import {ViewConfig} from './ViewCommunication';
import {Observable, concat} from 'rxjs';

export class HooksRunner {
  constructor(private rt: PipelineRuntime, private hooks: HookSpec[]) {}

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

  afterLoadRun(id: string, funcCall: DG.FuncCall, config: ViewConfig) {
    return this.execHooks('afterLoadRun', {id, funcCall, config});
  }

  beforeSaveRun(funcCall: DG.FuncCall, config: ViewConfig) {
    return this.execHooks('beforeSaveRun', {funcCall, config});
  }

  afterSaveRun(funcCall: DG.FuncCall, config: ViewConfig) {
    return this.execHooks('afterSaveRun', {funcCall, config});
  }

  onClose() {
    this.rt!.pipelineState.closed.next(true);
    return this.execHooks('onClose');
  }

  private execHooks(category: keyof PipelineHooks, additionalParams: Record<string, any> = {}) {
    const observables: Observable<any>[] = [];
    for (const {hooks, pipelinePath} of this.hooks!) {
      const items = hooks[category];
      for (const item of items ?? []) {
        const handler = item.handler!;
        const ctrlConf = new ControllerConfig(pipelinePath, item.from, item.to);
        const controller = new RuntimeControllerImpl(item.id, ctrlConf, this.rt!);
        const params = {...additionalParams, controller};
        const obs$ = callHandler(handler, params);
        observables.push(obs$);
      }
    }
    return concat(observables);
  }
}
