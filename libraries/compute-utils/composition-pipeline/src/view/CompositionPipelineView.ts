import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Observable} from 'rxjs';
import {PipelineView, RichFunctionView} from '../../../function-views';
import {takeUntil} from 'rxjs/operators';
import {ValidationResult} from '../../../shared-utils/validation';
import {ABILITY_STATE, VISIBILITY_STATE} from '../../../shared-utils/consts';
import {StepState} from '../../../function-views/src/pipeline-view';
import {ItemPath, PipelineHooks, ExportConfig, NqName} from '../PipelineConfiguration';
import {PathKey} from '../config-processing-utils';
import {ControllerConfig} from '../runtime/ControllerConfig';
import {PipelineRuntime} from '../runtime/PipelineRuntime';
import {RuntimeControllerImpl} from '../runtime/RuntimeControllerImpl';
import {callHandler} from '../utils';
import {RFVPopup} from './RFVPopup';

export interface StepSpec {
  id: string;
  funcName: string;
  friendlyName?: string;
  helpUrl?: string | HTMLElement;
}
export interface HookSpec {
  hooks: PipelineHooks;
  pipelinePath: ItemPath;
}

export interface ICompositionView {
  injectConfiguration(steps: StepSpec[], hooks: HookSpec[], rt: PipelineRuntime, exportConfig?: ExportConfig): void;
  getStateBindings<T = any>(viewId: PathKey, stateId: string): { changes: Observable<T>; setter: (x: any) => void; };
  enableStep(stepId: PathKey): void;
  getStepState(stepId: PathKey): VISIBILITY_STATE;
  getFollowingStep(): StepState | undefined;
  showSteps(...id: NqName[]): void;
  hideSteps(...id: NqName[]): void;
  getStepView<T = RichFunctionView>(viewId?: PathKey): T;
  setExternalValidationResults(viewId: PathKey, stateId: string, results: ValidationResult): void;
  isUpdating: BehaviorSubject<boolean>;
  getRunningUpdates(): string[];
  addMenuItem(name: string, action: () => void): void;
}

export class CompositionPipelineView extends PipelineView implements ICompositionView {
  private hooks: HookSpec[] = [];
  private rt?: PipelineRuntime;
  private actionsMenu?: DG.Menu;
  public customViews = new Map<string, RFVPopup>();

  constructor(funcName: string, options: {
    historyEnabled: boolean,
    isTabbed: boolean,
    skipInit?: boolean,
  } = {historyEnabled: true, isTabbed: false, skipInit: true}) {
    super(funcName, [], options);
  }

  public injectConfiguration(steps: StepSpec[], hooks: HookSpec[], rt: PipelineRuntime, exportConfig?: ExportConfig) {
    if (exportConfig)
      this.exportConfig = {...this.exportConfig, ...exportConfig};

    this.initialConfig = steps.map((step) => ({...step, customId: step.id}));
    this.hooks = hooks;
    this.rt = rt;
    this.rt.isUpdating.pipe(takeUntil(this.rt.pipelineState.closed)).subscribe((val) => this.isUpdating.next(val));
  }

  public getStateBindings(viewId: PathKey, stateId: string) {
    const stepView = this.getStepView<RichFunctionView>(viewId);
    const changes = stepView.getParamChanges(stateId);
    const setter = (x: any, inputState?: 'disabled' | 'restricted' | 'user input') => {
      stepView.setInput(stateId, x, inputState);
    };
    return {changes, setter};
  }

  public setExternalValidationResults(viewId: PathKey, inputName: string, results: ValidationResult) {
    const view = this.getStepView(viewId);
    if (view)
      view.setExternalValidationResults(inputName, results);
  }

  public addMenuItem(name: string, action: () => void) {
    this.actionsMenu?.item(name, action);
  }

  public enableStep(stepId: PathKey) {
    this.steps[stepId].ability.next(ABILITY_STATE.ENABLED);
  }

  public getStepState(stepId: PathKey): VISIBILITY_STATE {
    return this.steps[stepId].visibility.value;
  }

  public getFollowingStep() {
    const disabledSteps = Object.values(this.steps)
      .filter((stepData) => stepData.visibility.value === VISIBILITY_STATE.VISIBLE && stepData.ability.value === ABILITY_STATE.DISABLED);
    return disabledSteps[0];
  }

  override getRunningUpdates() {
    return this.rt!.getRunningUpdates();
  }

  override buildRibbonMenu() {
    super.buildRibbonMenu();
    this.actionsMenu = this.ribbonMenu.group('Actions');
  }

  override getStepView<T = RichFunctionView>(viewId: PathKey) {
    const view = super.getStepView<RichFunctionView>(viewId);
    if (!view)
      return this.customViews.get(viewId)! as T;

    return view as T;
  }

  override async init() {
    this.rt!.disableSubStepsIOSetters([]);
    await this.execHooks('beforeInit');
    await super.init();
    await this.execHooks('afterInit');
  }

  override async onBeforeStepFuncCallApply(nqName: string, scriptCall: DG.FuncCall, editorFunc: DG.Func) {
    await this.execHooks('beforeFuncCallReady', {nqName, scriptCall, editorFunc});
  }

  override async onAfterStepFuncCallApply(nqName: string, scriptCall: DG.FuncCall, view: RichFunctionView) {
    await this.execHooks('afterFuncCallReady', {nqName, scriptCall, view});
  }

  override async onBeforeLoadRun() {
    this.rt!.disableSubStepsIOSetters([]);
    await super.onBeforeLoadRun();
    await this.execHooks('beforeLoadRun');
  }

  override async onAfterLoadRun(run: DG.FuncCall) {
    await super.onAfterLoadRun(run);
    await this.rt!.awaitForAllUpdatesDone();
    this.rt!.enableSubStepsIOSetters([]);
    await this.execHooks('afterLoadRun', {run});
  }

  override async onBeforeSaveRun() {
    this.rt!.disableSubStepsIOSetters([]);
    await super.onBeforeLoadRun();
    await this.execHooks('beforeSaveRun');
  }

  override async onAfterSaveRun(run: DG.FuncCall) {
    await super.onAfterLoadRun(run);
    await this.rt!.awaitForAllUpdatesDone();
    this.rt!.enableSubStepsIOSetters([]);
    await this.execHooks('afterSaveRun', {run});
  }

  override build() {
    super.build();
    this.rt!.wireViews();
    this.rt!.wireLinks();
    this.rt!.enableSubStepsIOSetters([]);
    this.execHooks('onViewReady', {view: this});
  }

  override close() {
    this.rt!.pipelineState.closed.next(true);
    for (const v of this.customViews.values())
      v.customClose();

    super.close();
  }

  private async execHooks(category: keyof PipelineHooks, additionalParams: Record<string, any> = {}) {
    for (const {hooks, pipelinePath} of this.hooks!) {
      const items = hooks[category];
      for (const item of items ?? []) {
        const handler = item.handler!;
        const ctrlConf = new ControllerConfig(pipelinePath, item.from, item.to);
        const controller = new RuntimeControllerImpl(item.id, ctrlConf, this.rt!);
        const params = {...additionalParams, controller};
        await callHandler(handler, params);
      }
    }
  }
}
'';
