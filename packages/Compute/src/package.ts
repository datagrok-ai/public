/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';
import {
  ComputationView as ComputationViewInst,
  RichFunctionView as RichFunctionViewInst,
  PipelineView as PipelineViewInst,
  CustomFunctionView as CustomFunctionViewInst,
  UiUtils,
} from '@datagrok-libraries/compute-utils';
import {
  ValidationInfo,
  makeAdvice as makeAdviceInst,
  makeValidationResult as makeValidationResultInst,
  makeRevalidation as makeRevalidationInst,
  mergeValidationResults as mergeValidationResultsInst,
} from '@datagrok-libraries/compute-utils';
import {ModelHandler} from '@datagrok-libraries/compute-utils/model-catalog';
import {
  testPipeline as testPipelineInst,
} from '@datagrok-libraries/compute-utils';
import {
  deepCopy as  deepCopyInst,
} from '@datagrok-libraries/compute-utils';


export * from './package.g';
export const _package = new DG.Package();

// for compute-api pakage
export const testPipeline = testPipelineInst;
export const deepCopy = deepCopyInst;
export const CompView = ComputationViewInst;
export const RFV = RichFunctionViewInst;
export const CFV = CustomFunctionViewInst;
export const Pipeline = PipelineViewInst;
export const makeValidationResult = makeValidationResultInst;
export const makeAdvice = makeAdviceInst;
export const makeRevalidation = makeRevalidationInst;
export const mergeValidationResults = mergeValidationResultsInst;
export const fileInput = UiUtils.fileInput;
export const historyInput = UiUtils.historyInput;
export const historyInputJSON = UiUtils.historyInputJSON;
export const historyPanel = UiUtils.historyPanel;


export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
  }


  @grok.decorators.func()
  static openModelFromFuncall(funccall: DG.FuncCall) {
    ModelHandler.openModelFromFunccall(funccall as any);
  }


  @grok.decorators.func({
    name: 'OutliersSelectionViewer',
    description: 'Creates an outliers selection viewer',
    outputs: [{type: 'viewer', name: 'result'}],
    meta: {role: 'viewer'},
  })
  static OutliersSelection() {
    return new OutliersSelectionViewer();
  }


  @grok.decorators.editor({outputs: [{type: 'view', name: 'result'}]})
  static RichFunctionViewEditor(call: DG.FuncCall) {
    return RichFunctionViewInst.fromFuncCall(call as any, {historyEnabled: true, isTabbed: false});
  }


  @grok.decorators.editor({outputs: [{type: 'object', name: 'result'}]})
  static PipelineStepEditor(call: DG.FuncCall) {
    return RichFunctionViewInst.fromFuncCall(call as any, {historyEnabled: false, isTabbed: true});
  }

  //
  // Testing code
  //

  @grok.decorators.func({outputs: [{type: 'object', name: 'uploadedCalls'}]})
  static async CustomDataUploader(func: DG.Func) {
    await new Promise((r) => setTimeout(r, 1000));

    const dummyFunccall = await func.prepare({
      'ambTemp': 22,
      'initTemp': 100,
      'desiredTemp': 30,
      'area': 0.06,
      'heatCap': 4200,
      'heatTransferCoeff': 8.3,
      'simTime': 21600,
    }).call();

    return [dummyFunccall];
  }


  @grok.decorators.func({
    outputs: [
      {
        name: 'uploadWidget',
        type: 'widget',
      },
      {
        name: 'uploadFuncCall',
        type: 'funccall',
      },
    ],
  })
  static async CustomUploader(
    @grok.decorators.param({type: 'object'}) params: {func: DG.Func}) {
    const uploadFunc = await grok.functions.eval('Compute:CustomDataUploader') as DG.Func;
    const uploadFuncCall = uploadFunc.prepare({func: params.func});
    const uploadBtn = ui.bigButton('Click me to get mock calls', () => uploadFuncCall.call());

    const dummyWidget = DG.Widget.fromRoot(ui.panel([ui.divV([
      ui.label('This part of dialog comes from my custom data uploader'),
      ui.divH([uploadBtn], {style: {justifyContent: 'center'}}),
    ])]));

    const setLoadingSub = grok.functions.onBeforeRunAction.pipe(
      filter((call) => call.id === uploadFuncCall.id),
    ).subscribe(() => {
      ui.setUpdateIndicator(uploadBtn, true);
    });

    const unsetLoadingSub = grok.functions.onAfterRunAction.pipe(
      filter((call) => call.id === uploadFuncCall.id),
    ).subscribe(() => {
      ui.setUpdateIndicator(uploadBtn, false);
    });

    dummyWidget.subs.push(setLoadingSub, unsetLoadingSub);

    return {uploadWidget: dummyWidget, uploadFuncCall};
  }


  @grok.decorators.func()
  static CustomCustomizer(
    @grok.decorators.param({type: 'object'}) params: {defaultView: DG.TableView}) {
    const comparisonView = params.defaultView;
    comparisonView.scatterPlot({
      'xColumnName': 'Initial temperature',
      'yColumnName': 'Time to cool',
    });
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'validator'}]})
  static SimTimeValidator(
    @grok.decorators.param({type: 'object'}) params: any) {
    const {reasonableMin, reasonableMax} = params;
    return (val: number) => {
      return makeValidationResultInst({
        warnings: val < reasonableMin || val > reasonableMax ?
          [`Minimum reasonable time is ${reasonableMin}. Maximum reasonable time is ${reasonableMax}`]: undefined,
        errors: val < 0 ? [`Time should be strictly positive`]: undefined,
      });
    };
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'validator'}]})
  static DesiredTempValidator(
    @grok.decorators.param({type: 'object'}) params: any) {
    return (val: number, info: ValidationInfo) => {
      const ambTemp = info.funcCall.inputs['ambTemp'];
      const initTemp = info.funcCall.inputs['initTemp'];
      return makeValidationResultInst({
        errors: [
          ...(val < ambTemp) ?
            [makeAdviceInst(`Desired temperature cannot be less than ambient temperature (${ambTemp}). \n`, [
              {actionName: 'Set desired equal to ambient', action: () => info.funcCall.inputs['desiredTemp'] = ambTemp},
            ])]: [],
          ...(val > initTemp) ? [`Desired temperature cannot be higher than initial temperature (${initTemp})`]: [],
        ],
      });
    };
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'validator'}]})
  static InitialTempValidator(
    @grok.decorators.param({type: 'object'}) params: any) {
    return (val: number, info: ValidationInfo) => {
      const ambTemp = info.funcCall.inputs['ambTemp'];
      return makeValidationResultInst({
        errors: [
          ...(val < ambTemp) ? [`Initial temperature cannot be less than ambient temperature (${ambTemp}).`]: [],
        ],
      });
    };
  }

  @grok.decorators.func({outputs: [{type: 'object', name: 'validator'}]})
  static AmbTempValidator(
    @grok.decorators.param({type: 'object'}) params: any) {
    return (val: number, info: ValidationInfo) => {
      const initTemp = info.funcCall.inputs['initTemp'];
      return makeValidationResultInst({
        errors: [
          ...(val > initTemp) ? [`Ambient temperature cannot be higher than initial temperature (${initTemp})`]: [],
        ],
      });
    };
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'validator'}]})
  static HeatCapValidator(
    @grok.decorators.param({type: 'object'}) params: any) {
    return (val: number, info: ValidationInfo) => {
      return makeValidationResultInst({
        errors: [
          ...val <= 0 ? ['Heat capacity must be greater than zero.']: [],
        ],
        notifications: [
          makeAdviceInst(`Heat capacity is only dependent on the object material.`, [
            {actionName: 'Google it', action: () => {window.open(`http://google.com`);}},
          ]),
        ],
      });
    };
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static CustomStringInput(
    @grok.decorators.param({type: 'object'}) params: any) {
    const defaultInput = ui.input.string('Custom input', {value: ''});
    defaultInput.root.style.backgroundColor = 'aqua';
    defaultInput.input.style.backgroundColor = 'aqua';
    return defaultInput;
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static ObjectCoolingSelector(
    @grok.decorators.param({type: 'object'}) params: any) {
    return UiUtils.historyInputJSON(
      'Previous run',
      'ObjectCooling',
    );
  }


}
