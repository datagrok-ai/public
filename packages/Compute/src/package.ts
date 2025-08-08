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
} from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {ModelCatalogView, ModelHandler, startModelCatalog, makeModelTreeBrowser, renderRestPanel, setModelCatalogEventHandlers, setModelCatalogHandler} from '@datagrok-libraries/compute-utils/model-catalog';
import {
  testPipeline as testPipelineInst,
} from '@datagrok-libraries/compute-utils/shared-utils/function-views-testing';

import {FittingView} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';

export const _package = new DG.Package();

//name: openModelFromFuncall
//input: funccall funccall
export function openModelFromFuncall(funccall: DG.FuncCall) {
  ModelHandler.openModelFromFunccall(funccall);
}

//name: OutliersSelectionViewer
//description: Creates an outliers selection viewer
//tags: viewer
//output: viewer result
export function OutliersSelection() {
  return new OutliersSelectionViewer();
}

//name: RichFunctionViewEditor
//tags: editor
//input: funccall call
//output: view result
export function RichFunctionViewEditor(call: DG.FuncCall) {
  return RichFunctionViewInst.fromFuncCall(call, {historyEnabled: true, isTabbed: false});
}

//name: PipelineStepEditor
//tags: editor
//input: funccall call
//output: view result
export function PipelineStepEditor(call: DG.FuncCall) {
  return RichFunctionViewInst.fromFuncCall(call, {historyEnabled: false, isTabbed: true});
}

//name: renderRestPanel
//input: func func
//output: widget panel
export async function renderPanel(func: DG.Func): Promise<DG.Widget> {
  return renderRestPanel(func);
}

let startUriLoaded = false;
let initCompleted = false;

const options = {
  _package,
  ViewClass: ModelCatalogView,
  segment: 'Modelhub',
  viewName: 'Model Hub',
  funcName: 'modelCatalog',
  setStartUriLoaded: () => startUriLoaded = true,
  getStartUriLoaded: () => startUriLoaded,
};

//tags: init
export function init() {
  if (initCompleted)
    return;

  setModelCatalogHandler();
  setModelCatalogEventHandlers(options);

  initCompleted = true;
}

//name: Model Hub
//tags: app
//output: view v
//meta.browsePath: Compute
export function modelCatalog() {
  return startModelCatalog(options);
}

//input: dynamic treeNode
//input: view browseView
export function modelCatalogTreeBrowser(treeNode: DG.TreeViewGroup) {
  makeModelTreeBrowser(treeNode);
}

////
// Compute-utils API section
///

export const testPipeline = testPipelineInst;
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

// Compatibility for older compute-api
export const makeAdvice2 = makeAdviceInst;
export const makeValidationResult2 = makeValidationResultInst;


////
// For testing only
///


//name: CustomDataUploader
//input: func func
//output: object uploadedCalls
export async function CustomDataUploader(func: DG.Func) {
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

//name: CustomUploader
//input: object params
//output: widget uploadWidget
//output: funccall uploadFuncCall
export async function CustomUploader(params: {func: DG.Func}) {
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

//name: CustomCustomizer
//input: object params
export function CustomCustomizer(params: {defaultView: DG.TableView}) {
  const comparisonView = params.defaultView;
  comparisonView.scatterPlot({
    'xColumnName': 'Initial temperature',
    'yColumnName': 'Time to cool',
  });
}

//name: SimTimeValidator
//input: object params
//output: object validator
export function SimTimeValidator(params: any) {
  const {reasonableMin, reasonableMax} = params;
  return (val: number) => {
    return makeValidationResultInst({
      warnings: val < reasonableMin || val > reasonableMax ? [`Minimum reasonable time is ${reasonableMin}. Maximum reasonable time is ${reasonableMax}`]: undefined,
      errors: val < 0 ? [`Time should be strictly positive`]: undefined,
    });
  };
}

//name: DesiredTempValidator
//input: object params
//output: object validator
export function DesiredTempValidator(params: any) {
  return (val: number, info: ValidationInfo) => {
    const ambTemp = info.funcCall.inputs['ambTemp'];
    const initTemp = info.funcCall.inputs['initTemp'];
    return makeValidationResultInst({
      errors: [
        ...(val < ambTemp) ? [makeAdviceInst(`Desired temperature cannot be less than ambient temperature (${ambTemp}). \n`, [
          {actionName: 'Set desired equal to ambient', action: () => info.funcCall.inputs['desiredTemp'] = ambTemp},
        ])]: [],
        ...(val > initTemp) ? [`Desired temperature cannot be higher than initial temperature (${initTemp})`]: [],
      ],
    });
  };
}

//name: InitialTempValidator
//input: object params
//output: object validator
export function InitialTempValidator(params: any) {
  return (val: number, info: ValidationInfo) => {
    const ambTemp = info.funcCall.inputs['ambTemp'];
    return makeValidationResultInst({
      errors: [
        ...(val < ambTemp) ? [`Initial temperature cannot be less than ambient temperature (${ambTemp}).`]: [],
      ],
    });
  };
}

//name: AmbTempValidator
//input: object params
//output: object validator
export function AmbTempValidator(params: any) {
  return (val: number, info: ValidationInfo) => {
    const initTemp = info.funcCall.inputs['initTemp'];
    return makeValidationResultInst({
      errors: [
        ...(val > initTemp) ? [`Ambient temperature cannot be higher than initial temperature (${initTemp})`]: [],
      ],
    });
  };
}

//name: HeatCapValidator
//input: object params
//output: object validator
export function HeatCapValidator(params: any) {
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

//name: CustomStringInput
//input: object params
//output: object input
export function CustomStringInput(params: any) {
  const defaultInput = ui.input.string('Custom input', {value: ''});
  defaultInput.root.style.backgroundColor = 'aqua';
  defaultInput.input.style.backgroundColor = 'aqua';
  return defaultInput;
}

//name: ObjectCoolingSelector
//input: object params
//output: object input
export function ObjectCoolingSelector(params: any) {
  return UiUtils.historyInputJSON(
    'Previous run',
    'ObjectCooling',
  );
}

//name: fitTestFunc
//description: Test for optimization: multiple scalars output
//input: double x1 = 1 {caption: param1; min: -3; max: 3}
//input: double x2 = -1 {caption: param2; min: -3; max: 3}
//input: dataframe y {caption: table}
//input: bool bool
//output: int integer
//output: double float1
//output: double float2
//output: dataframe table1 {viewer: Line chart(block: 60) | Grid(block: 40) }
//output: dataframe table2 {viewer: Line chart(block: 60) | Grid(block: 40) }
//editor: Compute:RichFunctionViewEditor
//meta.features: {"fitting": true, "sens-analysis": true}
//meta.runOnOpen: true
//meta.runOnInput: true
export function fitTestFunc(x1: number, x2: number, y: DG.DataFrame, bool: boolean) {
  return {
    integer: x1**3 * (x1 - 1) * x2**3 * (x2 - 1),
    float1: (x2 - 1)**2 + (x1 - 1)**2,
    float2: (x2 - 1)**4 + (x1 - 1)**4,
    table1: DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [x1 + x2 + 1, 4, 9, 16, 25]),
    ]),
    table2: DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [x1 + x2 + 1, 8, 27, 64, 125]),
    ]),
  };
}

//name: testFittingOutputs
//description: Test for optimization: multiple scalars output
export async function testFittingOutputs() {
  const func = await grok.functions.find('Compute:fitTestFunc');

  if (func === null) {
    grok.shell.error('The function "Compute:fitTestFunc" not found!');
    return;
  }

  const targetDf1 = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [1, 4.3, 9.1, 16, 25]),
  ]);
  targetDf1.name = 'test-df1';

  const targetDf2 = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [1, 8.5, 27.6, 64.9, 125]),
  ]);
  targetDf2.name = 'test-df2';

  await FittingView.fromEmpty(func, {
    targets: {
      integer: {default: 123, enabled: true},
      float1: {default: 456.789, enabled: true},
      table1: {
        default: targetDf1,
        enabled: true,
        argumentCol: 'arg',
      },
      table2: {
        default: targetDf2,
        enabled: true,
        argumentCol: 'arg',
      },
    },
  });
}
