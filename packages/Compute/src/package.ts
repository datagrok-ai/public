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
import {ModelCatalogView,
  ModelHandler,
  startModelCatalog,
  makeModelTreeBrowser,
  renderRestPanel,
  setModelCatalogEventHandlers,
  setModelCatalogHandler} from '@datagrok-libraries/compute-utils/model-catalog';
import {
  testPipeline as testPipelineInst,
} from '@datagrok-libraries/compute-utils';
import {
  deepCopy as  deepCopyInst,
} from '@datagrok-libraries/compute-utils';


import {FittingView} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
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


export class PackageFunctions {
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


  @grok.decorators.func({name: 'renderRestPanel'})
  static async renderPanel(@grok.decorators.param({type: 'func'}) func: DG.Func) : Promise<DG.Widget> {
    return renderRestPanel(func as any) as any;
  }

  @grok.decorators.init()
  static init() {
    if (initCompleted)
      return;

    setModelCatalogHandler();
    setModelCatalogEventHandlers(options as any);

    initCompleted = true;
  }


  @grok.decorators.app({
    browsePath: 'Compute',
    name: 'Model Hub',
    outputs: [{type: 'view', name: 'result'}],
  })
  static modelCatalog() {
    return startModelCatalog(options as any);
  }


  @grok.decorators.func({
    meta: { role: ' ', app: ' '}
  })
  static modelCatalogTreeBrowser(treeNode: DG.TreeViewGroup, browseView: DG.ViewBase) {
    makeModelTreeBrowser(treeNode as any);
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


  @grok.decorators.func({
    meta: {
      features: '{"fitting": true, "sens-analysis": true}',
      runOnOpen: 'true',
      runOnInput: 'true',
    },
    outputs: [
      {
        name: 'integer',
        type: 'int',
      },
      {
        name: 'float1',
        type: 'double',
      },
      {
        name: 'float2',
        type: 'double',
      },
      {
        name: 'table1',
        type: 'dataframe',
        options: {viewer: 'Line chart(block:60) | Grid(block:40)'},
      },
      {
        name: 'table2',
        type: 'dataframe',
        options: {viewer: 'Line chart(block:60) | Grid(block:40)'},
      },
    ],
    description: 'Test for optimization: multiple scalars output',
    editor: 'Compute:RichFunctionViewEditor',
  })
  static fitTestFunc(
    @grok.decorators.param({options: {caption: 'param1', min: '-3', max: '3', initialValue: '1'}}) x1: number,
    @grok.decorators.param({options: {caption: 'param2', min: '-3', max: '3', initialValue: '-1'}}) x2: number,
    @grok.decorators.param({options: {caption: 'table'}}) y: DG.DataFrame,
      bool: boolean) {
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


  @grok.decorators.func({description: 'Test for optimization: multiple scalars output'})
  static async testFittingOutputs() {
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

    await FittingView.fromEmpty(func as any, {
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
}
