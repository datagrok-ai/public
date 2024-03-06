/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';
import {RichFunctionView} from "@datagrok-libraries/compute-utils";
import { ValidationInfo, makeAdvice, makeValidationResult } from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {ModelCatalogView, ModelHandler} from '@datagrok-libraries/compute-utils/model-catalog';
import './css/model-card.css';

let initCompleted: boolean = false;
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
  return RichFunctionView.fromFuncCall(call, {historyEnabled: true, isTabbed: false});
}

//name: PipelineStepEditor
//tags: editor
//input: funccall call
//output: view result
export function PipelineStepEditor(call: DG.FuncCall) {
  return RichFunctionView.fromFuncCall(call, {historyEnabled: false, isTabbed: true});
}

//name: renderRestPanel
//input: func func
//output: widget panel
export async function renderRestPanel(func: DG.Func): Promise<DG.Widget> {
  let params: object = {};
  func.inputs.forEach((i) => (<any>params)[i.name] = null);
let curl = `
curl --location --request POST '${(<any>grok.settings).apiUrl}/v1/func/${func.nqName}/run' \\
--header 'Authorization: ${getCookie('auth')}' \\
--header 'Content-Type: application/json' \\
--data-raw '${JSON.stringify(params)}'`
let js = `
var myHeaders = new Headers();
myHeaders.append("Authorization", "${getCookie('auth')}");
myHeaders.append("Content-Type", "application/json");

var raw = JSON.stringify(${JSON.stringify(params)});

var requestOptions = {
  method: 'POST',
  headers: myHeaders,
  body: raw,
  redirect: 'follow'
};

fetch("${(<any>grok.settings).apiUrl}/v1/func/${func.nqName}/run", requestOptions)
  .then(response => response.text())
  .then(result => console.log(result))
  .catch(error => console.log('error', error));`
  let tabs = ui.tabControl({'CURL': ui.div([ui.divText(curl)]), 'JS': ui.div([ui.divText(js)])})
  return DG.Widget.fromRoot(tabs.root);
}

function getCookie(name: string): string | undefined{
  let matches = document.cookie.match(new RegExp(
    "(?:^|; )" + name.replace(/([\.$?*|{}\(\)\[\]\\\/\+^])/g, '\\$1') + "=([^;]*)"
  ));
  return matches ? decodeURIComponent(matches[1]) : undefined;
}


//tags: autostart
//meta.autostartImmediate: true
export function autostart() {

}

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
    
  return [dummyFunccall]
}

//name: CustomUploader
//input: object params
//output: widget uploadWidget
//output: funccall uploadFuncCall
export async function CustomUploader(params: {func: DG.Func}) {
  const uploadFunc = await grok.functions.eval('Compute:CustomDataUploader') as DG.Func;
  const uploadFuncCall = uploadFunc.prepare({func: params.func})
  const uploadBtn = ui.bigButton('Click me to get mock calls', () => uploadFuncCall.call());

  const dummyWidget = DG.Widget.fromRoot(ui.panel([ui.divV([
    ui.label('This part of dialog comes from my custom data uploader'),
    ui.divH([uploadBtn], {style: {'justify-content': 'center'}})
  ])]));  

  const setLoadingSub = grok.functions.onBeforeRunAction.pipe(
    filter((call) => call.id === uploadFuncCall.id)
  ).subscribe(() => {
    ui.setUpdateIndicator(uploadBtn, true);
  })

  const unsetLoadingSub = grok.functions.onAfterRunAction.pipe(
    filter((call) => call.id === uploadFuncCall.id)
  ).subscribe(() => {
    ui.setUpdateIndicator(uploadBtn, false);
  })

  dummyWidget.subs.push(setLoadingSub, unsetLoadingSub)

  return {uploadWidget: dummyWidget, uploadFuncCall};
}

//name: CustomCustomizer
//input: object params
export function CustomCustomizer(params: {defaultView: DG.TableView}) {
  const comparisonView = params.defaultView;
  comparisonView.scatterPlot({
    "xColumnName": "Initial temperature",
    "yColumnName": "Time to cool",
  });
}

//tags: init
export function init() {
  if (initCompleted)
    return;

  DG.ObjectHandler.register(new ModelHandler());
  
  grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
    const ent = acc.context;
    if (ent == null)
      return;
    if (ent.type != 'script')
      return;
    let restPane = acc.getPane('REST');
    if (!restPane)
      acc.addPane('REST', () => ui.wait(async () => (await renderRestPanel(ent)).root));
  });

  initCompleted = true;
}

let startUriLoaded = false;

//name: Model Catalog
//tags: app
//sidebar: @compute
export function modelCatalog() {
  const view = ModelCatalogView.findOrCreateCatalogView('Model Catalog', 'modelCatalog', _package);
  
  const optionalPart = grok.shell.startUri.indexOf('?');
  const pathSegments = grok.shell.startUri
    .substring('https://'.length, optionalPart > 0 ? optionalPart: undefined)
    .split('/');
  if (!startUriLoaded && pathSegments.length > 3) {
    grok.dapi.functions.filter(`shortName = "${pathSegments[3]}" and #model`).list().then((lst) => {
      if (lst.length == 1)
        ModelHandler.openModel(lst[0]);
        startUriLoaded = true;
    });
  }

  return view;
}


//name: SimTimeValidator
//input: object params
//output: object validator
export function SimTimeValidator(params: any) {
  const {reasonableMin, reasonableMax} = params;
  return (val: number) => {
    return makeValidationResult({
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
    return makeValidationResult({
      errors: [
        ...(val < ambTemp) ? [makeAdvice(`Desired temperature cannot be less than ambient temperature (${ambTemp}). \n`, [
          {actionName: 'Set desired equal to ambient', action: () => info.funcCall.inputs['desiredTemp'] = ambTemp }
        ])]: [],
        ...(val > initTemp) ? [`Desired temperature cannot be higher than initial temperature (${initTemp})`]: [],
      ]
    });
  };
}

//name: InitialTempValidator
//input: object params
//output: object validator
export function InitialTempValidator(params: any) {
  return (val: number, info: ValidationInfo) => {
    const ambTemp = info.funcCall.inputs['ambTemp'];
    return makeValidationResult({
      errors: [
        ...(val < ambTemp) ? [`Initial temperature cannot be less than ambient temperature (${ambTemp}).`]: [],
      ]
    });
  };
}

//name: AmbTempValidator
//input: object params
//output: object validator
export function AmbTempValidator(params: any) {
  return (val: number, info: ValidationInfo) => {
    const initTemp = info.funcCall.inputs['initTemp'];
    return makeValidationResult({
      errors: [
        ...(val > initTemp) ? [`Ambient temperature cannot be higher than initial temperature (${initTemp})`]: [],
      ]
    });
  };
}

//name: HeatCapValidator
//input: object params
//output: object validator
export function HeatCapValidator(params: any) {
  return (val: number, info: ValidationInfo) => {
    return makeValidationResult({
      errors: [
        ...val <= 0 ? ['Heat capacity must be greater than zero.']: []
      ],
      notifications: [
        makeAdvice(`Heat capacity is only dependent on the object material.`, [
          {actionName: 'Google it', action: () => { window.open(`http://google.com`)}}
        ]),
      ]
    });
  };
}

//name: CustomStringInput
//input: object params
//output: object input
export function CustomStringInput(params: any) {
  const defaultInput = ui.stringInput('Custom input', '');
  defaultInput.root.style.backgroundColor = 'aqua';
  defaultInput.input.style.backgroundColor = 'aqua';
  return defaultInput;
}