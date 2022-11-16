/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import {_functionParametersGrid} from './function-parameters-grid';
import {ModelCatalogView} from './model-catalog-view';
import wu from 'wu';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';
//import {ModelsWidget} from './models-widget'
import {delay} from "@datagrok-libraries/utils/src/test";
import {ComputationView} from "@datagrok-libraries/utils/src/computation-view";
import {FunctionView} from "@datagrok-libraries/utils/src/function-view";
import './css/model-card.css';

let initCompleted: boolean = false;
export const _package = new DG.Package();

//name: OutliersSelectionViewer
//description: Creates an outliers selection viewer
//tags: viewer
//output: viewer
export function OutliersSelection() {
  return new OutliersSelectionViewer();
}

//name: ComputationView
//tags: model-editor, editor
//input: funccall call
//output: view result
export function ComputationViewEditor(call: DG.FuncCall) {
  let v = new ComputationView();
  v.linkFunccall(call);
  v.parentCall = call;
  v.name = call.func.friendlyName;
  v.parentView = grok.functions.getCurrentCall().parentCall.aux['view']; // modelCatalog view
  console.log('comp view editor called')
  //todo: parse url, set parameters to call
  return v;
}

//name: FunctionView
//tags: model-editor, editor
//input: funccall call
//output: view result
export function FunctionViewEditor(call: DG.FuncCall) {
  return new FunctionView(call);
}


/*//output: widget result
//tags: dashboard
export function modelsWidget(): DG.Widget {
  return new ModelsWidget();
}*/

/* eslint-disable */

//description: A spreadsheet that lets you interactively edit parameters and evaluate functions
//tags: functionAnalysis
//input: func f
//output: view result
export function functionParametersGrid(f: DG.Func): DG.View {
  return _functionParametersGrid(f);
}

//name: hof
//description: some description
//sidebar: @compute
export function hof() {
  console.log('hof');
  let f: DG.Func = DG.Func.byName('Sin');
  let v: DG.View = functionParametersGrid(f);
  v.parentCall = grok.functions.getCurrentCall();
  //v.parentView = v.parentCall?.aux['view'];
  grok.shell.addView(v);
}


//name: hof2
//tags: 
//description: some description 2 2 2
//sidebar: @compute
//meta.icon: package1.png
export function hof2() {
  console.log('hof2');

  let f: DG.Func = DG.Func.byName('Sin');
  let v: DG.View = functionParametersGrid(f);

  v.root.appendChild(ui.narrowForm([
    ui.switchInput('test', true),
    ui.boolInput('test', true)]));
  v.root.appendChild(ui.form([
    ui.switchInput('test', true),
    ui.boolInput('test', true)]));
  v.parentCall = grok.functions.getCurrentCall(); // hof2 call itself
  v.parentView = v.parentCall.parentCall?.aux['view']; // modelCatalog view
  let path = v.parentCall.parentCall?.aux['url']; // uri if called from model catalog
  // grok.shell.info(path);


  let path2 = v.parentCall?.aux['url']; // uri if called directly (if app)
  // grok.shell.info(path2);

  v.basePath = '/' + v.parentCall.func.name;
  v.path = '/';

  grok.shell.addView(v);
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
/*
  let modelsList = ui.waitBox(async () => {
    let models = await grok.dapi.scripts
      .filter('#model')
      .list();
    let list = ui.divV(models.map((model) => ui.render(model, {onClick: (_) => ModelHandler.openModel(model)})), {style: {lineHeight: '165%'}});

    let props = ['domain', 'modality'];
    let mtree: { model: DG.Func}[] = models.map((m) => { return {model: m}});
    mtree.forEach((m: {model: DG.Func}) => {
      props.forEach((k) => {
        (<any>m)[k] = m.model.options[k];
      });
    });

    let tree = DG.TreeViewNode.fromItemCategories(mtree,
      props,
      { itemToElement: (x) => ui.render(x.model, {onClick: (_) => ModelHandler.openModel(x.model)}), itemToValue: (x) => x.model, removeEmpty: true }).root

    return ui.tabControl({
      'LIST': list,
      'TREE': tree
    }).root;
  });*/

  initCompleted = true;
}

//name: Model Catalog
//tags: app
//sidebar: @compute
export function modelCatalog() {
  let modelsView = ModelHandler.findModelCatalogView();
  if (modelsView == null) {
    let view = new ModelCatalogView();
    view.name = 'ModelHub';
    let parser = document.createElement('a');
    parser.href = window.location.href;
    let pathSegments = parser.pathname.split('/');
    grok.shell.addView(view);
    // console.log(parser.href);
    if (pathSegments.length > 3) {
      let c = grok.functions.getCurrentCall();
      grok.dapi.functions.filter(`shortName = "${pathSegments[3]}" and #model`).list().then((lst) => {
        if (lst.length == 1)
          ModelHandler.openModel(lst[0], c);
      });
    }
  } else grok.shell.v = modelsView;
}


//name: computationTest
//input: int delayMs {description: Wait time; units: ms}
//input: string error {description: When specified, throws this error}
//output: dataframe result
export async function computationTest(delayMs: number, error: string): Promise<DG.DataFrame> {
  await delay(delayMs);
  if (error != null || error != '')
    throw error;
  return grok.data.demo.demog();
}

//name: testComputationView();
export function testComputationView() {
  let f = DG.Func.find({name: 'computationTest'})[0];
  //grok.shell.addView(new ComputationView(f));
}
