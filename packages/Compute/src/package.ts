/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import {_functionParametersGrid} from './function-views/function-parameters-grid';
import {ModelCatalogView} from './model-catalog-view';
import wu from 'wu';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';
//import {ModelsWidget} from './models-widget'
import {FunctionView} from "@datagrok-libraries/utils/src/function-view";
import {delay} from "@datagrok-libraries/utils/src/test";
import {ComputationView} from "@datagrok-libraries/utils/src/computation-view";
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

//name: ComputationViewTest
//description: Creates an outliers selection viewer
//tags: viewer
//input: funccall call
//output: view result
export function ComputationViewTest(call: DG.FuncCall) {
 // return new ComputationView(call.func);
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
  document.getElementById('root')?.appendChild(ui.h1('Hello World'));
}

//name: hof
//description: some description
//sidebar: @compute
export function hof() {
  grok.shell.info('hof');
  let f: DG.Func = DG.Func.byName('Sin');
  let v: DG.View = functionParametersGrid(f);
  v.parentCall = grok.functions.getCurrentCall();
  //v.parentView = v.parentCall?.aux['view'];
  grok.shell.addView(v);
}


//name: hof2
//description: some description 2 2 2
//sidebar: @compute
//meta.icon: package1.png
export function hof2() {
  grok.shell.info('hof2');
  let f: DG.Func = DG.Func.byName('Sin');
  let v: DG.View = functionParametersGrid(f);

  v.parentCall = grok.functions.getCurrentCall();
  v.parentView = v.parentCall.parentCall?.aux['view'];
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

//tags: init, autostart
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
  });

/*  grok.events.onViewAdding.subscribe((v: DG.ViewBase) => {
    if (v instanceof ComputationView && (v as ComputationView).func?.hasTag("model")) {
      let modelsView = wu(grok.shell.views).find((v) => v.parentCall?.func.name == 'modelCatalog');
      if (modelsView != undefined) {
        v.parentCall = modelsView!.parentCall;
        v.parentView = modelsView!;
        v.toolbox = modelsList;
      }
    }
  });*/

  initCompleted = true;
}

//name: Model Catalog
//tags: app
export function modelCatalog() {
  let modelsView = wu(grok.shell.views).find((v) => v.parentCall?.func.name == 'modelCatalog');
  if (modelsView == null) {
    let view = new ModelCatalogView();
    view.name = 'Models';
    let parser = document.createElement('a');
    parser.href = window.location.href;
    let pathSegments = parser.pathname.split('/');
    grok.shell.addView(view);
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
