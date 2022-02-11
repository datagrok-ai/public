/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import {exportFuncCall} from './export-funccall';
import {_functionParametersGrid} from './function-views/function-parameters-grid';
import {ModelCatalogView} from './model-catalog-view';
import wu from 'wu';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';

let initCompleted: boolean = false;
export const _package = new DG.Package();

//name: OutliersSelectionViewer
//description: Creates an outliers selection viewer
//tags: viewer
//output: viewer
export function OutliersSelection() {
  return new OutliersSelectionViewer();
}

//name: Export to Excel
//input: funccall call
//tags: export
export function exportToExcel(call: DG.FuncCall) {
  exportFuncCall(call);
}

/* eslint-disable */

//description: A spreadsheet that lets you interactively edit parameters and evaluate functions
//tags: functionAnalysis
//input: func f
//output: view result
export function functionParametersGrid(f: DG.Func): DG.View {
  return _functionParametersGrid(f);
}

//name: hof
export function hof() {
  let f: DG.Func = DG.Func.byName('Sin');
  let v: DG.View = functionParametersGrid(f);
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

  grok.events.onViewAdding.subscribe((v: DG.View) => {
    if (v instanceof DG.FunctionView && v.func?.hasTag("model")) {
      let modelsView = wu(grok.shell.views).find((v) => v.parentCall?.func.name == 'modelCatalog');
      if (modelsView != undefined) {
        v.parentCall = modelsView!.parentCall;
        v.parentView = modelsView!;
        v.toolbox = modelsList;
      }
    }
  });

  initCompleted = true;
}

//name: Model Catalog
//tags: app
export function modelCatalog() {
/*  let view = new DG.MultiView({
    viewFactories: {
      'Models': {factory: () => new ModelCatalogView(), allowClose: false}
    }
  });*/
  let view = new ModelCatalogView();
  view.name = 'Models';
  grok.shell.addView(view);
}
