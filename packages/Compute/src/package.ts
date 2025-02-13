/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {take} from 'rxjs/operators';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';
export {
  makeValidationResult as makeValidationResult2,
  makeAdvice as makeAdvice2,
  mergeValidationResults as mergeValidationResults2,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {ModelCatalogView, ModelHandler, makeModelCatalog, makeModelTreeBrowser} from '@datagrok-libraries/compute-utils/model-catalog';

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

//name: renderRestPanel
//input: func func
//output: widget panel
export async function renderRestPanel(func: DG.Func): Promise<DG.Widget> {
  const params: object = {};
  func.inputs.forEach((i) => (<any>params)[i.name] = null);
  const curl = `
curl --location --request POST '${(<any>grok.settings).apiUrl}/v1/func/${func.nqName}/run' \\
--header 'Authorization: ${getCookie('auth')}' \\
--header 'Content-Type: application/json' \\
--data-raw '${JSON.stringify(params)}'`;
  const js = `
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
  .catch(error => console.log('error', error));`;
  const tabs = ui.tabControl({'CURL': ui.div([ui.divText(curl)]), 'JS': ui.div([ui.divText(js)])});
  return DG.Widget.fromRoot(tabs.root);
}

function getCookie(name: string): string | undefined {
  const matches = document.cookie.match(new RegExp(
    '(?:^|; )' + name.replace(/([\.$?*|{}\(\)\[\]\\\/\+^])/g, '\\$1') + '=([^;]*)',
  ));
  return matches ? decodeURIComponent(matches[1]) : undefined;
}

//tags: init
export function init() {
  if (initCompleted)
    return;

  if (!(DG.ObjectHandler.list().find((handler) => handler.type === 'Model')))
    DG.ObjectHandler.register(new ModelHandler());


  grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
    const ent = acc.context;
    if (ent == null)
      return;
    if (ent.type != 'script')
      return;
    const restPane = acc.getPane('REST');
    if (!restPane)
      acc.addPane('REST', () => ui.wait(async () => (await renderRestPanel(ent)).root));
  });

  initCompleted = true;
}

let startUriLoaded = false;

//name: Model Catalog
//tags: app
//output: view v
//meta.browsePath: Compute
export function modelCatalog() {
  return makeModelCatalog({
    _package,
    ViewClass: ModelCatalogView,
    HandlerCass: ModelHandler,
    segment: 'Compute',
    viewName: 'Model Catalog',
    funcName: 'modelCatalog',
    startUriLoaded,
  }, () => startUriLoaded = true);
}

//input: dynamic treeNode
//input: view browseView
export async function modelCatalogTreeBrowser(treeNode: DG.TreeViewGroup, browseView: DG.BrowseView) {
  await makeModelTreeBrowser(treeNode, browseView);
}
