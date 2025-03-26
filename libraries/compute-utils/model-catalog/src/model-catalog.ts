/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {take} from 'rxjs/operators';
import {ModelCatalogView} from './model-catalog-view';
import {ModelHandler} from './model-handler';

export interface ModelCatalogConfig {
  segment: string,
  viewName: string,
  funcName: string,
  _package: DG.Package,
  ViewClass: typeof ModelCatalogView,
  setStartUriLoaded: () => void,
  getStartUriLoaded: () => boolean,
}

export function setModelCatalogEventHandlers(options: ModelCatalogConfig) {
  const {
    viewName, funcName, _package, ViewClass
  } = options;

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

  grok.functions.onBeforeRunAction.subscribe((fc) => {
    if (fc.func.hasTag('model')) {
      const view = options.ViewClass.findOrCreateCatalogView(viewName, funcName, _package);
      const currentView = [...grok.shell.views].find(v => v === view);
      if (!currentView)
        grok.shell.add(view);
      ViewClass.bindModel(view, fc);
    } else if (fc.inputs?.['call']?.func instanceof DG.Func && fc.inputs['call'].func.hasTag('model')) {
      const view = options.ViewClass.findOrCreateCatalogView(viewName, funcName, _package);
      const currentView = [...grok.shell.views].find(v => v === view);
      if (!currentView)
        grok.shell.add(view);
      ViewClass.bindModel(view, fc.inputs['call']);
    }
  });

  grok.events.onCurrentViewChanged.subscribe(async () => {
    let prevSyncState = true;
    if (grok.shell.v?.name === options.viewName) {
      prevSyncState = grok.shell.windows.help.syncCurrentObject;
      grok.shell.windows.help.syncCurrentObject = false;
    } else
      grok.shell.windows.help.syncCurrentObject = prevSyncState;
  });
}

export function setModelCatalogHandler() {
  if (!(DG.ObjectHandler.list().find((handler) => handler.type === 'Model')))
    DG.ObjectHandler.register(new ModelHandler());
}

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
myHeaders.append("Content-Type", "applicati4on/json");

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

export function startModelCatalog(options: ModelCatalogConfig) {
  const {
    segment, viewName, funcName, _package, ViewClass, setStartUriLoaded, getStartUriLoaded
  } = options;

  if (!getStartUriLoaded()) {
    const view = ViewClass.findOrCreateCatalogView(viewName, funcName, _package);
    handleInitialUri(segment);
    setStartUriLoaded();
    return view;
  }

  const view = ViewClass.findModelCatalogView(viewName);

  if (view) {
    grok.shell.v = view;
    return null;
  } else {
    const newView = ViewClass.createModelCatalogView(viewName, funcName, _package);
    return newView;
  }
}

function handleInitialUri(segment: string) {
  const url = new URL(grok.shell.startUri);
  const startPathSegments = url.pathname.split('/');
  const catlogUriSegmentIdx = startPathSegments.findIndex(x => x === segment);

  if (catlogUriSegmentIdx > 0) {
    if (catlogUriSegmentIdx + 1 < startPathSegments.length) {
      const shortName = startPathSegments[catlogUriSegmentIdx + 1];
      grok.dapi.functions.filter(`shortName = "${shortName}" and #model`).list().then((lst) => {
        if (lst.length == 1)
          ModelHandler.openModel(lst[0]);
      });
    }
  }
}

export async function makeModelTreeBrowser(treeNode: DG.TreeViewGroup) {
  const NO_CATEGORY = 'Uncategorized' as const;

  const modelList = await grok.dapi.functions.filter('(#model)').list();
  const departments = modelList.reduce((acc, model) => {
    if (model.options.department)
      acc.add(model.options.department);
    else
      acc.add(NO_CATEGORY);
    return acc;
  }, new Set([] as string[]));

  const hlProcesses = modelList.reduce((acc, model) => {
    if (model.options.HL_process)
      acc.add(model.options.HL_process);
    else
      acc.add(NO_CATEGORY);
    return acc;
  }, new Set([] as string[]));

  const processes = modelList.reduce((acc, model) => {
    if (model.options.process)
      acc.add(model.options.process);
    else
      acc.add(NO_CATEGORY);
    return acc;
  }, new Set([] as string[]));

  for (const department of departments) {
    const serverDep = (department !== NO_CATEGORY ? department: undefined);
    const hasModelsDep = modelList.find((model) =>
      model.options.department === serverDep,
    );

    if (!hasModelsDep) continue;

    const depNode = treeNode.getOrCreateGroup(department, null, false);
    for (const hlProcess of hlProcesses) {
      const serverHlProcess = (hlProcess !== NO_CATEGORY ? hlProcess: undefined);
      const hasModelsHl = modelList.find((model) =>
        model.options.department === serverDep &&
          model.options.HL_process === serverHlProcess,
      );

      if (!hasModelsHl) continue;

      const hlNode = hlProcess !== NO_CATEGORY ? depNode.getOrCreateGroup(hlProcess, null, false): depNode;
      for (const process of processes) {
        const serverProcess = (process !== NO_CATEGORY ? process: undefined);
        const hasModels = modelList.find((model) =>
          model.options.department === serverDep &&
            model.options.HL_process === serverHlProcess &&
            model.options.process === serverProcess,
        );

        if (!hasModels) continue;

        const processNode = process !== NO_CATEGORY ? hlNode.getOrCreateGroup(process, null, false): hlNode;
        processNode.onNodeExpanding.pipe(take(1)).subscribe(() => {
          const modelRule = `(#model)`;
          const depRules = department === NO_CATEGORY ? `(options.department = null)`: `(options.department in ("${department}"))`;
          const hlProcessRules = hlProcess === NO_CATEGORY ? `(options.HL_process = null)`: `(options.HL_process in ("${hlProcess}"))`;
          const processRules = process === NO_CATEGORY ? `(options.process = null)`: `(options.process in ("${process}"))`;
          const filteringRules = `(${[modelRule, depRules, hlProcessRules, processRules].join(' and ')})`;
          processNode.loadSources(grok.dapi.functions.filter(filteringRules));
        });
      }
    }
  }
}
