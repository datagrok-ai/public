/* eslint-disable max-len */
/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {welcomeView} from './welcome-view';
import {compareColumns} from './compare-columns';
import {AddNewColumnDialog} from './dialogs/add-new-column';
import {FormulaLinesDialog, DEFAULT_OPTIONS, EditorOptions} from './dialogs/formula-lines';
import {RecentProjectsWidget} from './widgets/recent-projects-widget';
import {CommunityWidget} from './widgets/community-widget';
import {WebWidget} from './widgets/web-widget';
import {appSearch, connectionsSearch,
  dockerSearch, entitySimilaritySearch, filesSearch, functionSearch, groupsSearch,
  helpSearch, jsSamplesSearch, pdbSearch, pubChemSearch, querySearch,
  scriptsSearch, usersSearch, wikiSearch} from './search/entity-search';
import {KpiWidget} from './widgets/kpi-widget';
import {CronInput} from './widgets/cron-input';
import {HtmlWidget} from './widgets/html-widget';
import {viewersDialog} from './viewers-gallery';
import {windowsManagerPanel} from './windows-manager';
import {initSearch, createFuncTableViewWidget} from './search/power-search';
import {newUsersSearch, registerDGUserHandler} from './dg-db';
import {merge} from 'rxjs';
import {HelpObjectHandler} from './search/help-entity';
import {ActivityDashboardWidget} from './widgets/activity-dashboard-widget';
import {DBExplorerEditor} from '@datagrok-libraries/db-explorer/src/editor';
import {setupDBQueryCellHandler, setupGlobalDBExplorer, runEnrichmentFromConfig} from './db-explorer';
export * from './package.g';
export const _package = new DG.Package();
export let _properties: { [propertyName: string]: any };


export class ExcelJSService {
  private _worker: Worker;
  private _isMainWorkerBusy = false;

  private createWorker(): Worker {
    return new Worker(new URL('./workers/exceljs-worker', import.meta.url));
  }

  protected constructor() {
    this._worker = this.createWorker();
  };

  private static _instance: ExcelJSService | null = null;
  public static getInstance() {
    return ExcelJSService._instance ?? (ExcelJSService._instance = new ExcelJSService());
  }

  public parse(bytes: Uint8Array, sheetName?: string): Promise<DG.DataFrame[]> {
    const worker = this._isMainWorkerBusy ? this.createWorker() : this._worker;
    const isUsingMainWorker = this._worker === worker;
    if (isUsingMainWorker)
      this._isMainWorkerBusy = true;
    return new Promise((resolve, reject) => {
      worker.postMessage({action: 'parse', data: bytes, sheetName});
      worker.onmessage = (e) => {
        try {
          if (e.data.error) {
            reject(e.data.error);
            return;
          }
          const results: (string | null)[][][] = e.data.result;
          const names: string[] = e.data.names;
          const dfResults = mapWithFor(results, (result, index) => {
            if (result.length === 0)
              return DG.DataFrame.fromCsv('');

            const namesUsed = new Set<string>();
            const headers = mapWithFor(result[0], (name, i) => {
              const nameToUse = namesUsed.has(name ?? `col${i + 1}`) ? name + ' (1)' : name ?? `col${i + 1}`;
              namesUsed.add(nameToUse);
              return nameToUse;
            });
            const columnCount = headers.length;
            const rowCount = result.length - 1;

            const columnData: string[][] = Array.from({length: columnCount}, () => new Array<string>(rowCount));

            for (let i = 1; i < result.length; i++) {
              const row = result[i];
              for (let j = 0; j < columnCount; j++)
                columnData[j][i - 1] = row[j] ?? '';
            }

            const columns: DG.Column[] = mapWithFor(headers, (name: string, i: number) =>
              DG.Column.fromStrings(name, columnData[i]));
            const df = DG.DataFrame.fromColumns(columns);
            df.name = names[index] || `Sheet ${index + 1}`;
            return df;
          });
          resolve(dfResults);
        } catch (e) {
          console.error(e);
          reject(e);
        } finally {
          if (!isUsingMainWorker)
            worker.terminate();
          else
            this._isMainWorkerBusy = false;
        }
      };
    });
  }

  public terminate(): void {
    this._worker.terminate();
  }
}

export function mapWithFor<T, K>(arr: ArrayLike<T>, mapFunc: (value: T, index: number) => K): K[] {
  const result: K[] = new Array<K>(arr.length);
  for (let i = 0; i < arr.length; i++)
    result[i] = mapFunc(arr[i], i);
  return result;
}

export class PackageFunctions {
  @grok.decorators.func({
    name: 'compareColumns',
    'top-menu': 'Data | Compare Columns...',
  })
  static _compareColumns(): void {
    compareColumns();
  }

  @grok.decorators.func({
    meta: {'autostartImmediate': 'true'},
    name: 'welcomeView',
    outputs: [{name: 'home', type: 'view'}],
  })
  static _welcomeView(): DG.View | undefined {
    return welcomeView();
  }

  @grok.decorators.dashboard({
    meta: {
      'showName': 'false',
    },
    order: '-1',
    name: 'Activity dashboard',
  })
  static activityDashboardWidget(): DG.Widget {
    return new ActivityDashboardWidget();
  }

  @grok.decorators.dashboard({
    order: '6',
    name: 'Community',
  })
  static communityWidget(): DG.Widget {
    return new CommunityWidget();
  }

  @grok.decorators.func()
  static webWidget(): DG.Widget {
    return new WebWidget();
  }

  @grok.decorators.func()
  static htmlWidget(): DG.Widget {
    return new HtmlWidget();
  }

  @grok.decorators.func()
  static kpiWidget(): DG.Widget {
    return new KpiWidget();
  }

  @grok.decorators.func({
    meta: {
      propertyType: 'string',
      semType: 'cron',
      role: 'valueEditor',
    },
    outputs: [{type: 'object', name: 'result'}],
  })
  static cronInput(): DG.InputBase {
    return new CronInput();
  }

  @grok.decorators.func({})
  static isFormulaColumn(
    col: DG.Column): boolean {
    return !!col.getTag(DG.Tags.Formula);
  }

  @grok.decorators.panel({
    name: 'Formula',
    condition: 'PowerPack:isFormulaColumn(col)',
  })
  static formulaWidget(
    col: DG.Column): DG.Widget {
    const expression = col.getTag(DG.Tags.Formula);
    const table = col.dataFrame;
    const f = DG.Func.byName('AddNewColumn');
    const fc = f.prepare({
      'table': table,
      'expression': expression,
      'name': col.name,
      'type': col.type,
    });
    fc.aux['addColumn'] = false;
    const widget = new DG.Widget(ui.div());
    new AddNewColumnDialog(fc, widget);
    return widget;
  }

  @grok.decorators.func({})
  static getFuncTableViewWidget(func: DG.Func, inputParams: Record<string, any>): DG.Widget {
    return DG.Widget.fromRoot(createFuncTableViewWidget(func, inputParams));
  }

  @grok.decorators.func({
    meta: {role: 'searchProvider'},
  })
  static powerPackSearchProvider(): DG.SearchProvider {
    const providers: DG.SearchProvider = {
      'home': [{
        name: 'Files', description: 'Files Search', options: {relatedViewName: 'files'},
        isApplicable: (s: string) => s.length > 2,
        search: (s: string) => filesSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'Functions', description: 'Functions Search', options: {relatedViewName: 'functions'},
        search: (s: string) => functionSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'Scripts', description: 'Scripts Search', options: {relatedViewName: 'scripts'},
        search: (s: string) => scriptsSearch(s).then((r) => ({priority: 10, results: r})),
      },
      {
        name: 'Similarity Search', description: 'Entity Similarity Search',
        search: (s: string) => entitySimilaritySearch(s).then((r) => ({priority: 10, results: r})),
      },
      {
        name: 'Samples', description: 'API Samples Search',
        search: (s: string) => jsSamplesSearch(s).then((r) => ({priority: 9, results: r})),
      }, {
        name: 'Queries', description: 'Queries Search', options: {relatedViewName: 'queries'},
        search: (s: string) => querySearch(s).then((r) => ({priority: 8, results: r})),
      }, {
        name: 'Users', description: 'Users Search', options: {relatedViewName: 'users'},
        search: (s: string) => usersSearch(s).then((r) => ({priority: 11, results: r})),
      }, {
        name: 'Groups', description: 'Groups Search', options: {relatedViewName: 'groups'},
        search: (s: string) => groupsSearch(s).then((r) => ({priority: 12, results: r})),
      }, {
        name: 'Dockers', description: 'Dockers Search', options: {relatedViewName: 'dockers'},
        search: (s: string) => dockerSearch(s).then((r) => ({priority: 13, results: r})),
      }, {
        name: 'Help', description: 'Help Search',
        search: (s: string) => helpSearch(s).then((r) => ({priority: 7, results: r})),
      }, {
        name: 'PDB', description: 'Protein Data Bank Search',
        getSuggestions: (s) => s?.length < 4 && s?.trim().length > 1 && s.toUpperCase() === s ?
          [{priority: 5, suggestionText: 'PDB ID, e.g. 4AKZ', suggestionValue: '4AKZ'}] : null,
        search: (s: string) => pdbSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'PubChem', description: 'PubChem Search', options: {widgetHeight: 500},
        getSuggestions: (s) => s.length > 1 && 'aspirin'.includes(s) ?
          [{priority: 5, suggestionText: 'Aspirin', suggestionValue: 'aspirin'}] : null,
        search: (s: string) => pubChemSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'Wiki', description: 'Wikipedia Search', options: {widgetHeight: 500},
        getSuggestions: (s) => null,
        search: (s: string) => wikiSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'Apps', description: 'Apps Search', options: {relatedViewName: 'apps'},
        search: (s) => appSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'Connections', description: 'Connections Search', options: {relatedViewName: 'connections'},
        search:
          (s) => connectionsSearch(s).then((r) => ({priority: 10, results: r})),
      }, {
        name: 'New Users', description: 'New Users Search',
        getSuggestions: (s) => s?.toLowerCase().startsWith('new') && !s?.toLowerCase()?.startsWith('new users ') ?
          [{suggestionText: 'New Users Today', priority: 19},
            {suggestionText: 'New users This Month', priority: 20},
            {suggestionText: 'New users This Year', priority: 21},
            {suggestionText: 'New users last 3 months', priority: 22},
            {suggestionText: 'New users yesterday', priority: 24},
            {suggestionText: 'New user last 7 days', priority: 23}] : null,
        search: (s: string) => newUsersSearch(s).then((r) => ({priority: 10, results: r})),
      },

      ],
    };
    return providers;
  }

  @grok.decorators.func({
    name: 'formulaLinesEditor',
    outputs: [],
  })
  static formulaLinesDialog(
    @grok.decorators.param({type: 'dataframe', options: {optional: true}}) src: DG.DataFrame | DG.Viewer,
    @grok.decorators.param({type: 'int', options: {optional: true}}) currentIndexToSet?: number,
    @grok.decorators.param({type: 'bool', options: {optional: true}}) isDataFrameValue?: boolean,
    @grok.decorators.param({type: 'bool', options: {optional: true}}) isAnnotationArea?: boolean,
  ): void {
    const options = Object.keys(_properties)
      .filter((k) => k in DEFAULT_OPTIONS)
      .reduce((opts, k) => (opts[k] = _properties[k], opts), <EditorOptions>{});
    //TODO: use property's 'category' or 'tags' to distinguish relevant properties

    new FormulaLinesDialog(src, options, {
      index: currentIndexToSet,
      isDataFrame: isDataFrameValue,
      isAnnotationArea: isAnnotationArea,
    });
  }

  @grok.decorators.init()
  static async powerPackInit() {
    DG.ObjectHandler.register(new HelpObjectHandler());
    setupGlobalDBExplorer(); // lazy without await
    setupDBQueryCellHandler(); // db-explorer for any query result - lazy without await
    initSearch();

    _properties = await _package.getProperties();
    registerDGUserHandler(); // lazy without await

    // saving and restoring the scrolls when changing views
    const maxDepth = 40;
    const elementMap = new Map<HTMLElement, number>();
    const getScrolledChild = (node: HTMLElement, depth = 0) => {
      if (depth > maxDepth)
        return;
      if (node.scrollTop)
        elementMap.set(node, node.scrollTop);
      node.children && Array.from(node.children)
        .forEach((child) => getScrolledChild(child as HTMLElement, depth + 1));
    };
    const setScrolls = () => {
      elementMap.forEach((scroll, node) => {
        if (document.body.contains(node)) {
          node.scrollTop = scroll;
          elementMap.delete(node);
        }
      });
    };
    merge(grok.events.onCurrentViewChanging, grok.events.onViewChanging)
      .subscribe((_) => getScrolledChild(document.body));
    DG.debounce(merge(grok.events.onCurrentViewChanged, grok.events.onViewChanged), 100).subscribe((_) => setScrolls());
  }

  @grok.decorators.autostart({description: 'Windows Manager'})
  static windowsManager() {
    windowsManagerPanel();
  }

  @grok.decorators.func({description: 'Open \'Viewer Gallery\' dialog'})
  static viewerDialog(
    tv: DG.TableView) {
    if (tv instanceof DG.TableView)
      return viewersDialog(tv, tv.table!);
  }

  @grok.decorators.autostart({description: 'ViewerGallery'})
  static viewerGallery(): void {
    grok.events.onViewAdded.subscribe((view) => _viewerGallery(view));
    _viewerGallery(grok.shell.v);
  }

  @grok.decorators.fileViewer({
    fileViewer: 'md,mdx',
  })
  static async markdownFileViewer(
    file: DG.FileInfo): Promise<DG.View> {
    const viewFile = DG.View.create();
    viewFile.name = file.name.slice(0, file.name.indexOf('.'));
    const mdText = await file.readAsString();
    const preview = await ui.input.markdownPreview(mdText);
    viewFile.append(preview);
    return viewFile;
  }

  @grok.decorators.fileHandler({
    ext: 'xlsx',
    description: 'Opens Excel file',
  })
  static async xlsxFileHandler(
    @grok.decorators.param({'type': 'list'}) bytes: Uint8Array,
    @grok.decorators.param({'options': {'optional': true}}) sheetName?: string): Promise<DG.DataFrame[]> {
    const XLSX_MAX_FILE_SIZE = 80 * 1024 * 1024; // 80 MB
    if (bytes.length > XLSX_MAX_FILE_SIZE)
      throw new Error('The file you are trying to open is too large. Excel max file size is 80MB.');
    const excelJSService = ExcelJSService.getInstance();
    return (await excelJSService.parse(bytes, sheetName));
  }

  @grok.decorators.func({
    meta: {role: 'transform'}
  })
  static async runEnrichment(conn: DG.DataConnection, schema: string, table: string, column: string, name: string, df: DG.DataFrame): Promise<void> {
    return runEnrichmentFromConfig(conn, schema, table, column, name, df);
  }
}

//name: addNewColumn
//input: funccall call {optional: true}
//editor-for: AddNewColumn
export function addNewColumnDialog(call: DG.FuncCall | null = null): AddNewColumnDialog {
  return new AddNewColumnDialog(call as any);
}

grok.events.onContextMenu.subscribe((args) => {
  const src = args.args.context;
  let menu;
  if (src instanceof DG.ScatterPlotViewer ||
      (src instanceof DG.Viewer && src.getOptions()['type'] == DG.VIEWER.LINE_CHART))
    menu = args.args.menu.find('Tools');

  if (src instanceof DG.Viewer && src.getOptions()['type'] == DG.VIEWER.TRELLIS_PLOT &&
      src.getOptions().look['viewerType'] == DG.VIEWER.SCATTER_PLOT)
    menu = args.args.menu.find(DG.VIEWER.SCATTER_PLOT).find('Tools');

  menu?.item('Formula Lines...', () => {
    PackageFunctions.formulaLinesDialog(src);
  });
});

function _viewerGallery(view: DG.ViewBase): void {
  if (view?.type == 'TableView') {
    const panels = view.getRibbonPanels();
    for (const p of panels) {
      for (const d of p) {
        if (d instanceof HTMLDivElement) {
          if (d.querySelector('.svg-add-viewer')) {
            const icon = ui.iconFA('',
              () => {
                viewersDialog(view as DG.TableView, (view as DG.TableView).table!);
              }, 'Add viewer');
            icon.className = 'grok-icon svg-icon svg-add-viewer';
            const btn = ui.div([icon]);
            btn.className = 'd4-ribbon-item';
            p[p.indexOf(d)] = btn;
          }
        }
      }
    }
    view.setRibbonPanels(panels);
  }
}
