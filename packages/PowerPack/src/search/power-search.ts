/* eslint-disable max-len */
/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {widgetHost} from '../utils';
import {initTemplates, templatesSearch} from './templates-search';
import {matchAndParseQuery, processPowerSearchTableView}
  from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';
// Power Search: community-curated, template-based, widget-driven search engine

const widgetFunctions = DG.Func.find({returnType: 'widget'});
const searchFunctions = DG.Func.find({tags: ['search'], returnType: 'list'});
const searchWidgetFunctions = DG.Func.find({tags: ['search'], returnType: 'widget'});
const tableQueriesSearchFunctions = DG.Func.find({meta: {searchPattern: null}, returnType: 'dataframe'})
  .filter((f) => f.options['searchPattern']);

export function initSearch() {
  //grok.dapi.queries.list().then((qs) => queries = qs);
  initTemplates();
}

export function powerSearch(s: string, host: HTMLDivElement): void {
  ui.empty(host);
  tableQueriesFunctionsSearch(s, host);
  jsEvalSearch(s, host) ||
  viewsSearch(s, host);

  regexEntitiesSearch(s, host);
  projectsSearch(s, host);
  functionEvaluationsSearch(s, host);
  searchFunctionsSearch(s, host);
  /*  queriesEntitiesSearch(s, host);
  queriesSearch(s, host);*/
  templatesSearch(s, host);
  widgetsSearch(s, host);
  specificWidgetsSearch(s, host);
}


/// Searches for default views that exactly match the search string
/// Example: functions
function viewsSearch(s: string, host: HTMLDivElement): boolean {
  if (DG.View.ALL_VIEW_TYPES.includes(s)) {
    try {
      host.appendChild(DG.View.createByType(s).root);
      return true;
    } catch (x) {
      return false;
    }
  }
  return false;
}

function regexEntitiesSearch(s: string, host: HTMLDivElement): void {
  const ids = s.split(/,\s*/);
  const semValues = ids.map((id) => DG.SemanticValue.parse(id));
  if (semValues.length < 1 || semValues[0]?.semType == null)
    return;

  if (semValues.every((sv) => sv?.semType && sv.semType == semValues[0].semType) &&
    DG.ObjectHandler.forEntity(semValues[0])) {
    const itemsPanel = ui.div([], {style: {display: 'flex', flexWrap: 'wrap'}});
    const panel = ui.divV([
      ui.span([
        `${semValues.length} ${semValues[0].semType} ${semValues.length == 1 ? 'object. ' : 'objects. '}`,
        ui.link('Open table', () => {
          const df = DG.DataFrame.create(semValues.length);
          df.columns.addNewString(semValues[0].semType).init((i) => semValues[i].value);
          setTimeout(() => {
            grok.shell.addTableView(df);
          });
        })], {style: {marginBottom: '10px'}}),
      itemsPanel,
    ]);
    for (const sv of semValues)
      itemsPanel.append(DG.ObjectHandler.forEntity(sv)!.renderCard(sv));
    host.append(panel);
  }
}


function projectsSearch(s: string, host: HTMLDivElement): void {
  if (s.length < 3)
    return;
  grok.dapi.projects.filter(s).list({pageSize: 5}).then((projects) => {
    if (projects.length > 0) {
      host.appendChild(
        ui.divV([
          ui.h3('Projects'),
          ui.divH(projects.map((p) => ui.renderCard(p))),
        ]),
      );
    }
  });
}

/// Evaluates some of the JavaScript code
/// Example: ui.label('produced by JavaScript')
function jsEvalSearch(s: string, host: HTMLDivElement): boolean {
  if (s.startsWith('grok.') || s.startsWith('DG.') || s.startsWith('ui.')) {
    try {
      const x = eval(s);
      host.appendChild(ui.render(x));
      return true;
    } catch (_) {}
  }
  return false;
}

function capitalizeFirstLetter(string: string): string {
  return string.charAt(0).toUpperCase() + string.slice(1);
}


let listSearchDockedView: DG.View | null = null;
/// Evaluates custom search functions, specifically those that return a list of items
function searchFunctionsSearch(s: string, host: HTMLDivElement): void {
  const listsHost = ui.divH([], {style: {flexWrap: 'wrap', width: '100%'}});
  host.appendChild(listsHost);
  for (const sf of searchFunctions) {
    sf.apply({s: s}).then((results: any[]) => {
      if (results.length > 0) {
        const header = ui.h3(sf.description ?? sf.name);
        const viewName = sf.options['relatedViewName'];
        if (viewName && DG.View.ALL_VIEW_TYPES.includes(viewName)) {
          header.classList.add('power-pack-search-list-has-preview');
          ui.tooltip.bind(header, `Click to open ${sf.options['relatedViewName']} view`);
          header.addEventListener('click', () => {
            if (listSearchDockedView) {
              try {
                grok.shell.dockManager.close(listSearchDockedView.root);
              } catch (e) {
                console.error('Error closing docked view:', e);
              }
              listSearchDockedView = null;
            }
            listSearchDockedView = DG.View.createByType(viewName);
            listSearchDockedView.name = sf.name;
            const node = grok.shell.dockManager.dock(listSearchDockedView.root, DG.DOCK_TYPE.DOWN, grok.shell.dockManager.findNode(grok.shell.v.root), capitalizeFirstLetter(viewName));
            setTimeout(() => {
              const inputElem: HTMLInputElement | null = node.container.containerElement.querySelector('input[type="text"]');
              if (inputElem) {
                inputElem.value = s;
                inputElem.dispatchEvent(new Event('input', {bubbles: true})); // Trigger input event to update the view
              }
            }, 300);
          });
        }

        header.classList.add('power-pack-search-list-header');
        const list = ui.list(results);
        list.classList.add('power-pack-search-list');

        listsHost.appendChild(ui.divV([
          header, list,
        ], {style: {width: 'fit-content'}, classes: 'power-pack-search-list-container'}));
      }
    });
  }
}


function tableQueriesFunctionsSearch(s: string, host: HTMLDivElement): void {
  for (const sf of tableQueriesSearchFunctions) {
    const matchStrings: string[] = (sf.options['searchPattern'] ?? '').split(',').map((s: string) => s.trim())
      .map((s: string) => {
        let ms = s;
        if (s.startsWith('"') || s.startsWith('\''))
          ms = ms.substring(1);
        if (s.endsWith('"') || s.endsWith('\''))
          ms = ms.substring(0, ms.length - 1);
        return ms;
      });
    const inputNames = sf.inputs.map((i) => i.name);
    for (const matchString of matchStrings) {
      const matches = matchAndParseQuery(matchString, s);
      if (!matches)
        continue;
      if (inputNames.length !== Object.entries(matches).length)
        continue;
      const inputParams: any = {};
      Object.entries(matches).forEach(([key, value]) => {
        inputParams[key] = value;
      });

      host.appendChild(widgetHost(DG.Widget.fromRoot(ui.wait(async () => {
        try {
          const fc = sf.prepare(inputParams);
          const resFuncCall = await fc.call();
          const v = resFuncCall.getResultViews();
          if (!v || !v[0] || v[0].type !== DG.VIEW_TYPE.TABLE_VIEW)
            return ui.divText('No result view produced');
          const tv = v[0] as DG.TableView;
          if ((tv.dataFrame?.rowCount ?? 0) === 0)
            return ui.divText('No results found');
          // comma separated list of values, with quotes if its string and without if its number
          const argumentList =
            sf.inputs.map((i) => i.propertyType === DG.TYPE.STRING ? `"${inputParams[i.name]}"` : inputParams[i.name]).join(',');
          tv.dataFrame.setTag(DG.Tags.CreationScript, `Result = ${sf.nqName}(${argumentList})`);
          setTimeout(() => {
            processPowerSearchTableView(tv);
            tv._onAdded();
          }, 200);
          return tv.root;
        } catch (e) {
          console.error(e);
          return ui.divText('Operation caused exception');
        }
      }))),
      );
    }
  }
}

/// Special widgets
function widgetsSearch(s: string, host: HTMLDivElement): void {
  for (const sf of searchWidgetFunctions) {
    const inputName = sf.inputs[0]?.name;
    if (!inputName)
      continue;
    sf.apply({[inputName]: s})
      .then((result: DG.Widget) => {
        if (result)
          host.appendChild(widgetHost(result));
      });
  }
}

/// Explicitly spelled widgets
/// Example: "kpiWidget"
function specificWidgetsSearch(s: string, host: HTMLDivElement): void {
  s = s.toLowerCase();
  for (const wf of widgetFunctions) {
    if (wf.name.toLowerCase() == s) {
      wf.apply().then((w: DG.Widget) => {
        if (w)
          host.appendChild(ui.div([widgetHost(w)]));
      });
    }
  }
}

export async function tableQueriesSearch() {


}

function functionEvaluationsSearch(s: string, host: HTMLDivElement): void {
  if (!s.includes('(') || !s.includes(')'))
    return;

  grok.functions
    .eval(s)
    .then((result) => host.appendChild(ui.span([s + ' = ', result], {style: {fontSize: '20px'}})))
    .catch(() => {});
}
