/* eslint-disable max-len */
/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {widgetHost} from '../utils';
import {initTemplates, templatesSearch} from './templates-search';
import {matchAndParseQuery, processPowerSearchTableView}
  from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';
import {exactAppFuncSearch} from './entity-search';
// Power Search: community-curated, template-based, widget-driven search engine

const widgetFunctions = DG.Func.find({returnType: 'widget'});
const searchFunctions = DG.Func.find({tags: ['search'], returnType: 'list'});
const searchWidgetFunctions = DG.Func.find({tags: ['search'], returnType: 'widget'});
const tableQueriesSearchFunctions = DG.Func.find({meta: {searchPattern: null}, returnType: 'dataframe'})
  .filter((f) => f.options['searchPattern']);

const searchProvidersFunctions = DG.Func.find({meta: {role: DG.FUNC_TYPES.SEARCH_PROVIDER}});

export function initSearch() {
  //grok.dapi.queries.list().then((qs) => queries = qs);
  initTemplates();
}

export function powerSearch(s: string, host: HTMLDivElement, inputElement: HTMLInputElement): void {
  ui.empty(host);
  tableQueriesFunctionsSearch(s, host);
  // tableQueriesFunctionsSearchLlm(s, host);
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
  exactAppViewSearch(s);
  searchProvidersSearch(s, host, inputElement);
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

function processListSearchResults(
  results: any[], host: HTMLDivElement, s: string,
  options: {name: string, description?: string, relatedViewName?: string},
  order?: number) {
  const listsHost = host.querySelector('.power-search-lists-host') ?? ui.divH([], {style: {flexWrap: 'wrap', width: '100%'}, classes: 'power-search-lists-host'});
  host.appendChild(listsHost);
  if (results && results.length > 0) {
    const header = ui.h3(options?.name ?? options?.description);
    const viewName = options.relatedViewName;
    if (viewName && DG.View.ALL_VIEW_TYPES.includes(viewName)) {
      header.classList.add('power-pack-search-list-has-preview');
      ui.tooltip.bind(header, `Click to open ${options.relatedViewName} view`);
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
        listSearchDockedView.name = options.name;
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
    ], {style: {width: 'fit-content', order: `${order ?? 9999}`}, classes: 'power-pack-search-list-container'}));
  }
}

/// Evaluates custom search functions, specifically those that return a list of items
function searchFunctionsSearch(s: string, host: HTMLDivElement): void {
  for (const sf of searchFunctions) {
    sf.apply({s: s}).then((results: any[]) => {
      processListSearchResults(results, host, s, {name: sf.name, description: sf.description, relatedViewName: sf.options['relatedViewName']});
    });
  }
}

const searchProviders: DG.SearchProvider[] = [];
let _searchProvidersPromise: Promise<void> | null = null;
async function _getSearchProviders(): Promise<void> {
  const providers = searchProvidersFunctions.map(async (f) => {
    try {
      return await f.apply({}) as Promise<DG.SearchProvider>;
    } catch (e) {
      console.error(`Error applying search provider function ${f.name}:`, e);
      return null; // Return null if the function fails
    }
  });
  const results = await Promise.all(providers);
  for (const provider of results) {
    if (provider)
      searchProviders.push(provider);
  }
}

let _suggestionsMenu: DG.Menu | null = null;
let _ctrlSpaceSubscribed = false;
async function searchProvidersSearch(s: string, host: HTMLDivElement, searchInput: HTMLInputElement): Promise<void> {
  _suggestionsMenu?.hide();
  _suggestionsMenu = null;
  if (!s?.trim())
    return;
  _searchProvidersPromise ??= _getSearchProviders();
  await _searchProvidersPromise;

  let showSuggestions = false;
  const addSuggestion = (text: string, value?: string, priority?: number, onValueEnter?: (s: string, view?: DG.ViewBase) => void) => {
    // make sure that its not an exact match, where it does not make sense to add it as a suggestion
    if ((value ?? text) === searchInput.value)
      return;
    _suggestionsMenu!.item(text, () => {
      if (value) {
        searchInput.value = value;
        searchInput.dispatchEvent(new Event('input', {bubbles: true})); // Trigger input event to update the view
      }
      searchInput.focus();
      if (onValueEnter)
        onValueEnter(searchInput.value, grok.shell.v);
    }, priority);
    showSuggestions = true;
  };

  // separate handler for suggestions
  const processSuggestions = () => {
    const si = searchInput.value?.trim() ?? '';
    _suggestionsMenu?.hide();
    _suggestionsMenu = null;
    _suggestionsMenu = DG.Menu.popup();
    showSuggestions = false;
    searchProviders.forEach((sp) => {
      if (!sp['home'])
        return;
      const homeSearchProviders = Array.isArray(sp.home) ? sp.home : [sp.home];
      homeSearchProviders.forEach((pp) => {
      // first handle the suggestions
        if (pp.getSuggestions) {
          const suggestions = pp.getSuggestions(si);
          if (suggestions && suggestions.length > 0) {
            suggestions.forEach((sug) => {
              addSuggestion(sug.suggestionText, sug.suggestionValue, sug.priority, pp.onValueEnter);
            });
          }
        }
      });
    });
    // handle automatic suggestions based on objecthandlers or dynamic queries
    DG.ObjectHandler.list()
      .filter((h) => h.regexpExample && h.regexpExample.nonVariablePart?.toLowerCase().includes(si.toLowerCase()))
      .forEach((h) => {
        const e = h.regexpExample!;
        const eg = !e.example || e.example == e.regexpMarkup ? '' : `(e.g. ${e.example})`;
        addSuggestion(`${e.regexpMarkup} ${eg}`, e.nonVariablePart, 0); // 0 order to make sure it is at the top
      });
    tableQueriesSearchFunctions.forEach((f) => {
      const pattern: string = f.options['searchPattern'] ?? '';
      if (pattern.toLowerCase().includes(si.toLowerCase()))
        addSuggestion(removeTrailingQuotes(pattern), removeTrailingQuotes(pattern), 1);
    });
    if (showSuggestions) {
      _suggestionsMenu.show({element: searchInput.parentElement!, x: searchInput.offsetLeft,
        y: searchInput.offsetTop + searchInput.offsetHeight});
    }
  };
  processSuggestions();
  if (!_ctrlSpaceSubscribed) {
    _ctrlSpaceSubscribed = true;
    searchInput.addEventListener('keydown', (event) => {
      if (event.ctrlKey && (event.key === ' ' || event.key === 'Spacebar')) {
        event.preventDefault();
        processSuggestions();
      }
    });
  }


  searchProviders.forEach((sp) => {
    if (!sp['home'])
      return; // Skip providers that do not have a home property
    const homeSearchProviders = Array.isArray(sp.home) ? sp.home : [sp.home];
    homeSearchProviders.forEach((pp) => {
      // handle search
      if (pp.isApplicable && !pp.isApplicable(s))
        return;
      pp.search(s, grok.shell.v).then((result) => {
        if (!result || !result.results)
          return;
        if (Array.isArray(result.results))
          processListSearchResults(result.results, host, s, {name: pp.name, description: pp.description, relatedViewName: pp.options?.relatedViewName}, result.priority);
        if (result.results instanceof DG.Widget) {
          const h = widgetHost(result.results);
          if (pp.options?.widgetHeight)
            h.style.setProperty('height', `${pp.options.widgetHeight}px`, 'important');
          host.appendChild(h);
        }
        // TODO: Add more return types handling
      }).catch((e) => {
        console.error(`Error searching with provider ${pp.name}:`, e);
      });
    });
  });
}

function removeTrailingQuotes(s: string): string {
  let ms = s;
  if (s.startsWith('"') || s.startsWith('\''))
    ms = ms.substring(1);
  if (s.endsWith('"') || s.endsWith('\''))
    ms = ms.substring(0, ms.length - 1);
  return ms;
}

function tableQueriesFunctionsSearch(s: string, host: HTMLDivElement): void {
  for (const sf of tableQueriesSearchFunctions) {
    const matchStrings: string[] = (sf.options['searchPattern'] ?? '').split(',').map((s: string) => s.trim())
      .map((s: string) => {
        return removeTrailingQuotes(s);
      });
    const inputNames = sf.inputs.map((i) => i.name);
    for (const matchString of matchStrings) {
      const matches = matchAndParseQuery(matchString, s);
      if (!matches)
        continue;
      if (inputNames.length !== Object.entries(matches).length)
        continue;
      const inputParams: Record<string, any> = {};
      if (Object.values(matches).some((v) => !v || v.startsWith('${') || v.endsWith('}')))
        continue;
      Object.entries(matches).forEach(([key, value]) => {
        inputParams[key] = value;
      });

      const outWidget = createFuncTableViewWidget(sf, inputParams);
      host.appendChild(outWidget);
    }
  }
}

export function createFuncTableViewWidget(sf: DG.Func, inputParams: Record<string, any>): HTMLElement {
  const outWidget = DG.Widget.fromRoot(ui.wait(async () => {
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
        sf.inputs.map((i: any) => i.propertyType === DG.TYPE.STRING ? `"${inputParams[i.name]}"` : inputParams[i.name]).join(',');
      if (inputParams && Object.keys(inputParams).length > 0)
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
  }));
  outWidget.addProperty('caption', DG.TYPE.STRING, sf.friendlyName ?? sf.name);
  outWidget.props.caption = sf.friendlyName ?? sf.name;
  outWidget.root.style.minHeight = '500px';
  return widgetHost(outWidget);
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

let _currentSearchAppView: DG.View | DG.ViewBase | null = null;
function exactAppViewSearch(s: string) {
  s = s.toLowerCase().trim();
  const appFunc = exactAppFuncSearch(s);
  if (appFunc) {
    appFunc.apply({}).then((v) => {
      if (v && (v instanceof DG.View || v instanceof DG.ViewBase)) {
        try {
          if (_currentSearchAppView && _currentSearchAppView.root && document.body.contains(_currentSearchAppView.root))
            grok.shell.dockManager.close(_currentSearchAppView.root);
        } catch (e) {
          console.error(`Error closing previous app view: ${e}`);
        }
        _currentSearchAppView = v;
        grok.shell.dockManager.dock(v.root, DG.DOCK_TYPE.DOWN, grok.shell.dockManager.findNode(grok.shell.v.root), v.name ?? capitalizeFirstLetter(v.type));
      }
    }).catch((e) => console.error(`Error opening app view for ${s}: ${e}`));
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

function functionEvaluationsSearch(s: string, host: HTMLDivElement): void {
  if (!s.includes('(') || !s.includes(')'))
    return;

  grok.functions
    .eval(s)
    .then((result) => host.appendChild(ui.span([s + ' = ', result], {style: {fontSize: '20px'}})))
    .catch(() => {});
}
