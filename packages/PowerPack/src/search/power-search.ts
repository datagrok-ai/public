/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from "../widgets/web-widget";
import {DataQuery} from "datagrok-api/dg";
import {widgetHost} from "../utils";
import {_package} from "../package";
import {tryParseJson} from "@datagrok-libraries/utils/src/string-utils";
import {initTemplates, templatesSearch} from "./templates-search";

// Power Search: community-curated, template-based, widget-driven search engine

let queries: DG.DataQuery[] = [];
let widgetFunctions = DG.Func.find({returnType: 'widget'});
let searchFunctions = DG.Func.find({tags: ['search'], returnType: 'list'});
let searchWidgetFunctions = DG.Func.find({tags: ['search'], returnType: 'widget'});

export function initSearch() {
  //grok.dapi.queries.list().then((qs) => queries = qs);
  initTemplates();
}

export function powerSearch(s: string, host: HTMLDivElement): void {
  ui.empty(host);

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

function queriesEntitiesSearch(s: string, host: HTMLDivElement): void {
  s = s.toLowerCase();
  host.appendChild(ui.list(queries.filter((q) => q.name.toLowerCase().includes(s))));
}


function regexEntitiesSearch(s: string, host: HTMLDivElement): void {
  var ids = s.split(/,\s*/);
  var semValues = ids.map(id => DG.SemanticValue.parse(id));
  if (semValues.every(sv => sv?.semType && sv.semType == semValues[0].semType) && DG.ObjectHandler.forEntity(semValues[0])) {
    for (let sv of semValues)
      host.append(DG.ObjectHandler.forEntity(sv)!.renderCard(sv));
  }
}


function projectsSearch(s: string, host: HTMLDivElement): void {
  if (s.length < 3)
    return;
  grok.dapi.projects.filter(s).list({pageSize: 5}).then((projects) => {
    if (projects.length > 0)
      host.appendChild(
        ui.divV([
          ui.h3('Projects'),
          ui.divH(projects.map((p) => ui.renderCard(p)))
        ])
      );
  });
}

/// Evaluates some of the JavaScript code
/// Example: ui.label('produced by JavaScript')
function jsEvalSearch(s: string, host: HTMLDivElement): boolean {
  if (s.startsWith('grok.') || s.startsWith('DG.') || s.startsWith('ui.')) {
    try {
      let x = eval(s);
      host.appendChild(ui.render(x));
      return true;
    }
    catch (_) {}
  }
  return false;
}

/// Evaluates custom search functions
function searchFunctionsSearch(s: string, host: HTMLDivElement): void {
  for (let sf of searchFunctions)
    sf.apply({s: s}).then((results: any[]) => {
      if (results.length > 0) {
        host.appendChild(ui.divV([
          ui.h3(sf.description ?? sf.name),
          ui.list(results)
        ]));
      }
    });
}

/// Special widgets
function widgetsSearch(s: string, host: HTMLDivElement): void {
  for (let sf of searchWidgetFunctions)
    sf.apply({s: s})
      .then((result: DG.Widget) => {
        if (result)
          host.appendChild(widgetHost(result));
      });
}

/// Explicitly spelled widgets
/// Example: "kpiWidget"
function specificWidgetsSearch(s: string, host: HTMLDivElement): void {
  s = s.toLowerCase();
  for (let wf of widgetFunctions)
    if (wf.name.toLowerCase() == s)
      wf.apply().then((w: DG.Widget) => host.appendChild(ui.div([widgetHost(w)])));
}


export function queriesSearch(s: string, host: HTMLDivElement): void {
  const dateRegExp = new RegExp('\\btoday\\b|\\bthis week\\b');
  const varRegExpStr = '@([a-zA-Z0-9_]+)';
  const varRegExp = new RegExp(varRegExpStr);
  const idRegExp = new RegExp('([a-zA-Z0-9_]+)');
  const byRegExp = new RegExp('by ([a-zA-Z0-9_]+)$');

  let byColumn: string | null = null;
  if (byRegExp.test(s)) {
    let byMatches = byRegExp.exec(s);
    if (byMatches != null && byMatches.length == 2) {
      s = s.replace(byRegExp, '');
      byColumn = byMatches[1];
    }
  }

  for (let q of queries.filter((q) => varRegExp.exec(q.friendlyName) !== null)) {
    let varMatches = varRegExp.exec(q.friendlyName);
    if (varMatches == null || varMatches.length != 2)
      continue;

    let qPattern = new RegExp(q.friendlyName.replace(varRegExp, '([a-zA-Z0-9_]+)'));
    let matches = qPattern.exec(s);
    if (matches != null && matches.length == 2) {
      let varName = varMatches[1];
      let varValue = matches[1];
      q.apply({ [varName] : varValue }).then((df: DG.DataFrame) => {
        let viewer: DG.Viewer = (byColumn == null)
          ? DG.Viewer.grid(df)
          : DG.Viewer.barChart(df, { split: byColumn });

        host.appendChild(ui.render(q));
        host.appendChild(viewer.root);
      });
    }
  }
}

function functionEvaluationsSearch(s: string, host: HTMLDivElement): void {
  if (!s.includes('(') && s.includes(')'))
    return;

  grok.functions
    .eval(s)
    .then((result) => host.appendChild(ui.span([s + ' = ', result], {style: {'font-size' : '20px'}})))
    .catch(() => {})
}
