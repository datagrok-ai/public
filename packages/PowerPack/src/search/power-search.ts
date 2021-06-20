/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from "../widgets/web-widget";
import {DataQuery} from "datagrok-api/dg";
import {widgetHost} from "../utils";

// Power Search: community-curated, template-based, widget-driven search engine

interface Template {
  template: string;
  url: string;
  regexp?: RegExp;   // cached regexp for the template
}

interface Card {
  id: string;
  name: string;
  widget: string;
  templates: Template[];
}

let queries: DG.DataQuery[] = [];
let widgetFunctions = DG.Func.find({returnType: 'widget'});
let searchFunctions = DG.Func.find({tags: ['search'], returnType: 'list'});
let searchWidgetFunctions = DG.Func.find({tags: ['search'], returnType: 'widget'});

export function initSearch() {
  grok.dapi.queries.list().then((qs) => queries = qs);
}

export function powerSearch(s: string, host: HTMLDivElement): void {
  ui.empty(host);

  jsEvalSearch(s, host) ||
  viewsSearch(s, host);

  projectsSearch(s, host);
  functionEvaluationsSearch(s, host);
  searchFunctionsSearch(s, host);
  queriesEntitiesSearch(s, host);
  queriesSearch(s, host);
  templatesSearch(s, host);
  widgetsSearch(s, host);
  specificWidgetsSearch(s, host);
}


/// Searches for default views that exactly match the search string
/// Example: functions
function viewsSearch(s: string, host: HTMLDivElement): boolean {
  if (DG.View.ALL_VIEW_TYPES.includes(s)) {
    host.appendChild(DG.View.createByType(s).root);
    return true;
  }
  return false;
}

function queriesEntitiesSearch(s: string, host: HTMLDivElement): void {
  s = s.toLowerCase();
  host.appendChild(ui.list(queries.filter((q) => q.name.toLowerCase().includes(s))));
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

/// Community-curated template collection
function templatesSearch(s: string, host: HTMLDivElement): void {
  for (let p of templates)
    for (let t of p.templates) {
      let x = <any>t;
      if (x.regexp == null)
        x.regexp = new RegExp(t.template);
      let matches = x.regexp.exec(s);

      if (matches !== null) {
        let widgetProperties: any = {};
        for (let [k, v] of Object.entries(t))
          if (k != 'template' && k != 'regexp') {
            for (let i = 1; i < matches.length; i++)
              v = v.replace('${' + i + '}', matches[i]);

            widgetProperties[k] = v;
          }

        DG.Func.byName(p.widget).apply().then((w: DG.Widget) => {
          w.props.setAll(widgetProperties);
          host.appendChild(w.root);
        });
      }
    }
}

//
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


const semTypes = [
  {
    name: 'CHEMBL_ID',
    description: 'ChEMBL compound identifier',
    template: '(CHEMBL[0-9]+)'
  }
];

const templates = [
  {
    id: 'kpi-test',
    name: 'KPI Test',
    widget: 'kpiWidget',
    templates: [
      {
        template: 'kpi ([0-9]+)',
        caption: '${1}'
      }
    ]
  },
  {
    id: 'ticket-stock-price',
    name: 'Stock Price',
    widget: 'webWidget',
    templates: [
      {
        template: 'kpi ([0-9]+)',
        caption: '${1}'
      }
    ]
  },
  {
    id: 'chembl-id-report-card',
    name: 'Report card',
    widget: 'webWidget',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/name_and_classification'
    }]
  },
  {
    id: 'chembl-id-representations',
    name: 'Representations',
    widget: 'webWidget',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/representations'
    }]
  },
  {
    id: 'chembl-id-alternative-forms',
    name: 'Alternative forms',
    widget: 'webWidget',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/alternate_forms'
    }]
  },
  {
    id: 'chembl-id-calculated-properties',
    name: 'Calculated properties',
    widget: 'webWidget',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/calculated_properties'
    }]
  },
  {
    id: 'chembl-id-cross-refs',
    name: 'Cross references',
    widget: 'webWidget',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/cross_refs'
    }]
  },
]

const widgetTemplates = [

];

// <object data="https://www.ebi.ac.uk/chembl/embed/#compound_report_card/CHEMBL1193654/name_and_classification" width="100%" height="100%"></object>