/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from "../widgets/web-widget";
import {DataQuery} from "datagrok-api/dg";
import {card} from "../utils";

// Power Search: community-curated, template-based, widget-driven search engine

interface Template {
  template: string;
  url: string;
  regexp?: RegExp;   // cached regexp for the template
}

interface Card {
  id: string;
  name: string;
  templates: Template[];
}

let queries: DG.DataQuery[] = [];
let searchFunctions = DG.Func.find({tags: ['search'], returnType: 'list'});
let searchWidgetFunctions = DG.Func.find({tags: ['search'], returnType: 'widget'});

export function initSearch() {
  grok.dapi.queries.list().then((qs) => queries = qs);
}

export function powerSearch(s: string, host: HTMLDivElement): void {
  viewsSearch(s, host);
  functionEvaluationsSearch(s, host);
  searchFunctionsSearch(s, host);
  queriesEntitiesSearch(s, host);
  queriesSearch(s, host);
  templatesSearch(s, host);
  widgetsSearch(s, host);
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
  host.appendChild(ui.list(queries.filter((q) => q.name.includes(s))));
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
          host.appendChild(card(result));
      });
}

/// Community-curated template collection
function templatesSearch(s: string, host: HTMLDivElement): void {
  for (let p of templates)
    for (let t of p.templates) {
      if (t.regexp == null)
        t.regexp = new RegExp(t.template);
      let matches = t.regexp.exec(s);
      let url = t.url;

      if (matches !== null) {
        for (let i = 1; i < matches.length; i++)
          url = url.replace('${' + i + '}', matches[i]);

        let widget = new WebWidget({
          src: url,
          width: '100%',
          height: '500px'
        });

        host.appendChild(widget.root);
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

  // let searchDateMatches = dateRegExp.exec(s);
  // if (searchDateMatches == null)
  //   return;
  //
  // function hasDateInput(q: DG.DataQuery): boolean {
  //   let matches = varRegExp.exec(q.friendlyName);
  //   if (matches == null)
  //     return false;
  //   return matches.length == 2 && matches[1] == '@date';
  // }
  //
  // // special case for date templates - do we even need it?
  // for (let q of queries.filter(hasDateInput)) {
  //   q.apply({date: searchDateMatches[1]}).then((df: DG.DataFrame) => {
  //     let grid = DG.Viewer.grid(df);
  //     host.appendChild(grid.root);
  //   });
  // }

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

const templates: Card[] = [
  {
    id: 'chembl-id-report-card',
    name: 'Report card',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/name_and_classification'
    }]
  },
  {
    id: 'chembl-id-representations',
    name: 'Representations',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/representations'
    }]
  },
  {
    id: 'chembl-id-alternative-forms',
    name: 'Alternative forms',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/alternate_forms'
    }]
  },
  {
    id: 'chembl-id-calculated-properties',
    name: 'Calculated properties',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/calculated_properties'
    }]
  },  {
    id: 'chembl-id-cross-refs',
    name: 'Cross references',
    templates: [{
      template: '(CHEMBL[0-9]+)',
      url: 'https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/cross_refs'
    }]
  },
]

const widgetTemplates = [

];

// <object data="https://www.ebi.ac.uk/chembl/embed/#compound_report_card/CHEMBL1193654/name_and_classification" width="100%" height="100%"></object>