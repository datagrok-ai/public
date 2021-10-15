/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as meta from './sdtm-meta';
import {study} from "./clinical-study";
import {StudySummaryView} from "./views/study-summary-view";
import {TimelinesView} from "./views/timelines-view";
import {PatientProfileView} from "./views/patient-profile-view";
import {AdverseEventsView} from "./views/adverse-events-view";
import {ValidationView} from './views/validation-view';
import {AdverseEventHandler} from './panels/adverse-event-handler';
import {LaboratoryView} from './views/laboratory-view';
import {AERiskAssessmentView} from './views/ae-risk-assessment-view';
import {SurvivalAnalysisView} from './views/survival-analysis-view';
import { BoxPlotsView } from './views/boxplots-view';
import { MatrixesView } from './views/matrixes-view';
import { createPropertyPanel } from './panels/panels-service';
import { TimeProfileView } from './views/time-profile-view';
import { AeBrowserView } from './views/adverse-events-browser';

export let _package = new DG.Package();

export let validationRulesList = null;


let links = {
  ae: { key: 'USUBJID', start: 'AESTDY', end: 'AEENDY', event: 'AETERM' },
  cm: { key: 'USUBJID', start: 'VISITDY', event: 'CMTRT' },
  ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT' },
  lb: { key: 'USUBJID', start: 'LBDY', event: 'LBTEST' }
};

let typeMap = { 'Char': 'string', 'Num': 'int' };
let datetimeFormat = 'ISO 8601';
let checkType = (column: DG.Column, variable) => (column.type === typeMap[variable.type] ||
  (column.type === 'datetime' && variable.format === datetimeFormat));

let applyRule = (column: DG.Column, variable) => {
  // let matcher = (variable.type === 'Char') ? DG.ValueMatcher.string(variable.rule) : DG.ValueMatcher.numerical(variable.rule);
  if (variable.type === 'Char') {
    let matcher = DG.ValueMatcher.string(variable.rule);

    for (let i = 0; i < column.length; i++) {
      if (!matcher.match(column.get(i))) return false;
    }
    return true;
  }
  return false;
};

let terminology: DG.DataFrame;
let submissionValueCol: DG.Column;
let submissionValues = [];

//name: SDTM Summary
//tags: panel
//output: widget result
//condition: true
export function sdtmSummaryPanel(): DG.Widget {
  return new DG.Widget(ui.divText('test'));
}

//name: SDTM Variable
//tags: panel, widgets
//input: column varCol
//output: view result
//condition: true
export function sdtmVariablePanel(varCol: DG.Column): DG.Widget {
  let domain = meta.domains[varCol.dataFrame.getTag('sdtm-domain')];
  let variable = domain[varCol.name];
  let text = `${varCol.name}\n${variable ?
    variable.label + '\nType: ' + (checkType(varCol, variable) ?
      'valid' : 'invalid') : 'Unknown variable'}\n`;
  let missingValueCount = varCol.stats.missingValueCount;
  let convertButton, outliers;

  let isTerm = submissionValues.includes(varCol.name);
  text += `CDISC Submission Value: ${isTerm}\n`;

  if (isTerm) {
    let match = terminology.rows.match({ 'CDISC Submission Value': varCol.name }).toDataFrame();
    text += `CDISC Synonym(s): ${match.get('CDISC Synonym(s)', 0)}\n`;
    text += `CDISC Definition: ${match.get('CDISC Definition', 0)}\n`;
    text += `NCI Preferred Term: ${match.get('NCI Preferred Term', 0)}\n`;

    let relatedRecords = terminology.rows.match({ 'Codelist Code': match.get('Code', 0) }).toDataFrame();
    let rowCount = relatedRecords.rowCount;
    if (rowCount) {
      let valueCol = relatedRecords.getCol('CDISC Submission Value');
      let synCol = relatedRecords.getCol('CDISC Synonym(s)');
      let nciTermCol = relatedRecords.getCol('NCI Preferred Term');
      let synonyms = {};

      for (let i = 0; i < rowCount; i++) {
        let submissionValue = valueCol.get(i);
        [...synCol.get(i).split('; '), nciTermCol.get(i), submissionValue].forEach(s => {
          if (s) synonyms[s.toLowerCase()] = submissionValue;
        });
      }

      let valuesToConvert = varCol.categories.filter(v => v && !valueCol.categories.includes(v));
      if (valuesToConvert.length) {
        outliers = ui.divText('Out-of-vocabulary values:\n' + valuesToConvert.join(', '));
        outliers.style = 'color: red';
        convertButton = ui.button('Convert', () => {
          varCol.init(i => synonyms[varCol.get(i).toLowerCase()] || varCol.get(i));
        }, 'Convert to CDISC submission values');
      }
    }
  }
  let container = [ui.divText(text)];
  if (missingValueCount) container.push(ui.divText(`Missing values: ${missingValueCount}`));
  if (outliers) container.push(outliers, convertButton);
  return new DG.Widget(ui.divV(container));
}

//name: Clinical Case
//tags: app
export async function clinicalCaseApp(): Promise<any> {
  let c: DG.FuncCall = grok.functions.getCurrentCall();
  validationRulesList = await grok.data.loadTable(`${_package.webRoot}tables/validation-rules.csv`);

  if (Object.keys(meta.domains).every((name) => grok.shell.table(name) == null))
    await grok.dapi.projects.open('clin-demo-files-2');

  study.initFromWorkspace();

  function addView(view: DG.ViewBase): DG.ViewBase {
    view.box = true;
    view.parentCall = c;
    view.path = '/'+view.name;
    grok.shell.addView(view);
    return view;
  }

  const views = [];

  views.push(<StudySummaryView>addView(new StudySummaryView('Summary')));
  views.push(<TimelinesView>addView(new TimelinesView('Timelines')));
  views.push(<PatientProfileView>addView(new PatientProfileView('Patient Profile')));
  views.push(<AdverseEventsView>addView(new AdverseEventsView('Adverse Events')));
  views.push(<LaboratoryView>addView(new LaboratoryView('Laboratory')));
  views.push(<AERiskAssessmentView>addView(new AERiskAssessmentView('AE Risk Assessment')));
  views.push(<SurvivalAnalysisView>addView(new SurvivalAnalysisView('Survival Analysis')));
  views.push(<BoxPlotsView>addView(new BoxPlotsView('Biomarkers distribution')));
  views.push(<MatrixesView>addView(new MatrixesView('Correlation Matrix')));
  views.push(<TimeProfileView>addView(new TimeProfileView('Time Profile')));
  views.push(<AeBrowserView>addView(new AeBrowserView('AE Browser')));
  DG.ObjectHandler.register(new AdverseEventHandler());

  let summary = views.find(it => it.name === 'Summary');
  summary.validationView = addView(new ValidationView(summary.errorsByDomain, 'Validation'));
  views.push(summary.validationView);

  grok.shell.v = summary;
  createPropertyPanel(summary);
  // showStudySummary();
  // showLabs();

  grok.events.onCurrentViewChanged.subscribe((layout) => {
    setTimeout(() => {
      const obj = views.find(it => it.name === grok.shell.v.name);
      createPropertyPanel(obj);
    }, 100)
  });
}

//tags: autostart
export async function clinicalCaseInit(): Promise<void> {
  terminology = await grok.data.loadTable(`${_package.webRoot}tables/sdtm-terminology.csv`);
  submissionValueCol = terminology.getCol('CDISC Submission Value');
  submissionValues = submissionValueCol.categories;

  grok.events.onTableAdded.subscribe(args => {
    let t = args.args.dataFrame;
    let domain = meta.domains[t.name.toLowerCase()];
    if (domain) {
      t.setTag('sdtm', 'true');
      t.setTag('sdtm-domain', t.name.toLowerCase());
      for (let variableName in domain)
        if (t.columns.contains(variableName)) {
          //t.col(variableName).semType = 'sdtm-' + t.name.toLowerCase() + '-' + variableName;
          t.col(variableName).setTag(DG.TAGS.DESCRIPTION, domain[variableName]['label']);
        }
    }

    if (Object.keys(links).every(key => grok.shell.tableByName(key))) {
      let clinMenu = grok.shell.topMenu.group('Clin');
      if (!clinMenu.find('Timelines'))
        clinMenu.item('Timelines', () => clinicalCaseTimelines());
      if (!clinMenu.find('Study Summary'))
        clinMenu.item('Study Summary', () => grok.shell.addView(new StudySummaryView('Summary')));
    }
  });
}

//name: clinicalCaseTimelines
export function clinicalCaseTimelines(): void {

  let result = null;

  let getTable = function (domain: string) {
    let info = links[domain];
    let t = grok.shell
      .tableByName(domain)
      .clone(null, Object.keys(info).map(e => info[e]));
    t.columns.addNew('domain', DG.TYPE.STRING).init(domain);
    for (let name in info)
      t.col(info[name]).name = name;
    return t;
  }

  for (let domain in links) {
    let t = getTable(domain);
    if (result == null)
      result = t;
    else
      result.append(t, true);
  }

  let v = grok.shell.addTableView(result);
  v.addViewer('TimelinesViewer');
}

function clearStdView(): void {
  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showProperties = false;
  windows.showConsole = false;
  windows.showHelp = false;
}

//name: showLabs
export function showLabs(): void {
  let lb = grok.shell.tableByName('lb');
  if (lb == null) return grok.shell.warning('Laboratory test results not found.');

  grok.shell.newView('Labs', [
    ui.divText('Labs view content'),
  ]);
}

