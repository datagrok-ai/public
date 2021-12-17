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
import { AEBrowserHelper } from './helpers/ae-browser-helper';
import { AE_END_DATE, AE_END_DAY, AE_SEVERITY, AE_START_DAY, AE_TERM, SUBJECT_ID } from './columns-constants';
import { STUDY_ID } from './columns-constants';
import { checkMissingColumns, checkMissingDomains } from './views/utils';
import { TreeMapView } from './views/tree-map-view';
import { MedicalHistoryView } from './views/medical-history-view';

export let _package = new DG.Package();

export let validationRulesList = null;

let domains = Object.keys(study.domains).map(it => `${it.toLocaleLowerCase()}.csv`);

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

  if (Object.keys(meta.domains).every((name) => grok.shell.table(name) == null)) {
    let demoFiles = await grok.dapi.projects.filter('clin-demo-files-2').list();
    if (demoFiles.length) {
      await grok.dapi.projects.open('clin-demo-files-2');
    } else {
      grok.shell.warning('Please load SDTM data or demo files');
    }
  }

  study.initFromWorkspace();

  function addView(view: DG.ViewBase): DG.ViewBase {
    view.box = true;
    view.parentCall = c;
    view.path = '/' + view.name;
    grok.shell.addView(view);
    return view;
  }

  const views = [];

  views.push(<StudySummaryView>addView(new StudySummaryView('Summary')));
  const timelinesView = new TimelinesView('Timelines');
  views.push(<TimelinesView>addView(timelinesView));
  views.push(<PatientProfileView>addView(new PatientProfileView('Patient Profile')));
  views.push(<AdverseEventsView>addView(new AdverseEventsView('Adverse Events')));
  views.push(<LaboratoryView>addView(new LaboratoryView('Laboratory')));
  //views.push(<AERiskAssessmentView>addView(new AERiskAssessmentView('AE Risk Assessment')));
  views.push(<SurvivalAnalysisView>addView(new SurvivalAnalysisView('Survival Analysis')));
  views.push(<BoxPlotsView>addView(new BoxPlotsView('Distributions')));
  views.push(<MatrixesView>addView(new MatrixesView('Correlations')));
  views.push(<TimeProfileView>addView(new TimeProfileView('Time Profile')));
  //views.push(<TreeMapView>addView(new TreeMapView('Tree map')));
  views.push(<MedicalHistoryView>addView(new MedicalHistoryView('Medical History')));

  let aeBrowserView;
  if (study.domains.ae) {
    aeBrowserView = DG.View.create();
    if (checkMissingColumns(aeBrowserView, ['ae'], { 'ae': {'req': [AE_TERM, AE_SEVERITY, AE_START_DAY, AE_END_DAY]} })) {
      const aeBrowserDf = study.domains.ae.clone();
      aeBrowserView = DG.TableView.create(aeBrowserDf);
      const aeBrowserHelper = new AEBrowserHelper(aeBrowserDf);
      timelinesView.aeBrowserHelper = aeBrowserHelper;
      aeBrowserDf.onCurrentRowChanged.subscribe(() => {
        aeBrowserHelper.currentSubjId = aeBrowserDf.get(SUBJECT_ID, aeBrowserDf.currentRowIdx);
        aeBrowserHelper.currentAeDay = aeBrowserDf.get(AE_START_DAY, aeBrowserDf.currentRowIdx);
        aeBrowserHelper.createAEBrowserPanel();
      })
    }
  } else {
    aeBrowserView = DG.View.create();
    checkMissingDomains({
      'req_domains': {
        'ae': { 
          'req': [AE_TERM, AE_SEVERITY, AE_START_DAY, AE_END_DAY] 
        } 
      }
    }, aeBrowserView);
  }
  aeBrowserView.name = 'AE browser';
  aeBrowserView.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/ae_browser.md';
  views.push(addView(aeBrowserView));

  DG.ObjectHandler.register(new AdverseEventHandler());

  let summary = views.find(it => it.name === 'Summary');
  summary.load();
  summary.loaded = true;
  let valView = addView(new ValidationView(summary.errorsByDomain, 'Validation'));
  summary.validationView = valView;
  views.push(valView);

  setTimeout(() => {
    grok.shell.v = summary;
  }, 1000);

  grok.events.onCurrentViewChanged.subscribe((v) => {
    setTimeout(() => {
      const obj = views.find(it => it.name === grok.shell.v.name);
      if (!obj.loaded) {
        obj.load();
        obj.loaded = true;
      }
      createPropertyPanel(obj);
    }, 100)
  });
}


//tags: folderViewer
//input: file folder
//input: list<file> files
//output: widget
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]): Promise<DG.Widget | undefined> {
  if (files.some((f) => f.fileName.toLowerCase() === 'dm.csv')) {
    let res = await grok.dapi.files.readAsText(`${folder.fullPath}/dm.csv`);
    let table = DG.DataFrame.fromCsv(res);
    let studyId = table.columns.names().includes(STUDY_ID) ? table.get(STUDY_ID, 0) : 'undefined';
    return DG.Widget.fromRoot(ui.div([
      ui.panel([
        ui.divText('Folder contains SDTM data'),
        ui.divText(`Study ID: ${studyId}`)]),
      ui.button('Run ClinicalCase', async () => {
        await Promise.all(files.map(async (file) => {
          if(domains.includes(file.fileName.toLowerCase())){
            let df = await grok.data.files.openTable(`${folder.fullPath}/${file.fileName.toLowerCase()}`);
            grok.shell.addTableView(df);
          }
        }));
        grok.functions.call("Clinicalcase:clinicalCaseApp");
      })
    ]));
  }
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

