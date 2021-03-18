/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as meta from './sdtm-meta';

export let _package = new DG.Package();

let links = {
  ae: { key: 'USUBJID', start: 'AESTDY', end: 'AEENDY', event: 'AETERM'},
  cm: { key: 'USUBJID', start: 'VISITDY', event: 'CMTRT'},
  ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT'},
  lb: { key: 'USUBJID', start: 'LBDY', event: 'LBTEST'}
};

let typeMap = { 'Char': 'string', 'Num': 'int' };

let terminology, submissionValueCol;
let submissionValues = [];

//name: SDTM Summary
//tags: panel, widgets
//input: dataframe df
//output: widget result
//condition: df.tags.get("sdtm")
export function sdtmSummaryPanel(df) {
  let domain = df.getTag('sdtm-domain');
  let text = `SDTM domain: ${domain.toUpperCase()}\n`;
  for (let column of df.columns) {
    let name = column.name;
    let variable = meta.domains[domain][name];
    text += `${name} ${variable ?
      typeMap[variable.type] === column.type ?
      'valid' : 'invalid' : 'unknown variable'}\n`;
  }
  return new DG.Widget(ui.divText(text));
}

//name: SDTM Variable
//tags: panel, widgets
//input: column varCol
//output: widget result
//condition: t.tags.get("sdtm")
export function sdtmVariablePanel(varCol) {
  let domain = meta.domains[varCol.dataFrame.getTag('sdtm-domain')];
  let variable = domain[varCol.name];
  let text = `${varCol.name}\n${variable ?
    variable.label + '\nType: ' + (typeMap[variable.type] === varCol.type ?
    'valid' : 'invalid') : 'Unknown variable'}\n`;
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
        [...synCol.get(i).split('; '), nciTermCol.get(i)].forEach(s => {
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
  if (outliers) container.push(outliers, convertButton);
  return new DG.Widget(ui.divV(container));
}


//name: Clinical Case
//tags: app
export function clinicalCaseApp() {
  grok.shell.info('This is clinical.');
}

//tags: autostart
export async function clinicalCaseInit() {
  terminology = await grok.data.loadTable(`${_package.webRoot}tables/sdtm-terminology.csv`);
  submissionValueCol = terminology.getCol('CDISC Submission Value');
  submissionValues = submissionValueCol.categories;

  grok.events.onTableAdded.subscribe(args => {
    let t = args.args.dataFrame;
    let domain = meta.domains[t.name.toLowerCase()];
    if (domain) {
      t.setTag('sdtm', true);
      t.setTag('sdtm-domain', t.name.toLowerCase());
      for (let variableName in domain)
        if (t.columns.contains(variableName)) {
          //t.col(variableName).semType = 'sdtm-' + t.name.toLowerCase() + '-' + variableName;
          t.col(variableName).setTag(DG.TAGS.DESCRIPTION, domain[variableName]['label']);
        }
    }

    if (Object.keys(links).every(key => grok.shell.tableByName(key))) {
      grok.shell.topMenu.group('Clin').item('Timelines', () => clinicalCaseTimelines());
    }
  });
}

//name: clinicalCaseTimelines
export function clinicalCaseTimelines() {

  let result = null;

  let getTable = function(domain) {
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
