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

//name: Clinical Case
//tags: app
export function clinicalCaseApp() {
  grok.shell.info('This is clinical.');
}

//tags: autostart
export function clinicalCaseInit() {
  grok.events.onTableAdded.subscribe((args) => {
    let t = args.args.dataFrame;
    let domain = meta.domains[t.name.toLowerCase()];
    if (domain) {
      for (let variableName in domain)
        if (t.columns.contains(variableName)) {
          t.col(variableName).semType = 'sdtm-' + t.name.toLowerCase() + '-' + variableName;
          t.col(variableName).setTag(DG.TAGS.DESCRIPTION, domain[variableName]);
        }
    }

    grok.shell.info('foo!');
    if (Object.keys(links).every((key) => grok.shell.tableByName(key))) {
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