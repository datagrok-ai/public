/** A program's data package manifest (F-43), as a view: what the program's gate case rests on, one
 * line each (criteria, documents, evidence tables, runs, datasets, snapshots), with status and where
 * each comes from. **Open as table** puts it in a table view; **Download CSV** saves it. Opened from
 * the program page's Data package button. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {dataPackageFrame, dataPackageOf} from '../data-package';
import {MANIFEST_COLUMNS, manifestCell, manifestCsv, manifestFileName, summarizeManifest} from '../model/data-package';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

export async function openDataPackageDialog(programId: string): Promise<DG.Dialog> {
  installStyles();
  const p = await dataPackageOf(programId);
  const t = named(document.createElement('table'), 'data-package-lines');
  t.className = 'pp-card-table';
  const h = t.insertRow();
  for (const c of MANIFEST_COLUMNS.filter((x) => x.key !== 'ref'))
    h.append(Object.assign(document.createElement('th'), {textContent: c.name}));
  p.lines.forEach((l, i) => {
    const tr = named(t.insertRow(), `data-package-line-${i}`);
    tr.dataset.section = l.section;
    tr.title = l.ref;
    for (const c of MANIFEST_COLUMNS.filter((x) => x.key !== 'ref'))
      tr.insertCell().textContent = manifestCell(l, c.key);
  });
  const csv = named(ui.button('Download CSV', () =>
    DG.Utils.download(manifestFileName(p.program.code, new Date()), manifestCsv(p.lines), 'text/csv')), 'data-package-csv');
  csv.title = 'The manifest as a CSV file, one line per item';
  const asTable = named(ui.button('Open as table', () => {
    grok.shell.addTableView(dataPackageFrame(p));
    dlg.close();
  }), 'data-package-open-table');
  asTable.title = 'The manifest in a table view: sort, filter, export';
  const dlg = ui.dialog({title: `Data package · ${p.program.code}`});
  dlg.root.dataset.u2Name = 'data-package-dialog';
  const scroll = divV([p.lines.length ? t : span('Nothing recorded yet', 'pp-muted')]);
  scroll.style.cssText = 'max-height: 460px; overflow: auto;';
  dlg.add(named(divV([
    named(span(`${p.program.code} · ${p.program.name}: ${summarizeManifest(p.lines)}`), 'data-package-summary'),
    span('Read from the portal\'s records now; hover a line for what it refers to.', 'pp-muted'),
    scroll,
    named(divH([asTable, csv], 'pp-actions'), 'data-package-actions'),
  ], 'pp-data-package'), 'data-package-content'));
  dlg.show({width: 900});
  return dlg;
}
