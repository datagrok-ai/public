/** A program's data package manifest (F-43), as a view: what the program's gate case rests on, one
 * line each (criteria, documents, evidence tables, runs, datasets, snapshots), with status and where
 * each comes from. **Open as table** puts it in a table view; **Download CSV** saves it. Opened from
 * the program page's Data package button. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {button, divH, divV, span} from '@datagrok-libraries/u2';
import {dataPackageFrame, dataPackageOf} from '../data-package';
import {MANIFEST_COLUMNS, manifestCell, manifestCsv, manifestFileName, summarizeManifest} from '../model/data-package';
import {installStyles} from './styles';
import {named, openDialog, PortalDialog, table} from './u2-controls';

export async function openDataPackageDialog(programId: string): Promise<PortalDialog> {
  installStyles();
  const p = await dataPackageOf(programId);
  const t = table('data-package-lines', MANIFEST_COLUMNS.filter((x) => x.key !== 'ref').map((c) =>
    ({header: c.name, render: (l: typeof p.lines[number]) => manifestCell(l, c.key)})), p.lines,
  {name: (_, i) => `data-package-line-${i}`, title: (l) => l.ref}).root;
  const csv = named(button('Download CSV', () =>
    DG.Utils.download(manifestFileName(p.program.code, new Date()), manifestCsv(p.lines), 'text/csv')), 'data-package-csv');
  csv.title = 'The manifest as a CSV file, one line per item';
  const asTable = named(button('Open as table', () => {
    grok.shell.addTableView(dataPackageFrame(p));
    dlg.close();
  }), 'data-package-open-table');
  asTable.title = 'The manifest in a table view: sort, filter, export';
  const scroll = divV([p.lines.length ? t : span('Nothing recorded yet', 'pp-muted')]);
  scroll.style.cssText = 'max-height: 460px; overflow: auto;';
  const content = named(divV([
    named(span(`${p.program.code} · ${p.program.name}: ${summarizeManifest(p.lines)}`), 'data-package-summary'),
    span('Read from the portal\'s records now; hover a line for what it refers to.', 'pp-muted'),
    scroll,
    named(divH([asTable, csv], 'pp-actions'), 'data-package-actions'),
  ], 'pp-data-package'), 'data-package-content');
  const dlg = openDialog({title: `Data package · ${p.program.code}`, name: 'data-package-dialog', width: 900},
    () => ({content}));
  return dlg;
}
