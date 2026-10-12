/** The provenance dialog (F-28): from one portal object, what it rests on and what rests on it, as
 * two indented trees. Each node opens what it is (a run's call, a snapshot, a document's file, the
 * program data, a model card, a table); Trace starts the walk again from that node. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {portalDb} from '../generated/db';
import {provenanceOf} from '../provenance';
import {
  describeNode, PROVENANCE_STARTS, Provenance, ProvenanceRecords, provenanceTree, ReachedNode,
} from '../model/provenance';
import {openSnapshotDialog} from './snapshot-dialog';
import {openDocumentFile} from './document-dialog';
import {openProgramTable} from '../program-data';
import {showModelCard} from './model-card';
import {installStyles} from './styles';
import {named, openDialog, PortalDialog} from './u2-controls';

export interface ProvenanceDialogOptions {
  /** Opens a program on the program page; without it a program node only traces. */
  openProgram?: (programId: string) => void | Promise<void>;
  depth?: number;
}

/** Opens what [n] is. Returns false when the node has nothing to open. */
async function openNode(n: ReachedNode, records: ProvenanceRecords, opts: ProvenanceDialogOptions): Promise<boolean> {
  switch (n.kind) {
  case 'model_run': {
    const runId = records.runs.find((r) => r.id === n.id)?.runId;
    const call = runId ? await grok.dapi.functions.calls.find(runId).catch(() => null) : null;
    if (!call)
      return false;
    grok.shell.o = call;
    return true;
  }
  case 'snapshot':
    await openSnapshotDialog(n.id);
    return true;
  case 'document': {
    const d = records.documents.find((x) => x.id === n.id);
    if (!d)
      return false;
    await openDocumentFile({title: d.title, path: d.path, doc_version: d.version, status: d.status});
    return true;
  }
  case 'dataset': {
    const program = await portalDb.programs.get(n.id.split(':')[0]);
    const portfolio = await portalDb.portfolios.get(program.portfolio_id);
    await openProgramTable(portfolio.code, program.code);
    return true;
  }
  case 'model':
    await showModelCard(n.id);
    return true;
  case 'table':
    grok.shell.addTableView(await grok.dapi.tables.getTable(n.id));
    return true;
  case 'url':
    window.open(n.id, '_blank', 'noopener');
    return true;
  case 'file':
    if (/\.(csv|tsv)$/i.test(n.id))
      grok.shell.addTableView(await grok.dapi.files.readCsv(n.id));
    else
      await openDocumentFile({title: n.name, path: n.id});
    return true;
  case 'program':
    if (!opts.openProgram)
      return false;
    await opts.openProgram(n.id);
    return true;
  default:
    return false;
  }
}

/** The start kind a node can be traced from, if any. */
const traceable = (n: ReachedNode): string | null => (PROVENANCE_STARTS as string[]).includes(n.kind) ? n.kind : null;

function tree(p: Provenance & {records: ProvenanceRecords}, direction: 'up' | 'down', opts: ProvenanceDialogOptions,
  retrace: (kind: string, id: string) => void, close: () => void): HTMLElement {
  const lines = provenanceTree(p, direction);
  if (!lines.length)
    return named(span('Nothing', 'pp-muted'), `provenance-${direction}-none`);
  return named(divV(lines.map((l, i) => {
    const name = `provenance-${direction}-${i}`;
    const label = named(span(describeNode(l.node), 'pp-link'), `${name}-open`);
    label.title = `${l.node.relation}: open it`;
    label.onclick = () => void openNode(l.node, p.records, opts).then((opened) => {
      if (!opened && traceable(l.node))
        retrace(l.node.kind, l.node.id);
      else if (opened && l.node.kind === 'program')
        close();
    }).catch((e) => grok.shell.error(e?.message ?? String(e)));
    const items: HTMLElement[] = [label, span(l.node.relation ?? '', 'pp-muted')];
    if (traceable(l.node)) {
      const trace = named(span('Trace', 'pp-link'), `${name}-trace`);
      trace.title = 'Trace provenance from here';
      trace.onclick = () => retrace(l.node.kind, l.node.id);
      items.push(trace);
    }
    const row = named(divH(items, 'pp-unit'), name);
    row.style.paddingLeft = `${l.level * 14}px`;
    row.dataset.kind = l.node.kind;
    return row;
  }), 'pp-provenance-tree'), `provenance-${direction}`);
}

/** Opens the provenance of [kind] [id] in a dialog. */
export async function openProvenanceDialog(kind: string, id: string,
  opts: ProvenanceDialogOptions = {}): Promise<PortalDialog> {
  installStyles();
  const body = named(divV([span('Reading the portal\'s records…', 'pp-muted')], 'pp-provenance'), 'provenance-content');
  const dlg = openDialog({title: 'Provenance', name: 'provenance-dialog', width: 620}, () => ({content: body}));
  const draw = async (k: string, i: string) => {
    try {
      const p = await provenanceOf(k, i, opts.depth ?? 4);
      const start = p.nodes.find((n) => n.key === p.start)!;
      const retrace = (kk: string, ii: string) => void draw(kk, ii);
      body.replaceChildren(
        named(span(describeNode(start), 'pp-title'), 'provenance-start'),
        span(`${p.nodes.length - 1} record${p.nodes.length === 2 ? '' : 's'} within ${opts.depth ?? 4} steps; click one to open it`, 'pp-muted'),
        span('Rests on', 'pp-heading'),
        tree(p, 'up', opts, retrace, () => dlg.close()),
        span('Used by', 'pp-heading'),
        tree(p, 'down', opts, retrace, () => dlg.close()),
      );
    } catch (e: any) {
      body.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'provenance-error'));
    }
  };
  await draw(kind, id);
  return dlg;
}
