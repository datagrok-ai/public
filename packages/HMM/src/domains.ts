/** Protein domain annotation: `hmmscan` of every sequence of a column against
 * a model library (Z = number of models), keeping reported domains (Pfam
 * gathering cutoffs or E-value thresholds), drawn as Bio region annotations
 * with a summary column. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {SearchOptions} from './hmmer/engine.ts';
import {HmmerPool} from './pool';
import {columnResidues} from './sequences';
import {type DomainHitSpan, setDomainAnnotations} from './annotations';
import {fetchPfam, loadBuiltIn, loadFile, type LoadedModels, parseAccessions} from './models';

export interface FoundDomain {
  row: number;
  model: string;
  accession: string;
  description: string;
  from: number;
  to: number;
  score: number;
  evalue: number;
}

/** Scan `column` against `models`; with gathering cutoffs, models without GA
 * cutoffs make HMMER stop, as `hmmscan --cut_ga` does. Coordinates are
 * 1-based residue indices of the gap-free sequence. */
export async function findDomains(pool: HmmerPool, models: LoadedModels, column: DG.Column<string>,
  options: SearchOptions): Promise<{found: FoundDomain[]; positions: number[][]}> {
  const {residues, positions} = await columnResidues(column);
  const rows: number[] = [];
  const sequences: {name: string; residues: string}[] = [];
  residues.forEach((r, i) => {
    if (!r) return;
    rows.push(i);
    sequences.push({name: String(i), residues: r});
  });
  const results = await pool.scan(models.key, sequences, options);
  const found: FoundDomain[] = [];
  results.forEach((result, q) => {
    if (result.status === 3) throw new Error('The scan stopped: a model lacks the requested cutoffs');
    for (const hit of result.hits) {
      const m = models.info.models[hit.target];
      for (const d of hit.domains) {
        if (!d.reported) continue;
        found.push({row: rows[q], model: m.name, accession: m.accession, description: m.description,
          from: d.envFrom, to: d.envTo, score: d.bitscore, evalue: d.iEvalue});
      }
    }
  });
  found.sort((a, b) => a.row - b.row || a.from - b.from);
  return {found, positions};
}

/** Annotate the column and add a `<column> domains` summary ("V-set 1-117; C1-set 140-223"). */
export function applyDomains(table: DG.DataFrame, column: DG.Column<string>, found: FoundDomain[],
  positions: number[][]): void {
  const perRow: DomainHitSpan[][] = Array.from({length: table.rowCount}, () => []);
  const summary = DG.Column.string(table.columns.getUnusedName(`${column.name} domains`), table.rowCount);
  const models = new Map<string, {id: string; name: string; description: string}>();
  const text: string[][] = Array.from({length: table.rowCount}, () => []);
  for (const d of found) {
    const id = d.accession || d.model;
    if (!models.has(id)) models.set(id, {id, name: d.model, description: d.description || d.model});
    perRow[d.row].push({model: id, from: d.from - 1, to: d.to - 1, score: d.score});
    text[d.row].push(`${d.model} ${d.from}-${d.to}`);
  }
  for (let i = 0; i < table.rowCount; i++) summary.set(i, text[i].join('; '), false);
  table.columns.insert(summary, table.columns.toList().indexOf(column) + 1);
  setDomainAnnotations(table, column, [...models.values()], perRow, positions);
}

/** `Bio | Annotate | Find Domains (HMMER)...` */
export async function showFindDomainsDialog(pool: HmmerPool, pkg: DG.Package): Promise<void> {
  const df = grok.shell.tv?.dataFrame;
  const seqCols = df?.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE) ?? [];
  if (!df || seqCols.length === 0) {
    grok.shell.warning('Open a table with a Macromolecule column');
    return;
  }
  const libraries = ['Pfam biologics library (built-in)', 'Pfam accessions (online)', 'HMM file'];
  const columnInput = ui.input.column('Sequence', {table: df, value: seqCols[0],
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MACROMOLECULE});
  const libraryInput = ui.input.choice('Library', {value: libraries[0], items: libraries});
  const accessionsInput = ui.input.string('Accessions', {value: 'PF07686, PF07654',
    tooltipText: 'Pfam accessions separated by commas'});
  const fileInput = ui.input.file('HMM file');
  const thresholds = ['Pfam gathering cutoffs (GA)', 'E-value'];
  const thresholdInput = ui.input.choice('Threshold', {value: thresholds[0], items: thresholds});
  const evalueInput = ui.input.float('E-value ≤',
    {value: 0.01, tooltipText: 'Sequence and domain E-value (-E, --domE)'});
  const refresh = () => {
    accessionsInput.root.style.display = libraryInput.value === libraries[1] ? '' : 'none';
    fileInput.root.style.display = libraryInput.value === libraries[2] ? '' : 'none';
    evalueInput.root.style.display = thresholdInput.value === thresholds[1] ? '' : 'none';
  };
  libraryInput.onChanged.subscribe(refresh);
  thresholdInput.onChanged.subscribe(refresh);
  refresh();

  ui.dialog({title: 'Find Domains (HMMER)'})
    .add(ui.inputs([columnInput, libraryInput, accessionsInput, fileInput, thresholdInput, evalueInput]))
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Finding domains...');
      try {
        let models: LoadedModels;
        if (libraryInput.value === libraries[0]) models = await loadBuiltIn(pool, pkg);
        else if (libraryInput.value === libraries[1])
          models = await fetchPfam(pool, parseAccessions(accessionsInput.value));
        else {
          if (!fileInput.value) throw new Error('Choose an HMM file');
          models = await loadFile(pool, fileInput.value);
        }
        const evalue = evalueInput.value ?? 0.01;
        const options: SearchOptions = thresholdInput.value === thresholds[0] ? {cutoff: 'ga'} :
          {evalue, domEvalue: evalue};
        const column = columnInput.value as DG.Column<string>;
        const {found, positions} = await findDomains(pool, models, column, options);
        applyDomains(df, column, found, positions);
        grok.shell.info(`${found.length} domains in ${new Set(found.map((d) => d.row)).size} sequences`);
      } catch (e) {
        grok.shell.error(`Find domains failed: ${e instanceof Error ? e.message : e}`);
      } finally {
        pi.close();
      }
    })
    .show();
}
