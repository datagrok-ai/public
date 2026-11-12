/** Profile HMM search over a sequence column ("the table is the database"):
 * `hmmsearch` of one model against the column's non-empty rows (Z = their
 * number), split across workers with results identical to one search. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitFlags, type SearchOptions} from './hmmer/engine.ts';
import {resolveClanOverlaps} from './hmmer/clans.ts';
import {HmmerPool, type Progress} from './pool';
import {taskProgress} from './progress';
import {columnResidues} from './sequences';
import {type DomainHitSpan, setDomainAnnotations} from './annotations';
import {
  fetchPfam, loadBuiltIn, loadFile, loadText, type LoadedModels, modelLabel, parseAccessions,
} from './models';
import {CLAN_NOTE, HMMER_NOTE, withNotes} from './notice';

export interface ColumnSearch {
  /** Score, E-value, domain count, significance and domain spans per row. */
  columns: DG.Column[];
  /** Reported domains per row (0-based inclusive envelope residue indices). */
  domains: DomainHitSpan[][];
  /** Monomer position of each residue (for annotations). */
  positions: number[][];
  reported: number;
  included: number;
}

/** Search `column` with model `index` of `models`. With `clanCompetition` and a
 * Pfam family in a clan, overlapping domains of a sequence are resolved as
 * Pfam does (the most significant is kept; hmmer/clans.ts). */
export async function searchColumn(pool: HmmerPool, models: LoadedModels, index: number, column: DG.Column<string>,
  options: SearchOptions = {}, clanCompetition = true, label = models.info.models[index].name, progress?: Progress):
  Promise<ColumnSearch> {
  const model = models.info.models[index];
  const {residues: all, positions} = await columnResidues(column);
  const sequences: {name: string; residues: string}[] = [];
  all.forEach((residues, i) => {
    if (residues) sequences.push({name: String(i), residues});
  });
  const merged = await pool.search(models.key, index, sequences, options, undefined, progress);
  const n = column.length;
  const double = (name: string) => DG.Column.fromFloat64Array(name, new Float64Array(n).fill(DG.FLOAT_NULL));
  const score = double(`${label} score`);
  const evalue = double(`${label} E-value`);
  evalue.meta.format = 'scientific';
  const count = DG.Column.int(`${label} domains`, n);
  const significant = DG.Column.bool(`${label} significant`, n);
  const spans = DG.Column.string(`${label} regions`, n);
  const domains: DomainHitSpan[][] = Array.from({length: n}, () => []);
  for (const hit of merged.hits) {
    if (!(hit.flags & HitFlags.reported)) continue;
    const row = Number(sequences[hit.target].name);
    let shown = hit.domains.filter((d) => d.reported);
    const clan = clanCompetition ? models.clans?.[index] ?? null : null;
    if (clan) {
      shown = resolveClanOverlaps(shown.map((d, domain) => ({d, clan, aliFrom: d.aliFrom, aliTo: d.aliTo,
        evalue: d.printed.iEvalue, hitRank: 0, domain}))).map(({d}) => d);
    }
    score.set(row, hit.score, false);
    evalue.set(row, hit.evalue, false);
    count.set(row, shown.length, false);
    significant.set(row, (hit.flags & HitFlags.included) !== 0, false);
    spans.set(row, shown.map((d) => `${d.envFrom}-${d.envTo}`).join('; '), false);
    domains[row] = shown.map((d) => ({model: model.accession || model.name, from: d.envFrom - 1, to: d.envTo - 1,
      score: d.bitscore}));
  }
  return {columns: [score, evalue, count, significant, spans], domains, positions, reported: merged.reported,
    included: merged.included};
}

/** Add search results to the table and annotate the domains on the sequence column. */
export function applySearch(table: DG.DataFrame, column: DG.Column<string>, models: LoadedModels, index: number,
  result: ColumnSearch, description = modelLabel(models.info, index)): void {
  for (const col of result.columns) {
    col.name = table.columns.getUnusedName(col.name);
    table.columns.add(col);
  }
  const m = models.info.models[index];
  const label = {id: m.accession || m.name, name: m.name, description};
  setDomainAnnotations(table, column, [label],
    result.domains, result.positions);
}

/** Resolve a model source: 'builtin:<index>', Pfam accessions, or HMM text. */
export async function resolveModels(pool: HmmerPool, pkg: DG.Package, source: string): Promise<LoadedModels> {
  if (source === 'builtin' || source.startsWith('builtin:')) return loadBuiltIn(pool, pkg);
  if (/^\s*(PF\d{5}[\s,;]*)+$/i.test(source)) return fetchPfam(pool, parseAccessions(source));
  return loadText(pool, source);
}

/** `Bio | Search | Profile HMM Search...` */
export async function showProfileSearchDialog(pool: HmmerPool, pkg: DG.Package): Promise<void> {
  const df = grok.shell.tv?.dataFrame;
  const seqCols = df?.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE) ?? [];
  if (!df || seqCols.length === 0) {
    grok.shell.warning('Open a table with a Macromolecule column');
    return;
  }
  // Profiles built in this table (Build Profile HMM) are kept with it.
  const {builtProfiles, profileModels} = await import('./logo');
  const built = builtProfiles(df);
  const sources = ['Pfam library (built-in)', 'Pfam accession (online)', 'HMM file',
    ...(built.length > 0 ? ['Profile built in this table'] : [])];
  const columnInput = ui.input.column('Sequence', {table: df, value: seqCols[0],
    filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MACROMOLECULE});
  const sourceInput = ui.input.choice('Profile from', {value: sources[0], items: sources,
    tooltipText: withNotes('The profile HMM to search the column with (hmmsearch).', HMMER_NOTE)});
  const familyInput = ui.input.choice<string>('Family', {value: '', items: [''], nullable: false});
  const accessionInput = ui.input.string('Accession', {value: 'PF07686', tooltipText: 'Pfam accession, e.g. PF07686'});
  const fileInput = ui.input.file('HMM file');
  const builtInput = ui.input.choice('Profile', {value: built[built.length - 1] ?? '', items: built,
    tooltipText: 'A profile built from rows of this table (Build Profile HMM)'});
  const evalueInput = ui.input.float('E-value ≤', {value: 10, tooltipText: 'Report sequences up to this E-value (-E)'});
  const cutoffInput = ui.input.bool('Pfam gathering cutoff', {value: false,
    tooltipText: 'Use the model\'s GA bit score thresholds (--cut_ga) instead of E-values'});
  const clanInput = ui.input.bool('Clan competition', {value: true,
    tooltipText: withNotes('Of overlapping domains of a Pfam family in a clan, keep only the most significant ' +
      '(as Pfam does).', CLAN_NOTE)});
  let library: LoadedModels | null = null;
  const refresh = async () => {
    const source = sourceInput.value;
    familyInput.root.style.display = source === sources[0] ? '' : 'none';
    accessionInput.root.style.display = source === sources[1] ? '' : 'none';
    fileInput.root.style.display = source === sources[2] ? '' : 'none';
    builtInput.root.style.display = source === sources[3] ? '' : 'none';
    if (source === sources[0] && !library) {
      library = await loadBuiltIn(pool, pkg);
      const labels = library.info.models.map((_, i) => modelLabel(library!.info, i));
      familyInput.items = labels;
      familyInput.value = labels[0];
    }
  };
  sourceInput.onChanged.subscribe(() => refresh());
  await refresh();

  ui.dialog({title: 'Profile HMM Search'})
    .add(ui.inputs([columnInput, sourceInput, familyInput, accessionInput, fileInput, builtInput, evalueInput,
      cutoffInput, clanInput]))
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Profile HMM search...');
      try {
        let models: LoadedModels;
        let index = 0;
        if (sourceInput.value === sources[0]) {
          models = library ?? await loadBuiltIn(pool, pkg);
          index = Math.max(0, models.info.models.findIndex((_, i) => modelLabel(models.info, i) === familyInput.value));
        } else if (sourceInput.value === sources[1])
          models = await fetchPfam(pool, parseAccessions(accessionInput.value));
        else if (sourceInput.value === sources[3])
          models = (await profileModels(pool, pkg, builtInput.value ?? '', df)).models;
        else {
          if (!fileInput.value) throw new Error('Choose an HMM file');
          models = await loadFile(pool, fileInput.value);
        }
        const options: SearchOptions = cutoffInput.value ? {cutoff: 'ga'} : {evalue: evalueInput.value ?? 10};
        const column = columnInput.value as DG.Column<string>;
        const result = await searchColumn(pool, models, index, column, options, clanInput.value, undefined,
          taskProgress(pi, 'Searching'));
        applySearch(df, column, models, index, result);
        grok.shell.info(`${models.info.models[index].name}: ${result.reported} sequences found, ` +
          `${result.included} significant`);
      } catch (e) {
        grok.shell.error(`Profile HMM search failed: ${e instanceof Error ? e.message : e}`);
      } finally {
        pi.close();
      }
    })
    .show();
}
