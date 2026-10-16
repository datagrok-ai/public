import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../../package';
import {isMolBlock} from '../../../utils/chem-common';
import {checkMoleculeValid} from '../../../utils/chem-common-rdkit';
import {linkStaged, MAX_MATRIX_CELLS, MAX_MATRIX_COLS, MAX_MATRIX_ROWS, rankByFrequency}
  from './sar-matrix-assemble';
import {ClusterDecomposition, PositionRecord} from './sar-matrix-decompose';
import {attachmentNumbers, fragmentLinks, planLink, PositionFills} from './sar-matrix-link';
import {CoreCluster} from '../sar-matrix-types';

/** Columns that already hold a decomposition, sorted into the axes they build. */
export interface SarFragmentColumns {
  core: DG.Column;
  rows: DG.Column[];
  column: DG.Column;
}

/** Numbered attachment points in every spelling in use: `[*:1]`, `[1*]`, `[1*:2]`, `[*1]`, `[R1]`, `[R:1]`. */
const R_LABEL = /\[(\d*)\*:(\d+)\]|\[(\d+)\*\]|\[\*(\d+)\]|\[R:?(\d+)\]/g;
const BARE_DUMMY = /(?<!\[)\*/g;
const ATOM_TOKEN = /\[[^\]]*\]|Br|Cl|[BCNOPSFIbcnops*]/g;

function relabel(smiles: string): string {
  return smiles.replace(R_LABEL, (_m, iso, map, iso2, star, r) => {
    const n = map !== undefined ? (Number(map) > 0 ? map : iso) : iso2 ?? star ?? r;
    return n ? `[*:${Number(n)}]` : '*';
  });
}

/** CXSMILES `*C |$_R1;$|` names its dummies in atom order. */
function applyCxLabels(text: string): string {
  const m = text.match(/^(\S+)\s+\|(.*)\|$/);
  if (m === null)
    return text;
  const labels = m[2].match(/\$([^$]*)\$/)?.[1].split(';') ?? [];
  let atom = 0;
  return m[1].replace(ATOM_TOKEN, (token) => {
    const label = (labels[atom++] ?? '').match(/^_?R(\d+)$/);
    return label !== null && (token === '*' || token === '[*]') ? `[*:${label[1]}]` : token;
  });
}

/**
 * A core or R-group as canonical SMILES with its attachment points written `[*:n]`, whatever notation
 * it came in: SMILES or CXSMILES labels, a molblock with `R#` atoms, or an isotope-labelled dummy.
 * Text RDKit cannot read is a label and comes back as written.
 */
export function standardizeFragment(text: string): {value: string, structure: boolean} {
  const molblock = isMolBlock(text);
  const written = molblock ? text : text.trim();
  const mol = checkMoleculeValid(molblock ? text : relabel(applyCxLabels(written)));
  let smiles: string;
  try {
    if (!mol?.is_valid())
      return {value: written, structure: false};
    smiles = mol.get_smiles();
  } finally {
    mol?.delete();
  }
  const labelled = relabel(smiles);
  // An isotope dummy orders the atoms differently from a numbered one, so the relabelled form is canonicalized again.
  return {value: labelled === smiles ? smiles : standardizeFragment(labelled).value, structure: true};
}

/** The one attachment number most of a column's fragments carry, else the number in its name. */
function columnSite(values: string[], name: string): number | null {
  const counts = new Map<number, number>();
  let numbered = 0;
  for (const value of values) {
    const numbers = attachmentNumbers(value);
    if (numbers.size > 0)
      numbered++;
    for (const n of numbers)
      counts.set(n, (counts.get(n) ?? 0) + 1);
  }
  for (const [n, count] of counts) {
    if (count * 2 > numbered)
      return n;
  }
  const fromName = name.match(/r[\s_\-:]*(\d+)/i);
  return fromName ? Number(fromName[1]) : null;
}

/** Numbers most of a position's fragments carry. */
function majoritySites(values: string[]): number[] {
  const counts = new Map<number, number>();
  let filled = 0;
  for (const value of values) {
    if (value === '')
      continue;
    filled++;
    for (const n of attachmentNumbers(value))
      counts.set(n, (counts.get(n) ?? 0) + 1);
  }
  return [...counts].filter(([, count]) => count * 2 > filled).map(([n]) => n).sort((a, b) => a - b);
}

/** Reads fragment columns in the standard notation, memoized over the distinct strings they hold. */
class FragmentReader {
  readonly structures = new Set<string>();
  private readonly cache = new Map<string, string>();

  read(column: DG.Column, i: number): string {
    const text = column.isNone(i) ? '' : column.getString(i);
    if (text.trim() === '')
      return '';
    let value = this.cache.get(text);
    if (value === undefined) {
      const standard = standardizeFragment(text);
      value = standard.value;
      if (standard.structure)
        this.structures.add(value);
      this.cache.set(text, value);
    }
    return value;
  }

  /** A column read whole, a lone unnumbered `*` numbered after the column's own site. */
  readColumn(column: DG.Column, rowCount: number, numberBareDummies: boolean): string[] {
    const values = Array.from({length: rowCount}, (_, i) => this.read(column, i));
    if (!numberBareDummies)
      return values;
    const site = columnSite(values, column.name);
    if (site === null)
      return values;
    const numbered = new Map<string, string>();
    return values.map((value) => {
      if (value === '' || attachmentNumbers(value).size > 0 || (value.match(BARE_DUMMY) ?? []).length !== 1)
        return value;
      let fixed = numbered.get(value);
      if (fixed === undefined) {
        fixed = standardizeFragment(value.replace(BARE_DUMMY, `[*:${site}]`)).value;
        this.structures.add(fixed);
        numbered.set(value, fixed);
      }
      return fixed;
    });
  }
}

/** Whether a column holds cores or R-groups: values carrying attachment points. */
export function holdsFragments(column: DG.Column): boolean {
  if (column.type !== DG.COLUMN_TYPE.STRING)
    return false;
  const labelled = /\[\d*\*|\*:\d|\[R:?\d|M {2}RGP|R#|\$_?R\d/;
  let seen = 0;
  for (const value of column.categories) {
    if (value === '')
      continue;
    if (labelled.test(value) || (!/\s/.test(value) && value.includes('*')))
      return true;
    if (++seen >= 20)
      return false;
  }
  return false;
}

/** The default matrix columns: R-groups rather than linkers, filled in the most compounds, then the
 *  highest attachment point. */
export function defaultAxis(columns: DG.Column[]): DG.Column | undefined {
  const scored = columns.map((column) => {
    const sample = column.categories.filter((v) => v !== '').slice(0, 20).map((v) => standardizeFragment(v).value);
    const linkers = sample.filter((v) => attachmentNumbers(v).size > 1).length;
    return {column, terminal: linkers * 2 <= sample.length ? 1 : 0,
      filled: column.length - column.stats.missingValueCount, site: columnSite(sample, column.name) ?? -1};
  });
  scored.sort((a, b) => b.terminal - a.terminal || b.filled - a.filled || b.site - a.site);
  return scored[0]?.column;
}

/** A series too large for one matrix, and how much of it the matrix shows. */
export interface SeriesCut {
  clusterId: string;
  rows: number;
  columns: number;
  keptRows: number;
  keptColumns: number;
  compounds: number;
  leftOut: number;
}

/** Records cut to what one matrix may hold, keeping the rows and columns with the most compounds. */
function boundedRecords(records: PositionRecord[], spec: SarFragmentColumns):
  {records: PositionRecord[], cut: Omit<SeriesCut, 'clusterId'> | null} {
  const rowKey = (r: PositionRecord): string =>
    [r.coreSmiles, ...spec.rows.map((c) => r.values[c.name] ?? '')].join('\0');
  const rows = rankByFrequency(records.map(rowKey));
  const columns = rankByFrequency(records.map((r) => r.values[spec.column.name]));
  // Spent against the trimmed other axis: dividing both by the untrimmed one shrinks the matrix twice.
  const maxColumns = Math.min(MAX_MATRIX_COLS, columns.length,
    Math.floor(MAX_MATRIX_CELLS / Math.min(MAX_MATRIX_ROWS, rows.length)));
  const maxRows = Math.min(MAX_MATRIX_ROWS, rows.length, Math.floor(MAX_MATRIX_CELLS / maxColumns));
  if (rows.length <= maxRows && columns.length <= maxColumns)
    return {records, cut: null};
  const keptRows = new Set(rows.slice(0, maxRows));
  const keptColumns = new Set(columns.slice(0, maxColumns));
  const kept = records.filter((r) => keptRows.has(rowKey(r)) && keptColumns.has(r.values[spec.column.name]));
  return {records: kept, cut: {rows: rows.length, columns: columns.length, keptRows: maxRows,
    keptColumns: maxColumns, compounds: records.length, leftOut: records.length - kept.length}};
}

/** The warning for series cut to fit one matrix, naming them as the navigator does. */
export function cutWarning(cuts: SeriesCut[], labelOf: (clusterId: string) => string | undefined): string | null {
  if (cuts.length === 0)
    return null;
  const n = (value: number): string => value.toLocaleString('en-US');
  const advice = 'Pick an R-group with fewer values for Matrix columns, or split the table with a Series column, ' +
    'to see them all.';
  if (cuts.length > 1) {
    const names = cuts.map((cut) => labelOf(cut.clusterId)).filter((name) => name);
    const leftOut = cuts.reduce((sum, cut) => sum + cut.leftOut, 0);
    return `SAR Matrix: ${cuts.length} series${names.length ? ` (${names.join(', ')})` : ''} are too large for ` +
      `one matrix, so they show only the rows and columns with the most compounds and leave out ${n(leftOut)} ` +
      `compounds. ${advice}`;
  }
  const cut = cuts[0];
  const shown = [cut.keptRows < cut.rows ? `${n(cut.keptRows)} rows` : '',
    cut.keptColumns < cut.columns ? `${n(cut.keptColumns)} columns` : ''].filter((part) => part).join(' and ');
  return `SAR Matrix: ${labelOf(cut.clusterId) ?? 'A series'} is too large for one matrix ` +
    `(${n(cut.rows)} × ${n(cut.columns)} R-group combinations), so it shows only the ${shown} with the most ` +
    `compounds and leaves out ${n(cut.leftOut)} of its ${n(cut.compounds)} compounds. ${advice}`;
}

/**
 * Clusters and decompositions read from columns that already hold a decomposition. Each distinct core
 * is a series unless a series column groups them; with no R-group on the rows, the cores are the rows.
 * Compounds without a core or a series value are left out.
 */
export function decomposeByColumns(spec: SarFragmentColumns, series: (string | null)[] | null,
  activities: Float32Array): {clusters: CoreCluster[], decomps: ClusterDecomposition[], cuts: SeriesCut[]} {
  const rowCount = activities.length;
  const reader = new FragmentReader();
  const coreValues = reader.readColumn(spec.core, rowCount, false);
  const axisValues = reader.readColumn(spec.column, rowCount, true);
  const rowValues = spec.rows.map((c) => reader.readColumn(c, rowCount, true));
  const positions = [spec.column.name, ...spec.rows.map((c) => c.name)];
  const splitByCore = series === null && spec.rows.length > 0;
  const byGroup = new Map<string, PositionRecord[]>();
  const cells = new Map<string, PositionRecord>();
  let duplicates = 0;

  for (let i = 0; i < rowCount; i++) {
    const value = series === null ? '' : series[i];
    const coreSmiles = coreValues[i];
    if (value === null || coreSmiles === '')
      continue;
    const group = splitByCore ? coreSmiles : value;
    const values: {[position: string]: string} = {[spec.column.name]: axisValues[i]};
    spec.rows.forEach((c, k) => values[c.name] = rowValues[k][i]);
    const key = [group, coreSmiles, ...positions.map((p) => values[p])].join('\0');
    const claimed = cells.get(key);
    if (claimed !== undefined) {
      duplicates++;
      if (!Number.isFinite(activities[claimed.molIdx]) && Number.isFinite(activities[i]))
        claimed.molIdx = i;
      continue;
    }
    const record = {molIdx: i, coreSmiles, values};
    cells.set(key, record);
    if (!byGroup.has(group))
      byGroup.set(group, []);
    byGroup.get(group)!.push(record);
  }
  if (duplicates > 0) {
    const message = `SAR Matrix: ${duplicates} compounds repeat the core and R-groups of another; one of ` +
      'each is shown, a measured one where there is one.';
    _package.logger.warning(message);
    grok.shell.warning(message);
  }

  const clusters: CoreCluster[] = [];
  const decomps: ClusterDecomposition[] = [];
  const cuts: SeriesCut[] = [];
  const unfilled = new Set<number>();
  for (const [label, records] of byGroup) {
    const {records: bounded, cut} = boundedRecords(records, spec);
    const clusterId = `f${clusters.length}`;
    if (cut !== null)
      cuts.push({clusterId, ...cut});
    clusters.push({id: clusterId, series: [], siteKey: '', level: 2, label: series === null ? '' : label});
    const fills: PositionFills = {};
    const sites: PositionFills = {};
    for (const position of positions) {
      const values = bounded.map((record) => record.values[position] ?? '');
      const numbers = new Set<number>();
      for (const value of values) {
        for (const n of attachmentNumbers(value))
          numbers.add(n);
      }
      fills[position] = [...numbers].sort((a, b) => a - b);
      sites[position] = majoritySites(values);
    }
    // A label column fills nothing, and would report the very point it occupies.
    if (positions.every((p) => fills[p].length > 0)) {
      const covered = new Set(Object.values(fills).flat());
      for (const core of new Set(bounded.map((r) => r.coreSmiles))) {
        for (const n of attachmentNumbers(core)) {
          if (!covered.has(n))
            unfilled.add(n);
        }
      }
    }
    decomps.push({records: bounded, positions, links: fragmentLinks(fills, sites, reader.structures)});
  }
  if (unfilled.size > 0) {
    const points = [...unfilled].sort((a, b) => a - b).map((n) => `[*:${n}]`).join(', ');
    const message = `SAR Matrix: no R-group column fills ${points} on the core, so predicted compounds ` +
      'there have no structure.';
    _package.logger.warning(message);
    grok.shell.warning(message);
  }
  return {clusters, decomps, cuts};
}

const CHECK_SAMPLE = 300;

function canonicalSmiles(molecule: string): string {
  const mol = checkMoleculeValid(molecule);
  try {
    return mol?.is_valid() ? mol.get_smiles() : '';
  } finally {
    mol?.delete();
  }
}

/** The largest component without stereo, canonical. */
function comparable(molecule: string): string {
  const smiles = canonicalSmiles(molecule).split('.').reduce((a, b) => b.length > a.length ? b : a, '');
  return smiles ? canonicalSmiles(smiles.replace(/@+/g, '').replace(/[/\\]/g, '')) : '';
}

/**
 * How many of a sample of compounds their own core and R-groups rebuild into something else, which
 * means the columns belong to another table or another molecule column. Stereo and counter-ions are
 * not compared: a decomposition may drop them.
 */
export async function checkAgainstMolecules(decomps: ClusterDecomposition[], molecules: string[]):
  Promise<{checked: number, mismatched: number}> {
  const all = decomps.flatMap((decomp) => decomp.records.map((record) => ({decomp, record})));
  const step = Math.max(1, Math.floor(all.length / CHECK_SAMPLE));
  const sample = all.filter((_t, i) => i % step === 0).slice(0, CHECK_SAMPLE);
  const plans = sample.map(({decomp, record}) =>
    planLink(record.coreSmiles, record.values, decomp.positions, decomp.links!, [], true));
  const built = await linkStaged(sample.map((t) => t.record.coreSmiles), plans,
    (i, position) => sample[i].record.values[position] ?? '');
  let checked = 0;
  let mismatched = 0;
  sample.forEach((t, i) => {
    const smiles = built[i];
    const molecule = molecules[t.record.molIdx];
    if (smiles === null || smiles === '' || smiles.includes('[*:') || !molecule)
      return;
    checked++;
    if (comparable(smiles) !== comparable(molecule))
      mismatched++;
  });
  return {checked, mismatched};
}
