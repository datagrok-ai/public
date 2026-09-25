import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {isMolBlock} from '../../utils/chem-common';
import {checkMoleculeValid} from '../../utils/chem-common-rdkit';
import {MAX_MATRIX_CELLS, MAX_MATRIX_COLS, MAX_MATRIX_ROWS, rankByFrequency}
  from './sar-matrix-assemble';
import {ClusterDecomposition, PositionRecord} from './sar-matrix-decompose';
import {attachmentNumbers, PositionFills} from './sar-matrix-link';
import {CoreCluster} from './sar-matrix-types';

/** Columns that already hold a decomposition, sorted into the axes they build. */
export interface SarFragmentColumns {
  core: DG.Column;
  rows: DG.Column[];
  column: DG.Column;
}

/** Distinct attachment points a column's fragments carry, ascending. */
function columnAttachments(column: DG.Column): number[] {
  const sites = new Set<number>();
  for (const value of column.categories) {
    for (const n of attachmentNumbers(value))
      sites.add(n);
  }
  return [...sites].sort((a, b) => a - b);
}

/**
 * Which R-group the matrix should enumerate across, by default.
 *
 * The one filling the LAST attachment runs along the top and the earlier ones fold into the row, so a
 * core with R1 and R2 opens as R1 down the side and R2 across. Which attachment a fragment fills
 * decides this, not the order the columns were picked in.
 */
export function defaultAxis(columns: DG.Column[]): DG.Column | undefined {
  const sites = columns.map((c) => columnAttachments(c)[0] ?? -1);
  return columns[sites.lastIndexOf(Math.max(...sites))];
}

/**
 * Canonicalizing reader for a fragment column, memoized over the distinct strings it holds.
 *
 * R-Group Analysis writes its core as molblocks whose coordinates differ per compound, so keying raw
 * text would make one scaffold read as one row per molecule. Attachment points are rewritten to the
 * `[*:n]` form the linker matches on, since a molblock round-trip can come back isotope-labelled.
 * Text RDKit cannot read is a label, kept as written — which is what lets component names be an axis.
 */
function fragmentReader(): {read: (column: DG.Column, i: number) => string, structures: Set<string>} {
  const canonical = new Map<string, string>();
  const structures = new Set<string>();
  const read = (column: DG.Column, i: number): string => {
    const text = column.isNone(i) ? '' : column.getString(i);
    if (text.trim() === '')
      return '';
    // A molblock opens with an empty title line, so trimming it shifts every header line up.
    const raw = isMolBlock(text) ? text : text.trim();
    const cached = canonical.get(raw);
    if (cached !== undefined)
      return cached;
    const mol = checkMoleculeValid(raw);
    let value = raw;
    if (mol?.is_valid()) {
      value = (mol.get_smiles() || raw).replace(/\[(\d+)\*\]/g, '[*:$1]');
      structures.add(value);
    }
    mol?.delete();
    canonical.set(raw, value);
    return value;
  };
  return {read, structures};
}

/** Records cut to what one matrix may hold, least-populated lines first so a trim keeps the data.
 *  Nothing else bounds this path: the size caps sit in the single-position assembler, and the gate on
 *  cluster size belongs to the decomposition that columns mode skips. */
function boundedRecords(records: PositionRecord[], spec: SarFragmentColumns): PositionRecord[] {
  const rowKey = (r: PositionRecord): string =>
    [r.coreSmiles, ...spec.rows.map((c) => r.values[c.name] ?? '')].join('\0');
  const rows = rankByFrequency(records.map(rowKey));
  const columns = rankByFrequency(records.map((r) => r.values[spec.column.name]));
  // Spent against what each axis is cut to: dividing both by the untrimmed other axis shrinks the
  // matrix quadratically.
  const maxColumns = Math.min(MAX_MATRIX_COLS, columns.length,
    Math.floor(MAX_MATRIX_CELLS / Math.min(MAX_MATRIX_ROWS, rows.length)));
  const maxRows = Math.min(MAX_MATRIX_ROWS, rows.length, Math.floor(MAX_MATRIX_CELLS / maxColumns));
  if (rows.length <= maxRows && columns.length <= maxColumns)
    return records;
  const keptRows = new Set(rows.slice(0, maxRows));
  const keptColumns = new Set(columns.slice(0, maxColumns));
  const bounded = records.filter((r) => keptRows.has(rowKey(r)) && keptColumns.has(r.values[spec.column.name]));
  grok.shell.warning(`SAR Matrix: the fragment columns describe ${rows.length} rows × ${columns.length} ` +
    `columns, which is more than one matrix can hold. It was cut to ${maxRows} × ${maxColumns}, ` +
    'keeping the rows and columns carrying the most compounds.');
  return bounded;
}

/**
 * Turn already-decomposed columns into the clusters and decompositions assembly consumes, skipping
 * fragmentation entirely.
 *
 * A series is the compounds sharing a core, so each distinct core value makes its own matrix — three
 * linkers are three series. Two things override that: a series column the user gave, which is their
 * own grouping and replaces this one, and having nothing on the row axis, where a matrix per core
 * would hold a single row and so be no matrix at all — there the cores are the rows instead.
 * A compound with no series value is left out, as it is when a series column groups fragmented ones.
 */
export function decomposeByColumns(spec: SarFragmentColumns, series: (string | null)[] | null,
  rowCount: number): {clusters: CoreCluster[], decomps: ClusterDecomposition[]} {
  const {read, structures} = fragmentReader();
  const positions = [spec.column.name, ...spec.rows.map((c) => c.name)];
  const splitByCore = series === null && spec.rows.length > 0;
  const byGroup = new Map<string, PositionRecord[]>();
  const claimed = new Set<string>();
  let duplicates = 0;

  for (let i = 0; i < rowCount; i++) {
    const value = series === null ? '' : series[i];
    if (value === null)
      continue;
    // No core means the scaffold never matched. A blank axis value is not that — it is the
    // unsubstituted parent, and it takes its own column.
    const axis = read(spec.column, i);
    const coreSmiles = read(spec.core, i);
    if (coreSmiles === '')
      continue;
    const group = splitByCore ? coreSmiles : value;
    const values: {[position: string]: string} = {[spec.column.name]: axis};
    for (const c of spec.rows)
      values[c.name] = read(c, i);
    // Assembly resolves a collision by overwriting, so replicates would vanish uncounted.
    const cell = [group, coreSmiles, ...spec.rows.map((c) => values[c.name]), axis].join('\0');
    if (claimed.has(cell)) {
      duplicates++;
      continue;
    }
    claimed.add(cell);
    if (!byGroup.has(group))
      byGroup.set(group, []);
    byGroup.get(group)!.push({molIdx: i, coreSmiles, values});
  }
  if (duplicates > 0) {
    const message = `SAR Matrix: ${duplicates} compounds share a core, row fragments and ` +
      `"${spec.column.name}" value with another compound and cannot have their own cell. The first ` +
      'of each set is shown — add the column that tells them apart to the row fragments.';
    _package.logger.warning(message);
    grok.shell.warning(message);
  }

  const clusters: CoreCluster[] = [];
  const decomps: ClusterDecomposition[] = [];
  /** Attachment points the cores carry that no picked column fills. Nothing can complete a proposal
   *  over one, so those cells come out with a value and no structure — which needs saying, or it
   *  reads as a bug. */
  const unfilled = new Set<number>();
  for (const [label, records] of byGroup) {
    const bounded = boundedRecords(records, spec);
    clusters.push({id: `f${clusters.length}`, series: [], siteKey: '', level: 2,
      label: series === null ? '' : label});
    const fills: PositionFills = {};
    for (const position of positions) {
      const numbers = new Set<number>();
      for (const record of bounded) {
        for (const n of attachmentNumbers(record.values[position] ?? ''))
          numbers.add(n);
      }
      fills[position] = [...numbers].sort((a, b) => a - b);
    }
    // Skipped where a position fills nothing: scanning against a label column would report the very
    // point it occupies.
    if (positions.every((p) => fills[p].length > 0)) {
      const covered = new Set(Object.values(fills).flat());
      for (const core of new Set(bounded.map((r) => r.coreSmiles))) {
        for (const n of attachmentNumbers(core)) {
          if (!covered.has(n))
            unfilled.add(n);
        }
      }
    }
    decomps.push({records: bounded, positions, links: {fills, structures}});
  }
  if (unfilled.size > 0) {
    const points = [...unfilled].sort((a, b) => a - b).map((n) => `[*:${n}]`).join(', ');
    const message = `SAR Matrix: the core carries ${points}, which none of the picked columns fills. ` +
      'Measured compounds are unaffected, but no predicted structure can be completed over an open ' +
      'attachment point, so predicted cells show a value and no structure. Add the column that fills ' +
      'it to the fragment columns.';
    _package.logger.warning(message);
    grok.shell.warning(message);
  }
  return {clusters, decomps};
}
