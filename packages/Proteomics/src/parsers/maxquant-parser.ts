import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {log2TransformColumns, addPrimaryColumnIfNeeded} from './shared-utils';
import {resolveGeneLabels} from '../utils/gene-label-resolver';

/** Intensity column prefixes to detect, checked in this order.
 * "intensity" must be last since it is a substring of "lfq intensity". */
const INTENSITY_PREFIXES = ['lfq intensity', 'ibaq', 'reporter intensity', 'intensity'];

/** Column name variants for filter columns across MaxQuant versions. */
const FILTER_COLUMNS = {
  contaminant: ['Potential contaminant', 'Potential.contaminant'],
  reverse: ['Reverse'],
  onlyBySite: ['Only identified by site', 'Only.identified.by.site'],
} as const;

/** Column name variants for protein IDs. */
const PROTEIN_ID_COLUMNS = ['Protein IDs', 'Majority protein IDs'];

/** Column name variants for gene names. */
const GENE_NAME_COLUMNS = ['Gene names', 'Gene name'];

/** Finds the first existing column from a list of name candidates. */
function findCol(df: DG.DataFrame, names: readonly string[]): DG.Column | null {
  for (const name of names) {
    const col = df.col(name);
    if (col)
      return col;
  }
  return null;
}

/** Marks rows for exclusion where the given column has value "+". */
function filterByMarker(df: DG.DataFrame, colNames: readonly string[]): void {
  const col = findCol(df, colNames);
  if (!col) return;
  for (let i = 0; i < df.rowCount; i++) {
    if (col.get(i) === '+')
      df.filter.set(i, false, false);
  }
}

/** Marks rows for exclusion when any protein-id column starts with CON__/REV__.
 * Scans every column in PROTEIN_ID_COLUMNS rather than the first match — MaxQuant
 * sometimes drops the prefix from `Protein IDs` but retains it in `Majority protein IDs`. */
function filterByIdPrefix(df: DG.DataFrame): void {
  const cols: DG.Column[] = [];
  for (const name of PROTEIN_ID_COLUMNS) {
    const c = df.col(name);
    if (c) cols.push(c);
  }
  if (cols.length === 0) return;
  for (let i = 0; i < df.rowCount; i++) {
    for (const col of cols) {
      const val = col.get(i);
      if (typeof val === 'string' && (val.startsWith('CON__') || val.startsWith('REV__'))) {
        df.filter.set(i, false, false);
        break;
      }
    }
  }
}

/** Finds the source column for primary extraction from MaxQuant-specific name variants,
 * then delegates to shared addPrimaryColumnIfNeeded which handles semicolon detection. */
function addPrimaryColumnFromVariants(df: DG.DataFrame, sourceNames: readonly string[],
  newName: string, semType: string): void {
  const srcCol = findCol(df, sourceNames);
  if (!srcCol) return;
  addPrimaryColumnIfNeeded(df, srcCol.name, newName, semType);
}

/** Scans header row to find intensity column names and build columnImportOptions
 * that force them to double. Prevents fromCsv from mistyping large-integer
 * intensity columns as boolean. */
function buildIntensityColumnOptions(text: string): DG.CsvImportColumnOptions[] {
  const newlineIdx = text.indexOf('\n');
  const firstLine = (newlineIdx >= 0 ? text.substring(0, newlineIdx) : text).replace(/\r$/, '');
  const headers = firstLine.split('\t');
  const opts: DG.CsvImportColumnOptions[] = [];
  for (const header of headers) {
    const lower = header.toLowerCase();
    for (const prefix of INTENSITY_PREFIXES) {
      if (lower.startsWith(prefix)) {
        opts.push({name: header, type: 'double'});
        break;
      }
    }
  }
  return opts;
}

/** Detects intensity columns by MaxQuant-specific prefix matching and
 * delegates to shared log2TransformColumns for the actual transformation. */
function processIntensityColumns(df: DG.DataFrame): void {
  const intensityNames: string[] = [];
  for (const name of df.columns.names()) {
    const col = df.col(name);
    if (!col) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT)
      continue;
    const lowerName = name.toLowerCase();
    for (const prefix of INTENSITY_PREFIXES) {
      if (lowerName.startsWith(prefix)) {
        intensityNames.push(name);
        break;
      }
    }
  }
  log2TransformColumns(df, intensityNames);
}

/** Assigns semantic types to protein ID and gene name columns. */
function assignSemanticTypes(df: DG.DataFrame): void {
  const proteinCol = findCol(df, PROTEIN_ID_COLUMNS);
  if (proteinCol)
    proteinCol.semType = SEMTYPE.PROTEIN_ID;

  const geneCol = findCol(df, GENE_NAME_COLUMNS);
  if (geneCol)
    geneCol.semType = SEMTYPE.GENE_SYMBOL;
}

/** Parses MaxQuant proteinGroups.txt TSV text into a filtered, typed DataFrame.
 *
 * Applies standard proteomics filtering:
 * - Removes contaminant rows ('+' marker or CON__ prefix)
 * - Removes reverse hits ('+' marker or REV__ prefix)
 * - Removes only-identified-by-site rows
 *
 * Adds log2-transformed intensity columns and semantic type annotations. */
export async function parseMaxQuantText(text: string): Promise<DG.DataFrame> {
  const columnImportOptions = buildIntensityColumnOptions(text);
  const raw = DG.DataFrame.fromCsv(text, {delimiter: '\t', columnImportOptions});

  // Apply row filters
  raw.filter.init((_i) => true);
  filterByMarker(raw, FILTER_COLUMNS.contaminant);
  filterByMarker(raw, FILTER_COLUMNS.reverse);
  filterByMarker(raw, FILTER_COLUMNS.onlyBySite);
  filterByIdPrefix(raw);
  raw.filter.fireChanged();

  // Clone filtered rows
  const df = raw.clone(raw.filter);

  // Assign semantic types to ID and gene columns
  assignSemanticTypes(df);

  // Parse semicolon-delimited fields into primary columns (only if semicolons present)
  addPrimaryColumnFromVariants(df, PROTEIN_ID_COLUMNS, 'Primary Protein ID', SEMTYPE.PROTEIN_ID);
  addPrimaryColumnFromVariants(df, GENE_NAME_COLUMNS, 'Primary Gene Name', SEMTYPE.GENE_SYMBOL);

  // Detect and log2-transform intensity columns
  processIntensityColumns(df);

  df.setTag('proteomics.source', 'maxquant');

  await resolveGeneLabels(df);

  return df;
}
