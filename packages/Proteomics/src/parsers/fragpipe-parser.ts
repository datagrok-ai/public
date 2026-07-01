import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {log2TransformColumns, addPrimaryColumnIfNeeded} from './shared-utils';
import {resolveGeneLabels} from '../utils/gene-label-resolver';

/** Per-sample intensity column suffixes in FragPipe's combined_protein.tsv,
 * in order of quantitative preference. MaxLFQ is the FragPipe-recommended
 * LFQ quant; bare " Intensity" is the summed peptide intensity. Razor is
 * the razor-protein attribution and not normally used as the primary quant. */
const INTENSITY_SUFFIXES = [' MaxLFQ Intensity', ' Intensity', ' Razor Intensity'] as const;

/** Column name variants for protein-level identifiers in FragPipe output.
 * `Protein ID` is the UniProt accession; `Protein` is the full FASTA header
 * (e.g., `sp|P12345|GENE_HUMAN`) which we prefer only if `Protein ID` is missing. */
const PROTEIN_ID_COLUMNS = ['Protein ID', 'Protein'];

/** Column name variants for gene symbols. */
const GENE_NAME_COLUMNS = ['Gene', 'Gene Names', 'Gene Symbol'];

/** Decoy / contaminant prefixes to filter out at the row level.
 * `contam_` is FragPipe/Philosopher's convention; `rev_` is Philosopher's decoy
 * prefix (usually filtered upstream, but defensive); CON__/REV__ cover the
 * case where the FASTA was MaxQuant-formatted. */
const CONTAMINANT_PREFIXES = ['contam_', 'rev_', 'CON__', 'REV__'];

/** Finds the first existing column from a list of name candidates. */
function findCol(df: DG.DataFrame, names: readonly string[]): DG.Column | null {
  for (const name of names) {
    const col = df.col(name);
    if (col)
      return col;
  }
  return null;
}

/** Checks whether a header is a per-sample intensity column.
 * Returns the matched suffix (lowercased) or null. */
function matchIntensitySuffix(header: string): string | null {
  const lower = header.toLowerCase();
  for (const suffix of INTENSITY_SUFFIXES) {
    if (lower.endsWith(suffix.toLowerCase()) && lower.length > suffix.length)
      return suffix.toLowerCase();
  }
  return null;
}

/** Marks rows for exclusion when their primary protein identifier starts with
 * any contaminant/decoy prefix. Checks both `Protein ID` and `Protein` because
 * Philosopher contam_/rev_ prefixes appear on whichever column is populated. */
function filterContaminantRows(df: DG.DataFrame): void {
  const candidates = ['Protein ID', 'Protein'];
  const cols: DG.Column[] = [];
  for (const name of candidates) {
    const c = df.col(name);
    if (c) cols.push(c);
  }
  if (cols.length === 0) return;

  for (let i = 0; i < df.rowCount; i++) {
    for (const col of cols) {
      const val = col.get(i);
      if (typeof val !== 'string') continue;
      let matched = false;
      for (const prefix of CONTAMINANT_PREFIXES) {
        if (val.startsWith(prefix)) {
          matched = true;
          break;
        }
      }
      if (matched) {
        df.filter.set(i, false, false);
        break;
      }
    }
  }
}

/** Scans the header row to find intensity columns and forces them to `double`
 * during CSV import. Without this, fromCsv can mistype large-integer intensities
 * (or all-zero / all-blank samples) as boolean or int. Mirrors the same guard
 * used by the MaxQuant parser. */
function buildIntensityColumnOptions(text: string): DG.CsvImportColumnOptions[] {
  const newlineIdx = text.indexOf('\n');
  const firstLine = (newlineIdx >= 0 ? text.substring(0, newlineIdx) : text).replace(/\r$/, '');
  const headers = firstLine.split('\t');
  const opts: DG.CsvImportColumnOptions[] = [];
  for (const header of headers) {
    if (matchIntensitySuffix(header))
      opts.push({name: header, type: 'double'});
  }
  return opts;
}

/** Returns the set of intensity columns to log2-transform. If MaxLFQ Intensity
 * columns exist, prefer them and drop the redundant bare " Intensity" / " Razor
 * Intensity" columns for the same sample. If MaxLFQ is absent (e.g. spectral-
 * count-only run, or LFQ not enabled), fall back to " Intensity". Razor is only
 * used as a last resort. Returned columns still get SEMTYPE.INTENSITY below. */
function pickQuantColumns(df: DG.DataFrame): string[] {
  type Bucket = {maxLfq?: string; intensity?: string; razor?: string};
  const samples = new Map<string, Bucket>();

  for (const name of df.columns.names()) {
    const col = df.col(name);
    if (!col) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT &&
        col.type !== DG.COLUMN_TYPE.BIG_INT)
      continue;
    const lower = name.toLowerCase();
    let sample: string | null = null;
    let kind: keyof Bucket | null = null;
    if (lower.endsWith(' maxlfq intensity')) {
      sample = name.substring(0, name.length - ' MaxLFQ Intensity'.length);
      kind = 'maxLfq';
    } else if (lower.endsWith(' razor intensity')) {
      sample = name.substring(0, name.length - ' Razor Intensity'.length);
      kind = 'razor';
    } else if (lower.endsWith(' intensity')) {
      sample = name.substring(0, name.length - ' Intensity'.length);
      kind = 'intensity';
    }
    if (!sample || !kind) continue;
    if (!samples.has(sample)) samples.set(sample, {});
    samples.get(sample)![kind] = name;
  }

  const chosen: string[] = [];
  for (const bucket of samples.values()) {
    if (bucket.maxLfq) chosen.push(bucket.maxLfq);
    else if (bucket.intensity) chosen.push(bucket.intensity);
    else if (bucket.razor) chosen.push(bucket.razor);
  }
  return chosen;
}

/** Assigns semantic types to ID, gene, and the chosen per-sample intensity columns.
 * Only `chosenQuant` columns receive SEMTYPE.INTENSITY — the redundant Razor /
 * bare-Intensity columns we dropped from log2 transform must not surface as
 * intensities in downstream pickers. */
function assignSemanticTypes(df: DG.DataFrame, chosenQuant: string[]): void {
  const proteinCol = findCol(df, PROTEIN_ID_COLUMNS);
  if (proteinCol)
    proteinCol.semType = SEMTYPE.PROTEIN_ID;

  const geneCol = findCol(df, GENE_NAME_COLUMNS);
  if (geneCol)
    geneCol.semType = SEMTYPE.GENE_SYMBOL;

  for (const name of chosenQuant) {
    const col = df.col(name);
    if (col) col.semType = SEMTYPE.INTENSITY;
  }
}

/** Finds the source column for primary extraction from FragPipe-specific name
 * variants and delegates to shared addPrimaryColumnIfNeeded. */
function addPrimaryColumnFromVariants(df: DG.DataFrame, sourceNames: readonly string[],
  newName: string, semType: string): void {
  const srcCol = findCol(df, sourceNames);
  if (!srcCol) return;
  addPrimaryColumnIfNeeded(df, srcCol.name, newName, semType);
}

/** Parses a FragPipe `combined_protein.tsv` into a filtered, typed DataFrame.
 *
 * Applies:
 * - Contaminant/decoy row filtering (`contam_`, `rev_`, plus MaxQuant-style
 *   CON__/REV__ as a fallback)
 * - Semantic type assignment (PROTEIN_ID, GENE_SYMBOL, INTENSITY)
 * - Primary protein/gene split for semicolon-delimited identifier lists
 * - Log2 transform of the preferred quantitative intensity columns
 *   (MaxLFQ Intensity > Intensity > Razor Intensity, picked per sample)
 *
 * The returned DataFrame is ready for the standard pipeline: Annotate
 * Experiment -> Normalize -> Impute -> Differential Expression -> Viewers. */
export async function parseFragPipeText(text: string): Promise<DG.DataFrame> {
  const columnImportOptions = buildIntensityColumnOptions(text);
  const raw = DG.DataFrame.fromCsv(text, {delimiter: '\t', columnImportOptions});

  raw.filter.init((_i) => true);
  filterContaminantRows(raw);
  raw.filter.fireChanged();

  const df = raw.clone(raw.filter);

  // Choose one quant column per sample for log2 transform.
  const chosen = pickQuantColumns(df);
  assignSemanticTypes(df, chosen);

  addPrimaryColumnFromVariants(df, PROTEIN_ID_COLUMNS, 'Primary Protein ID', SEMTYPE.PROTEIN_ID);
  addPrimaryColumnFromVariants(df, GENE_NAME_COLUMNS, 'Primary Gene Name', SEMTYPE.GENE_SYMBOL);

  log2TransformColumns(df, chosen);

  df.setTag('proteomics.source', 'fragpipe');

  await resolveGeneLabels(df);

  return df;
}
