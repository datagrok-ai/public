/** A program's governed data one click away (F-09): the portfolio's program dataset and any extra
 * datasets it names, what each is called, and which viewers open with it. Pure; the program page,
 * openProgramData and programData all read the same list. */
import type {Setting} from './config';

/** The viewer types a dataset may name; the platform's own names. */
export const DATASET_VIEWERS = ['Scatter plot', 'Histogram', 'Bar chart', 'Line chart', 'Box plot', 'Pie chart',
  'Statistics', 'Correlation plot'] as const;
export type DatasetViewerType = typeof DATASET_VIEWERS[number];

/** A viewer opened with a dataset: its type and the platform's viewer properties
 * (xColumnName, yColumnName, colorColumnName, valueColumnName, splitColumnName, ...). */
export interface DatasetViewer {
  type: DatasetViewerType;
  options?: Record<string, unknown>;
}

/** A further dataset of every program of a portfolio: a SAR table beside the process dataset,
 * assay results beside the compound table. One SQL on a named connection, `@program` for the code. */
export interface ProgramDatasetConfig {
  code: string;
  name: string;
  connection: string;
  sql: string;
  description?: string;
  /** Viewers to open with it; default: by the portfolio's setting and the columns. */
  viewers?: DatasetViewer[];
}

/** One dataset as the program page offers it. `program` is the program dataset itself. */
export interface ProgramDataset {
  code: string;
  name: string;
  connection: string;
  sql: string;
  description?: string;
  viewers?: DatasetViewer[];
}

export const MAIN_DATASET = 'program';
const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

/** The datasets of a portfolio in order: the program dataset first (when declared), then the extras. */
export function programDatasets(main: {connection: string, sql: string, name?: string, viewers?: DatasetViewer[]} | null | undefined,
  extras: ProgramDatasetConfig[] | null | undefined): ProgramDataset[] {
  const out: ProgramDataset[] = [];
  if (main?.sql)
    out.push({code: MAIN_DATASET, name: main.name || 'Program data', connection: main.connection, sql: main.sql,
      viewers: main.viewers});
  return [...out, ...(extras ?? [])];
}

/** Problems with a portfolio's extra datasets and the viewers it names, as sentences. */
export function validateProgramDatasets(extras: ProgramDatasetConfig[] | undefined,
  mainViewers?: DatasetViewer[]): string[] {
  const errors: string[] = [];
  const viewerErrors = (where: string, viewers: DatasetViewer[] | undefined) => {
    for (const v of viewers ?? []) {
      if (!(DATASET_VIEWERS as readonly string[]).includes(v?.type))
        errors.push(`${where}: viewer "${v?.type}" is not one of ${DATASET_VIEWERS.join(', ')}`);
    }
  };
  viewerErrors('programData', mainViewers);
  const codes = new Set<string>([MAIN_DATASET]);
  for (const d of extras ?? []) {
    const where = `programDatasets ${d.code ?? '?'}`;
    if (!d.code || !CODE.test(d.code))
      errors.push(`${where}: code is required and may hold only letters, digits, '_', '.', '-'`);
    else if (codes.has(d.code))
      errors.push(`${where}: code is used twice${d.code === MAIN_DATASET ? ` ("${MAIN_DATASET}" is the program dataset)` : ''}`);
    codes.add(d.code);
    if (!d.name)
      errors.push(`${where}: name is required`);
    if (!d.connection)
      errors.push(`${where}: connection is required`);
    if (!d.sql)
      errors.push(`${where}: sql is required`);
    else if (!d.sql.includes('@program'))
      errors.push(`${where}: sql must select one program's rows with @program`);
    viewerErrors(where, d.viewers);
  }
  return errors;
}

/** The table's name: the program, and the dataset when it is not the program dataset. */
export function datasetTableName(program: {code: string, name?: string | null}, dataset: {code: string, name: string}): string {
  const label = program.name && program.name !== program.code ? `${program.code} · ${program.name}` : program.code;
  return dataset.code === MAIN_DATASET ? label : `${label} · ${dataset.name}`;
}

export interface ColumnInfo {
  name: string;
  /** 'number' for int, double, bigint; 'string' for string; anything else is left out. */
  kind: 'number' | 'string' | 'other';
  /** Distinct values, for the categorical columns a split or colour can use. */
  categories?: number;
}

/** Viewers to open with a dataset that names none. Discovery looks at structure-activity: the
 * first two measures against each other, coloured by the first category, and the first measure's
 * distribution. Development and rollout look at counts: a bar chart over the first category, and
 * the first measure's distribution. The program column, the same on every row, is never used. */
export function defaultViewers(setting: Setting | string | null | undefined, _modality: string | null | undefined,
  columns: ColumnInfo[]): DatasetViewer[] {
  const usable = columns.filter((c) => c.name.toLowerCase() !== 'program');
  const numbers = usable.filter((c) => c.kind === 'number');
  const categories = usable.filter((c) => c.kind === 'string' && (c.categories ?? 0) > 1 && (c.categories ?? 0) <= 20);
  const out: DatasetViewer[] = [];
  if (setting === 'discovery') {
    if (numbers.length >= 2) {
      out.push({type: 'Scatter plot', options: {xColumnName: numbers[0].name, yColumnName: numbers[1].name,
        ...(categories[0] ? {colorColumnName: categories[0].name} : {})}});
    }
  } else if (categories[0])
    out.push({type: 'Bar chart', options: {splitColumnName: categories[0].name}});
  if (numbers[0])
    out.push({type: 'Histogram', options: {valueColumnName: numbers[0].name}});
  return out;
}
