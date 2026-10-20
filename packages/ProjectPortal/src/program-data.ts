/** A program's governed data against the platform (F-09): the portfolio's queries run for one
 * program, as a dataframe named after it, opened as a table view with its viewers. The query
 * runs as the user, on a connection shared with them, so the data stays governed by the
 * connection's grants; the program itself is read under row security first. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {portfolioByCode, queryFrame, sqlLiteral} from './portal';
import {
  ColumnInfo, datasetTableName, DatasetViewer, defaultViewers, MAIN_DATASET, ProgramDataset, ProgramDatasetConfig,
  programDatasets,
} from './model/program-data';
import {ProgramDataConfig} from './model/conditions';
import {filterValue, ROW_LIMIT} from './rows';

/** The datasets a portfolio offers for each of its programs, the program dataset first. */
export function datasetsOf(portfolio: PortfolioRow): ProgramDataset[] {
  return programDatasets(portfolio.program_data as ProgramDataConfig | undefined,
    (portfolio.program_datasets as {datasets?: ProgramDatasetConfig[]} | undefined)?.datasets);
}

async function programOf(portfolio: PortfolioRow, code: string): Promise<ProgramRow> {
  const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(code)}`, limit: ROW_LIMIT}))[0];
  if (!p)
    throw new Error(`${portfolio.code} has no program ${code}, or it is not shared with you`);
  return p;
}

export interface ProgramTable {
  df: DG.DataFrame;
  dataset: ProgramDataset;
  program: ProgramRow;
  portfolio: PortfolioRow;
}

/** One program's rows of a dataset (default: the program dataset), named after the program. */
export async function programTable(portfolioCode: string, programCode: string, dataset?: string): Promise<ProgramTable> {
  const portfolio = await portfolioByCode(portfolioCode);
  const program = await programOf(portfolio, programCode);
  const all = datasetsOf(portfolio);
  const code = dataset || MAIN_DATASET;
  const d = all.find((x) => x.code === code);
  if (!d) {
    throw new Error(all.length === 0 ? `${portfolio.code} declares no program data` :
      `${portfolio.code} has no dataset ${code}; it has ${all.map((x) => x.code).join(', ')}`);
  }
  const df = await queryFrame(d.connection, d.sql.replace(/@program\b/g, () => sqlLiteral(program.code)));
  df.name = datasetTableName(program, d);
  df.setTag('portal.portfolio', portfolio.code);
  df.setTag('portal.program', program.code);
  df.setTag('portal.dataset', d.code);
  return {df, dataset: d, program, portfolio};
}

const columnInfo = (df: DG.DataFrame): ColumnInfo[] => df.columns.toList().map((c) => ({name: c.name,
  kind: c.isNumerical ? 'number' : c.type === DG.COLUMN_TYPE.STRING ? 'string' : 'other',
  categories: c.type === DG.COLUMN_TYPE.STRING ? c.categories.length : undefined}));

/** The viewers a dataset opens with: its own, else the defaults for the portfolio's setting. A
 * viewer naming a column the table lacks is skipped, not an error. */
export function viewersFor(t: ProgramTable): DatasetViewer[] {
  const viewers = t.dataset.viewers?.length ? t.dataset.viewers :
    defaultViewers(t.portfolio.setting, t.portfolio.modality_class, columnInfo(t.df));
  const names = new Set(t.df.columns.names());
  return viewers.filter((v) => Object.entries(v.options ?? {})
    .every(([k, val]) => !k.endsWith('ColumnName') || typeof val !== 'string' || names.has(val)));
}

/** Opens the dataset as a table view with its viewers. */
export async function openProgramTable(portfolioCode: string, programCode: string, dataset?: string): Promise<DG.TableView> {
  const t = await programTable(portfolioCode, programCode, dataset);
  const tv = grok.shell.addTableView(t.df);
  for (const v of viewersFor(t)) {
    try {
      tv.addViewer(v.type, v.options ?? {});
    } catch (e: any) {
      grok.log.warning(`Project Portal: viewer ${v.type} not added to ${t.df.name}: ${e?.message ?? e}`);
    }
  }
  return tv;
}
