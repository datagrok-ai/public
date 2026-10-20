/** Compounds needing attention against the platform (F-17): the portfolio's rules evaluated over one
 * program's dataset, as a table for scripts, and a flagged compound opened in the program's data.
 * The rules live in ./model/attention; the query runs as the user, like every program dataset. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {portfolioByCode, programData} from './portal';
import {AttentionRule, evaluateAttention, ProgramAttention} from './model/attention';
import {DataRow, DEFAULT_LEAD_ROLES, ProgramDataConfig, rowId} from './model/conditions';
import {openProgramTable} from './program-data';
import {filterValue, ROW_LIMIT} from './rows';

/** The attention rules a portfolio declares; none when absent. */
export function attentionRulesOf(portfolio: PortfolioRow): AttentionRule[] {
  return (portfolio.attention as {rules?: AttentionRule[]} | undefined)?.rules ?? [];
}

/** The rules of [portfolio] over [rows]; null when it declares no rules. */
export function attentionFor(portfolio: PortfolioRow, rows: DataRow[] | null): ProgramAttention | null {
  const rules = attentionRulesOf(portfolio);
  if (rules.length === 0)
    return null;
  const leadRoles = (portfolio.program_data as ProgramDataConfig | undefined)?.leadRoles ?? DEFAULT_LEAD_ROLES;
  return evaluateAttention(rules, rows ?? [], leadRoles);
}

async function programOf(portfolio: PortfolioRow, code: string): Promise<ProgramRow> {
  const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(code)}`, limit: ROW_LIMIT}))[0];
  if (!p)
    throw new Error(`${portfolio.code} has no program ${code}, or it is not shared with you`);
  return p;
}

/** One program's attention: the program is read under row security first, then its dataset. */
export async function programAttentionOf(portfolioCode: string, programCode: string):
  Promise<{portfolio: PortfolioRow, program: ProgramRow, attention: ProgramAttention}> {
  const portfolio = await portfolioByCode(portfolioCode);
  if (attentionRulesOf(portfolio).length === 0)
    throw new Error(`${portfolio.code} declares no attention rules`);
  const program = await programOf(portfolio, programCode);
  const rows = await programData(portfolio, program.code);
  return {portfolio, program, attention: attentionFor(portfolio, rows)!};
}

/** One row per item and rule it trips: rule, rule name, severity, item, message. */
export function attentionFrame(a: ProgramAttention, rules: AttentionRule[], name: string): DG.DataFrame {
  const nameOf = new Map(rules.map((r) => [r.code, r.name]));
  const S = DG.COLUMN_TYPE.STRING;
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(S, 'item', a.flags.map((f) => f.item)),
    DG.Column.fromList(S, 'rule', a.flags.map((f) => f.rule)),
    DG.Column.fromList(S, 'rule name', a.flags.map((f) => nameOf.get(f.rule) ?? f.rule)),
    DG.Column.fromList(S, 'severity', a.flags.map((f) => f.severity)),
    DG.Column.fromList(S, 'message', a.flags.map((f) => f.message)),
  ]);
  df.name = name;
  return df;
}

/** Opens the program's dataset with [item]'s row current and every flagged row selected. */
export async function openAttentionItem(portfolioCode: string, programCode: string, item: string,
  flagged: string[] = [item]): Promise<DG.TableView> {
  const tv = await openProgramTable(portfolioCode, programCode);
  const df = tv.dataFrame;
  const cols = df.columns.toList();
  const idOf = (i: number) => rowId(Object.fromEntries(cols.map((c) => [c.name, c.isNone(i) ? null : c.get(i)])), i);
  const wanted = new Set(flagged);
  let current = -1;
  for (let i = 0; i < df.rowCount; i++) {
    const id = idOf(i);
    if (wanted.has(id))
      df.selection.set(i, true, false);
    if (id === item && current < 0)
      current = i;
  }
  df.selection.fireChanged();
  if (current >= 0)
    df.currentRowIdx = current;
  else
    grok.shell.warning(`${item} is not in ${programCode}'s data any more`);
  return tv;
}
