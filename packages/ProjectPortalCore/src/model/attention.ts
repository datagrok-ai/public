/** Compounds needing attention (F-17, FR-104): rules a portfolio declares over its program dataset,
 * in the same condition language gate criteria use (scope, where, measured). A gate condition
 * counts the rows that qualify; an attention rule names the rows that trip it. Pure: the caller
 * runs the program dataset query and hands the rows in. */
import {
  columnValue, ConditionScope, DataCondition, DataRow, DEFAULT_LEAD_ROLES, holds, isBlank, rowId, rowsInScope,
  validateCondition,
} from './conditions';

export type AttentionSeverity = 'critical' | 'warning' | 'info';
export const ATTENTION_SEVERITIES: AttentionSeverity[] = ['critical', 'warning', 'info'];

export interface AttentionRule {
  code: string;
  name: string;
  /** Default `warning`. */
  severity?: AttentionSeverity;
  /** A row trips the rule when the scope admits it, every `where` predicate holds and, when the
   * rule names `measured` columns, it lacks at least one of them. `minCount` does not apply. */
  condition: DataCondition;
  /** One line per flagged row. `{column}` stands for the row's value, `{missing}` for the measured
   * columns it lacks. Default: the rule's name, the values it tested and what is missing. */
  message?: string;
}

/** One row that trips one rule. */
export interface AttentionFlag {
  rule: string;
  severity: AttentionSeverity;
  /** The row's compound, id, code, lot or name. */
  item: string;
  message: string;
  /** The values the rule tested: its where columns and its measured columns. */
  values: Record<string, unknown>;
  /** Measured columns the row lacks. */
  missing: string[];
}

export interface RuleAttention {
  code: string;
  name: string;
  severity: AttentionSeverity;
  /** Distinct items that trip it. */
  items: string[];
}

export interface AttentionItem {
  item: string;
  /** The rules it trips, most severe first. */
  rules: string[];
  severity: AttentionSeverity;
  messages: string[];
}

export interface ProgramAttention {
  /** Every rule, in declared order, with what trips it (none is a count of 0). */
  rules: RuleAttention[];
  /** Every flagged item once, most severe first, then by name. */
  items: AttentionItem[];
  flags: AttentionFlag[];
  /** Distinct flagged items. */
  count: number;
  worst: AttentionSeverity | null;
  summary: string;
}

const rank = (s: AttentionSeverity) => ATTENTION_SEVERITIES.indexOf(s);
const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

const show = (v: unknown): string => isBlank(v) ? '—' : typeof v === 'number' ? String(+v.toPrecision(4)) : String(v);

/** The rule's line for one row: its template filled in, or name, tested values and what is missing. */
export function attentionMessage(rule: AttentionRule, row: DataRow, missing: string[]): string {
  if (rule.message) {
    return rule.message.replace(/\{([^{}]+)\}/g, (_m, key: string) =>
      key === 'missing' ? missing.join(', ') || '—' : show(columnValue(row, key)));
  }
  const tested = (rule.condition.where ?? []).map((p) => `${p.column} ${show(columnValue(row, p.column))}`);
  return [rule.name, ...tested, ...(missing.length ? [`no ${missing.join(', ')}`] : [])].join('; ');
}

/** Which rows of the program dataset trip which rule. */
export function evaluateAttention(rules: AttentionRule[], rows: DataRow[], leadRoles = DEFAULT_LEAD_ROLES): ProgramAttention {
  const flags: AttentionFlag[] = [];
  for (const rule of rules) {
    const c = rule.condition;
    const severity = rule.severity ?? 'warning';
    const measured = c.measured ?? [];
    for (const {row, i} of rowsInScope((c.scope ?? 'program') as ConditionScope, rows, leadRoles)) {
      if (!(c.where ?? []).every((p) => holds(row, p)))
        continue;
      const missing = measured.filter((m) => isBlank(columnValue(row, m)));
      if (measured.length > 0 && missing.length === 0)
        continue;
      const values = Object.fromEntries([...(c.where ?? []).map((p) => p.column), ...measured]
        .map((col) => [col, columnValue(row, col) ?? null]));
      flags.push({rule: rule.code, severity, item: rowId(row, i), message: attentionMessage(rule, row, missing), values, missing});
    }
  }
  const byRule: RuleAttention[] = rules.map((r) => ({code: r.code, name: r.name, severity: r.severity ?? 'warning',
    items: [...new Set(flags.filter((f) => f.rule === r.code).map((f) => f.item))]}));
  const itemMap = new Map<string, AttentionItem>();
  for (const f of flags) {
    const it = itemMap.get(f.item) ?? {item: f.item, rules: [], severity: f.severity, messages: []};
    if (!it.rules.includes(f.rule))
      it.rules.push(f.rule);
    it.messages.push(f.message);
    if (rank(f.severity) < rank(it.severity))
      it.severity = f.severity;
    itemMap.set(f.item, it);
  }
  const sevOf = new Map(byRule.map((r) => [r.code, r.severity]));
  const items = [...itemMap.values()].map((it) => ({...it, rules: it.rules.sort((a, b) => rank(sevOf.get(a)!) - rank(sevOf.get(b)!))}))
    .sort((a, b) => rank(a.severity) - rank(b.severity) || a.item.localeCompare(b.item));
  const worst = items.length ? items[0].severity : null;
  const tripped = byRule.filter((r) => r.items.length > 0);
  return {
    rules: byRule, items, flags, count: items.length, worst,
    summary: items.length === 0 ? 'Nothing needs attention' :
      `${items.length} ${items.length === 1 ? 'item needs' : 'items need'} attention: ` +
      tripped.map((r) => `${r.name} ${r.items.length}`).join(', '),
  };
}

/** Every problem with a portfolio's attention rules, as sentences. */
export function validateAttentionRules(rules: AttentionRule[] | undefined): string[] {
  const errors: string[] = [];
  const codes = new Set<string>();
  for (const r of rules ?? []) {
    const where = `attention rule ${r?.code ?? '?'}`;
    if (!r?.code || !CODE.test(r.code))
      errors.push(`${where}: code is required and may hold only letters, digits, '_', '.', '-'`);
    else if (codes.has(r.code))
      errors.push(`${where}: code is used twice`);
    codes.add(r?.code);
    if (!r?.name)
      errors.push(`${where}: name is required`);
    if (r?.severity !== undefined && !ATTENTION_SEVERITIES.includes(r.severity))
      errors.push(`${where}: severity "${r.severity}" is not one of ${ATTENTION_SEVERITIES.join(', ')}`);
    if (!r?.condition) {
      errors.push(`${where}: condition is required`);
      continue;
    }
    const c = r.condition;
    if (c.minCount !== undefined || c.minSeries !== undefined)
      errors.push(`${where}: minCount and minSeries count rows for a gate; a rule flags rows one by one`);
    if (!c.where?.length && !c.measured?.length)
      errors.push(`${where}: give it where or measured, else it flags every row`);
    // the gate's own checks, less the one about counts that a rule does not have
    errors.push(...validateCondition({...c, minCount: undefined, minSeries: undefined})
      .filter((e) => !e.startsWith('the condition holds for any program')).map((e) => `${where}: ${e}`));
  }
  return errors;
}
