/** Work units beneath a program (FR-06): campaign, series, study, process stage. A portfolio enables
 * the subset it uses and may feed them from a source system; units are keyed by (program, code).
 * Pure: the caller reads the feed rows and the units the portal holds, and applies the plan. Like
 * programs, a unit its feed stops delivering is flagged, never deleted. */

export type WorkUnitKind = 'campaign' | 'series' | 'study' | 'process_stage';

export const WORK_UNIT_KINDS: WorkUnitKind[] = ['campaign', 'series', 'study', 'process_stage'];

export const WORK_UNIT_LABELS: Record<WorkUnitKind, string> = {
  campaign: 'Campaigns', series: 'Series', study: 'Studies', process_stage: 'Process stages',
};

export interface WorkUnitFeedConfig {
  code: string;
  kind: WorkUnitKind;
  /** The source system; defaults to `<portfolio>-<kind>`. */
  system?: string;
  connection: string;
  /** One row per unit: `program` (the program's code), `code` and `name` are required; the kind's
   * own fields are mapped when present; anything else lands in `details`. */
  sql: string;
}

type FieldKind = 'string' | 'int' | 'datetime';
interface FieldSpec {
  type: FieldKind;
  choices?: string[];
}

/** The columns a feed may set, per kind; the schema's choices repeated so a plan can be checked
 * without the server. */
export const WORK_UNIT_FIELDS: Record<WorkUnitKind, Record<string, FieldSpec>> = {
  campaign: {
    type: {type: 'string', choices: ['screening', 'hit_triage', 'dmta_cycle', 'library', 'other']},
    status: {type: 'string', choices: ['planned', 'active', 'paused', 'closed']},
    outcome: {type: 'string', choices: ['open', 'succeeded', 'abandoned']},
    objective: {type: 'string'},
    started_on: {type: 'datetime'},
  },
  series: {
    maturity: {type: 'string'},
    status: {type: 'string', choices: ['active', 'backup', 'deprioritized', 'dropped']},
    scaffold: {type: 'string'},
  },
  study: {
    type: {type: 'string'},
    status: {type: 'string', choices: ['planned', 'running', 'completed', 'cancelled']},
    planned_end: {type: 'datetime'},
    completed_on: {type: 'datetime'},
    result_summary: {type: 'string'},
  },
  process_stage: {
    process_type: {type: 'string'},
    ordinal: {type: 'int'},
    status: {type: 'string', choices: ['not_started', 'in_development', 'locked', 'transferred']},
  },
};

const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

/** Problems with a portfolio's work-unit declarations, as sentences. */
export function validateWorkUnits(enabled: string[] | undefined, feeds: WorkUnitFeedConfig[] | undefined,
  seriesMaturity: string[] | undefined): string[] {
  const errors: string[] = [];
  const kinds = new Set<string>();
  for (const k of enabled ?? []) {
    if (!(WORK_UNIT_KINDS as string[]).includes(k))
      errors.push(`work unit "${k}" is not one of ${WORK_UNIT_KINDS.join(', ')}`);
    else if (kinds.has(k))
      errors.push(`work unit "${k}" is enabled twice`);
    kinds.add(k);
  }
  const codes = new Set<string>();
  for (const f of feeds ?? []) {
    if (!f.code || !CODE.test(f.code))
      errors.push(`work-unit feed code "${f.code ?? ''}" may hold only letters, digits, '_', '.', '-'`);
    else if (codes.has(f.code))
      errors.push(`work-unit feed code "${f.code}" is declared twice`);
    codes.add(f.code);
    if (!kinds.has(f.kind))
      errors.push(`work-unit feed ${f.code} feeds ${f.kind}, which the portfolio does not enable`);
    if (!f.connection)
      errors.push(`work-unit feed ${f.code}: connection is required`);
    if (!f.sql)
      errors.push(`work-unit feed ${f.code}: sql is required`);
  }
  if (seriesMaturity !== undefined) {
    if (!kinds.has('series'))
      errors.push('seriesMaturity needs the series work unit enabled');
    const seen = new Set<string>();
    for (const m of seriesMaturity) {
      if (!m || seen.has(m.toLowerCase()))
        errors.push(`series maturity "${m}" is blank or declared twice`);
      seen.add((m ?? '').toLowerCase());
    }
  }
  return errors;
}

export interface WorkUnitRef {
  id: string;
  programId: string;
  code: string;
  source?: string | null;
  missingFromFeed?: boolean | null;
  values: Record<string, unknown>;
}

export interface WorkUnitValues {
  program_id: string;
  code: string;
  name: string;
  details: Record<string, unknown> | null;
  [field: string]: unknown;
}

export interface WorkUnitPlan {
  portfolio: string;
  kind: WorkUnitKind;
  feed: string;
  inserts: WorkUnitValues[];
  updates: {id: string, key: string, changes: Record<string, unknown>}[];
  unchanged: string[];
  /** Units this feed used to deliver and no longer does: flagged, kept. */
  missing: {id: string, key: string}[];
  /** Flagged units the feed delivers again. */
  returned: {id: string, key: string}[];
  rejected: {row: number, reason: string}[];
  warnings: string[];
}

export interface WorkUnitPlanInput {
  portfolio: string;
  kind: WorkUnitKind;
  feed: string;
  rows: Record<string, unknown>[];
  /** The portfolio's programs by code. */
  programs: Map<string, string>;
  /** Units of this kind under the portfolio's programs. */
  existing: WorkUnitRef[];
  seriesMaturity?: string[];
}

const blank = (v: unknown) => v === null || v === undefined || (typeof v === 'string' && v.trim() === '');

/** A datetime as the server returns it (dayjs, Date or text) in one comparable form. */
export function normalizeDate(v: unknown): string | null {
  if (blank(v))
    return null;
  const iso = (v as any)?.toISOString?.() ?? new Date(String(v)).toISOString?.();
  return typeof iso === 'string' && !iso.startsWith('Invalid') ? iso : null;
}

function normalize(spec: FieldSpec, v: unknown): unknown {
  if (blank(v))
    return null;
  switch (spec.type) {
  case 'int': {
    const n = Number(v);
    return Number.isInteger(n) ? n : null;
  }
  case 'datetime': {
    const d = new Date(String((v as any)?.toISOString?.() ?? v));
    return Number.isNaN(d.getTime()) ? null : d.toISOString();
  }
  default: return String(v).trim();
  }
}

const same = (a: unknown, b: unknown) => JSON.stringify(a ?? null) === JSON.stringify(b ?? null);

export function planWorkUnits(input: WorkUnitPlanInput): WorkUnitPlan {
  const plan: WorkUnitPlan = {portfolio: input.portfolio, kind: input.kind, feed: input.feed, inserts: [], updates: [],
    unchanged: [], missing: [], returned: [], rejected: [], warnings: []};
  const fields = WORK_UNIT_FIELDS[input.kind];
  const maturity = new Set((input.seriesMaturity ?? []).map((m) => m.toLowerCase()));
  const existing = new Map(input.existing.map((u) => [`${u.programId}/${u.code}`, u]));
  const delivered = new Set<string>();
  const label = (program: string, code: string) => `${program}/${code}`;

  input.rows.forEach((row, i) => {
    const program = blank(row.program) ? '' : String(row.program).trim();
    const code = blank(row.code) ? '' : String(row.code).trim();
    if (!program || !code) {
      plan.rejected.push({row: i + 1, reason: 'program and code are required'});
      return;
    }
    const programId = input.programs.get(program);
    if (!programId) {
      plan.rejected.push({row: i + 1, reason: `${label(program, code)}: ${program} is not a program of ${input.portfolio}`});
      return;
    }
    const key = `${programId}/${code}`;
    if (delivered.has(key)) {
      plan.rejected.push({row: i + 1, reason: `${label(program, code)} is delivered twice`});
      return;
    }
    delivered.add(key);
    const values: WorkUnitValues = {program_id: programId, code, name: blank(row.name) ? code : String(row.name).trim(),
      details: null};
    const details: Record<string, unknown> = {};
    for (const [col, v] of Object.entries(row)) {
      if (['program', 'code', 'name'].includes(col))
        continue;
      const spec = fields[col];
      if (!spec) {
        if (!blank(v))
          details[col] = (v as any)?.toISOString?.() ?? v;
        continue;
      }
      const value = normalize(spec, v);
      if (value !== null && spec.choices && !spec.choices.includes(String(value))) {
        plan.warnings.push(`${label(program, code)}: ${col} "${value}" is not one of ${spec.choices.join(', ')}; left unset`);
        continue;
      }
      if (value !== null && col === 'maturity' && maturity.size > 0 && !maturity.has(String(value).toLowerCase()))
        plan.warnings.push(`${label(program, code)}: maturity "${value}" is not in the portfolio's series maturity`);
      values[col] = value;
    }
    values.details = Object.keys(details).length ? details : null;

    const was = existing.get(key);
    if (!was) {
      plan.inserts.push(values);
      return;
    }
    const changes: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(values)) {
      if (k === 'program_id' || k === 'code')
        continue;
      const before = fields[k]?.type === 'datetime' ? normalizeDate(was.values[k]) : was.values[k];
      if (!same(before, v))
        changes[k] = v;
    }
    if (was.missingFromFeed)
      plan.returned.push({id: was.id, key: label(program, code)});
    if (Object.keys(changes).length)
      plan.updates.push({id: was.id, key: label(program, code), changes});
    else if (!was.missingFromFeed)
      plan.unchanged.push(label(program, code));
  });

  const codeOf = new Map([...input.programs.entries()].map(([c, id]) => [id, c]));
  for (const u of input.existing) {
    if (u.source === 'feed' && !u.missingFromFeed && !delivered.has(`${u.programId}/${u.code}`))
      plan.missing.push({id: u.id, key: label(codeOf.get(u.programId) ?? '?', u.code)});
  }
  return plan;
}

export function summarizeWorkUnitPlan(p: WorkUnitPlan): string {
  const parts = [`${p.portfolio} ${p.kind} (${p.feed}): ${p.inserts.length} new, ${p.updates.length} updated, ` +
    `${p.unchanged.length} unchanged, ${p.missing.length} missing from feed`];
  if (p.returned.length)
    parts.push(`  back in the feed: ${p.returned.map((r) => r.key).join(', ')}`);
  for (const r of p.rejected)
    parts.push(`  rejected row ${r.row}: ${r.reason}`);
  for (const w of p.warnings)
    parts.push(`  ${w}`);
  return parts.join('\n');
}
