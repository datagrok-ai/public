/** Profile tooling (F-18): portfolio files compared and bundled. A loaded portfolio exported back to
 * its file format, two portfolios diffed in readable lines, and a bundle that moves a set of
 * portfolios between instances. Pure: the platform side reads the rows and loads the set. */
import {PortfolioConfig, validatePortfolioSet} from './config';

export const PROFILE_FORMAT = 'project-portal-profile';
export const PROFILE_VERSION = 1;

/** A set of portfolio files, as one instance exports them and another imports them. */
export interface ProfileBundle {
  format: typeof PROFILE_FORMAT;
  version: number;
  exportedOn: string;
  portfolios: PortfolioConfig[];
}

export function makeBundle(portfolios: PortfolioConfig[], exportedOn: string): ProfileBundle {
  return {format: PROFILE_FORMAT, version: PROFILE_VERSION, exportedOn, portfolios};
}

/** The bundle in [text] (or an object already parsed), or the reasons it is not one. */
export function readBundle(input: string | unknown): {bundle: ProfileBundle | null, errors: string[]} {
  let b: any;
  try {
    b = typeof input === 'string' ? JSON.parse(input) : input;
  } catch (e: any) {
    return {bundle: null, errors: [`not JSON: ${e?.message ?? e}`]};
  }
  const errors: string[] = [];
  if (!b || typeof b !== 'object' || Array.isArray(b))
    errors.push('a profile bundle is an object with format, version and portfolios');
  else {
    if (b.format !== PROFILE_FORMAT)
      errors.push(`format is "${b.format}", not "${PROFILE_FORMAT}"`);
    if (typeof b.version !== 'number' || b.version > PROFILE_VERSION)
      errors.push(`version ${b.version} is not one this package reads (up to ${PROFILE_VERSION})`);
    if (!Array.isArray(b.portfolios) || b.portfolios.length === 0)
      errors.push('the bundle holds no portfolios');
  }
  return {bundle: errors.length ? null : b as ProfileBundle, errors};
}

// ---------------------------------------------------------------------------------------------
// Normal form: what loading makes of a file

const empty = (v: unknown) => v === undefined || v === null || v === '' || (Array.isArray(v) && v.length === 0);

/** Drops empty values and sorts object keys, recursively: a loaded portfolio does not keep an empty
 * list apart from none, and a json column comes back with its keys in the server's order. */
function clean(v: any): any {
  if (Array.isArray(v))
    return v.map(clean);
  if (v && typeof v === 'object') {
    const out: any = {};
    for (const [k, x] of Object.entries(v).sort(([a], [b]) => a.localeCompare(b))) {
      const c = clean(x);
      if (!empty(c) && !(c && typeof c === 'object' && !Array.isArray(c) && Object.keys(c).length === 0))
        out[k] = c;
    }
    return out;
  }
  return v;
}

/** The portfolio with every default loading applies written out and empty values dropped, so a
 * file and the same portfolio exported from an instance compare equal. */
export function normalizePortfolio(c: PortfolioConfig): PortfolioConfig {
  const n: any = JSON.parse(JSON.stringify(c));
  n.programNoun ??= 'Program';
  for (const s of n.stageModel?.stages ?? []) {
    s.isGate ??= true;
    s.isTerminal ??= false;
    for (const cr of s.criteria ?? [])
      cr.required ??= true;
  }
  for (const f of n.workUnitFeeds ?? [])
    f.system ??= `${n.code}-${f.kind}`;
  for (const r of n.catalogue?.runTypes ?? [])
    r.isExecutionMode ??= true;
  for (const m of n.catalogue?.models ?? []) {
    m.package ??= String(m.code ?? '').split(':')[0];
    m.status ??= 'active';
    if (m.card) {
      m.card.approval ??= 'draft';
      if (m.card.approvedOn)
        m.card.approvedOn = new Date(m.card.approvedOn).toISOString().slice(0, 10);
    }
  }
  return clean(n);
}

// ---------------------------------------------------------------------------------------------
// Diff

export interface PortfolioDiff {
  /** Stages, criteria, work units, series maturity, catalogue codes. */
  vocabulary: string[];
  /** Workspace root, workstreams, stage folders. */
  workspace: string[];
  /** Everything else: feeds, groups, datasets, rules, schedule, names. */
  fields: string[];
}

const show = (v: unknown): string => {
  if (v === undefined)
    return '—';
  const s = typeof v === 'string' ? `"${v}"` : JSON.stringify(v);
  return s.length > 80 ? `${s.slice(0, 77)}…` : s;
};

/** One line for a value that changed; long text is named, not printed. */
function changed(path: string, a: unknown, b: unknown): string {
  if (a === undefined)
    return `${path}: set to ${show(b)}`;
  if (b === undefined)
    return `${path}: removed (was ${show(a)})`;
  if (typeof a === 'string' && typeof b === 'string' && (a.length > 60 || b.length > 60 || a.includes('\n') || b.includes('\n')))
    return `${path}: changed (${a.length} → ${b.length} characters)`;
  return `${path}: ${show(a)} → ${show(b)}`;
}

const keyOf = (x: any): string | null =>
  x && typeof x === 'object' && !Array.isArray(x) ? (typeof x.code === 'string' ? x.code : typeof x.name === 'string' ? x.name : null) : null;

/** Lines for every difference between [a] and [b] under [path]: objects key by key, lists of
 * coded (or named) items item by item, anything else as a whole. */
export function deepDiff(a: unknown, b: unknown, path: string, out: string[] = []): string[] {
  if (JSON.stringify(a) === JSON.stringify(b))
    return out;
  if (Array.isArray(a) && Array.isArray(b) && [...a, ...b].every((x) => keyOf(x) !== null)) {
    const ka = new Map(a.map((x) => [keyOf(x)!, x]));
    const kb = new Map(b.map((x) => [keyOf(x)!, x]));
    for (const [k] of ka)
      if (!kb.has(k)) out.push(`${path} ${k}: removed`);
    for (const [k] of kb)
      if (!ka.has(k)) out.push(`${path} ${k}: added`);
    for (const [k, x] of ka)
      if (kb.has(k)) deepDiff(x, kb.get(k), `${path} ${k}`, out);
    return out;
  }
  const isObj = (x: unknown) => x !== null && typeof x === 'object' && !Array.isArray(x);
  if (isObj(a) && isObj(b)) {
    const keys = [...new Set([...Object.keys(a as any), ...Object.keys(b as any)])];
    for (const k of keys)
      deepDiff((a as any)[k], (b as any)[k], path ? `${path}.${k}` : k, out);
    return out;
  }
  if (Array.isArray(a) && Array.isArray(b) && [...a, ...b].every((x) => typeof x !== 'object')) {
    const added = b.filter((x) => !a.includes(x));
    const removed = a.filter((x) => !b.includes(x));
    if (added.length || removed.length)
      out.push(`${path}: ${[added.length ? `added ${added.join(', ')}` : '', removed.length ? `removed ${removed.join(', ')}` : '']
        .filter((s) => s).join('; ')}`);
    else
      out.push(`${path}: reordered ${a.join(', ')} → ${b.join(', ')}`);
    return out;
  }
  out.push(changed(path, a, b));
  return out;
}

const omit = (o: any, keys: string[]) => Object.fromEntries(Object.entries(o ?? {}).filter(([k]) => !keys.includes(k)));

/** What changes from [a] to [b], in readable lines by kind. Both are normalized first, so a default
 * written out on one side and left out on the other is no difference. */
export function diffPortfolios(a: PortfolioConfig, b: PortfolioConfig): PortfolioDiff {
  const x: any = normalizePortfolio(a);
  const y: any = normalizePortfolio(b);
  const vocabulary: string[] = [];
  const workspace: string[] = [];
  const fields: string[] = [];

  // the stage model: stages and criteria are vocabulary, stage folders are the workspace tree
  const ma = x.stageModel ?? {};
  const mb = y.stageModel ?? {};
  deepDiff(omit(ma, ['stages']), omit(mb, ['stages']), 'stage model', vocabulary);
  const sa = new Map<string, any>((ma.stages ?? []).map((s: any) => [s.code, s]));
  const sb = new Map<string, any>((mb.stages ?? []).map((s: any) => [s.code, s]));
  for (const [code, s] of sa)
    if (!sb.has(code)) vocabulary.push(`stage ${code} (${s.name}): removed`);
  for (const [code, s] of sb) {
    if (!sa.has(code)) {
      vocabulary.push(`stage ${code} (${s.name}): added at ordinal ${s.ordinal}` +
        `${s.criteria?.length ? ` with criteria ${s.criteria.map((c: any) => c.code).join(', ')}` : ''}`);
      if (s.folders?.length)
        workspace.push(`stage ${code} folders: added ${s.folders.join(', ')}`);
    }
  }
  for (const [code, s] of sa) {
    const t = sb.get(code);
    if (!t)
      continue;
    deepDiff(omit(s, ['criteria', 'folders']), omit(t, ['criteria', 'folders']), `stage ${code}`, vocabulary);
    deepDiff(s.criteria ?? [], t.criteria ?? [], `criterion ${code}/`, vocabulary);
    deepDiff(s.folders ?? [], t.folders ?? [], `stage ${code} folders`, workspace);
  }

  deepDiff(x.workUnits ?? [], y.workUnits ?? [], 'work units', vocabulary);
  deepDiff(x.seriesMaturity, y.seriesMaturity, 'series maturity', vocabulary);
  const ca = x.catalogue ?? {};
  const cb = y.catalogue ?? {};
  deepDiff(ca.unitOperations ?? [], cb.unitOperations ?? [], 'unit operation', vocabulary);
  deepDiff(ca.runTypes ?? [], cb.runTypes ?? [], 'run type', vocabulary);
  deepDiff(ca.models ?? [], cb.models ?? [], 'model', vocabulary);
  deepDiff(x.documentKinds ?? [], y.documentKinds ?? [], 'document kind', vocabulary);

  deepDiff(x.workspaceRoot, y.workspaceRoot, 'workspace root', workspace);
  deepDiff(x.workstreams ?? [], y.workstreams ?? [], 'workstreams', workspace);

  const handled = ['stageModel', 'workUnits', 'seriesMaturity', 'catalogue', 'documentKinds', 'workspaceRoot', 'workstreams'];
  deepDiff(omit(x, handled), omit(y, handled), '', fields);
  // `criterion H2L/ potency` reads better as `criterion H2L/potency`
  const tidy = (ls: string[]) => ls.map((l) => l.replace(/\/ /g, '/'));
  return {vocabulary: tidy(vocabulary), workspace: tidy(workspace), fields: tidy(fields)};
}

export const diffCount = (d: PortfolioDiff) => d.vocabulary.length + d.workspace.length + d.fields.length;

/** The diff as text: a heading per kind that has lines. */
export function summarizeDiff(title: string, d: PortfolioDiff): string {
  if (diffCount(d) === 0)
    return `${title}: no differences`;
  const lines = [`${title}: ${diffCount(d)} difference(s)`];
  for (const [heading, ls] of [['Vocabulary', d.vocabulary], ['Workspace tree', d.workspace], ['Other fields', d.fields]] as const) {
    if (ls.length)
      lines.push(`  ${heading}:`, ...ls.map((l) => `    ${l}`));
  }
  return lines.join('\n');
}

// ---------------------------------------------------------------------------------------------
// Import plan

export interface ImportPlan {
  /** Why the bundle cannot be loaded; empty when it can. */
  errors: string[];
  portfolios: {code: string, status: 'new' | 'changed' | 'unchanged', diff: PortfolioDiff | null}[];
}

/** What importing [bundle] would do on an instance that has [loaded] (its exported portfolios, by
 * code). The bundle is checked as a set; a handoff may name a portfolio the instance already has. */
export function planImport(bundle: ProfileBundle, loaded: Map<string, PortfolioConfig>, loadedCodes: string[]): ImportPlan {
  const errors = validatePortfolioSet(bundle.portfolios, loadedCodes);
  const portfolios = bundle.portfolios.map((p) => {
    const was = loaded.get(p.code);
    if (!was)
      return {code: p.code, status: 'new' as const, diff: null};
    const diff = diffPortfolios(was, p);
    return {code: p.code, status: diffCount(diff) ? 'changed' as const : 'unchanged' as const, diff};
  });
  return {errors, portfolios};
}

export function summarizeImport(plan: ImportPlan, dryRun: boolean): string {
  if (plan.errors.length)
    return [`Refused: the bundle does not validate (${plan.errors.length} problem(s))`, ...plan.errors.map((e) => `  ${e}`)].join('\n');
  const lines = [`${dryRun ? 'Would import' : 'Imported'} ${plan.portfolios.length} portfolio(s): ` +
    ['new', 'changed', 'unchanged'].map((s) => `${plan.portfolios.filter((p) => p.status === s).length} ${s}`).join(', ')];
  for (const p of plan.portfolios)
    lines.push(p.status === 'new' ? `${p.code}: new` : p.diff ? summarizeDiff(p.code, p.diff) : `${p.code}: ${p.status}`);
  if (plan.portfolios.some((p) => p.diff && diffCount(p.diff) > 0 && p.diff && [...p.diff.vocabulary, ...p.diff.workspace,
    ...p.diff.fields].some((l) => /: removed/.test(l))))
    lines.push('Loading adds and updates; it never removes what the bundle leaves out (stages, criteria, feeds, catalogue entries, folders).');
  return lines.join('\n');
}
