/// The question set (core/docs/knowledge-graph/questions/): one YAML file per question a developer or an agent
/// asks the graph, its Cypher, its parameters and what it needs. `grok kg ask` and the browser run them; the
/// benchmarks reuse them as the yardstick of what the graph can answer.
import * as fs from 'fs';
import * as path from 'path';
import * as yaml from 'js-yaml';
import {KuzuConnection, run, QueryRows} from './kuzu';
import {sourceCaveats} from './ops';

export const QUESTIONS_DIR = 'questions';
const TYPES = ['string', 'number', 'date'] as const;
const RELATIVE = /^([+-]\d+)d$/;

export interface Param {
  type: typeof TYPES[number];
  default?: string | number;
  description?: string;
}

export interface Question {
  id: string;
  file: string;
  question: string;
  why: string;
  tags: string[];
  params: Record<string, Param>;
  cypher: string;
  /** Columns that carry node ids, for the browser to light up. */
  highlight: string[];
  /** Manifest sources the answer depends on; a caveat is attached when one is not `ok`. */
  needs: string[];
  status: 'ok' | 'blocked';
  blocked_by?: string;
}

export interface Answer extends QueryRows {
  question: Question;
  params: Record<string, string | number>;
  ms: number;
  notes: string[];
}

/** Every `<kgRoot>/questions/*.yaml`, by id; a malformed file is an error the caller reports. */
export function loadQuestions(kgRoot: string): {questions: Question[], errors: string[]} {
  const dir = path.join(kgRoot, QUESTIONS_DIR);
  const questions: Question[] = [];
  const errors: string[] = [];
  if (!fs.existsSync(dir)) return {questions, errors: [`${dir}: no questions folder`]};
  for (const file of fs.readdirSync(dir).filter((f) => f.endsWith('.yaml')).sort()) {
    const full = path.join(dir, file);
    try {
      const doc = yaml.load(fs.readFileSync(full, 'utf8')) as Record<string, unknown>;
      questions.push(parse(doc, file, full));
    }
    catch (e: any) {
      errors.push(`${file}: ${e.message}`);
    }
  }
  return {questions, errors};
}

function parse(doc: Record<string, unknown>, file: string, full: string): Question {
  const id = file.replace(/\.yaml$/, '');
  if (doc.id !== id) throw new Error(`id '${doc.id}' does not match the file name`);
  for (const key of ['question', 'why', 'cypher']) if (typeof doc[key] !== 'string' || !String(doc[key]).trim()) throw new Error(`${key} is required`);
  const params: Record<string, Param> = {};
  for (const [name, spec] of Object.entries((doc.params ?? {}) as Record<string, Param>)) {
    if (!/^[a-z][a-z0-9_]*$/.test(name)) throw new Error(`param '${name}': not a name`);
    if (!TYPES.includes(spec?.type)) throw new Error(`param '${name}': type must be ${TYPES.join(', ')}`);
    params[name] = spec;
  }
  for (const name of Object.keys(params))
    if (!doc.cypher || !String(doc.cypher).includes(`$${name}`)) throw new Error(`param '${name}' is not used in the cypher`);
  const status = doc.status ?? 'ok';
  if (status !== 'ok' && status !== 'blocked') throw new Error(`status must be ok or blocked`);
  if (status === 'blocked' && typeof doc.blocked_by !== 'string') throw new Error('a blocked question says what blocks it in blocked_by');
  return {id, file: full, question: String(doc.question), why: String(doc.why), tags: list(doc.tags), params, cypher: String(doc.cypher),
    highlight: doc.highlight === undefined ? ['id'] : list(doc.highlight), needs: list(doc.needs), status, blocked_by: doc.blocked_by as string | undefined};
}

function list(value: unknown): string[] {
  return Array.isArray(value) ? value.map(String) : [];
}

/** The parameters a run uses: what the caller gave, typed, over the defaults; a relative date resolves now. */
export function resolveParams(question: Question, given: Record<string, unknown>, now = new Date()): Record<string, string | number> {
  const out: Record<string, string | number> = {};
  for (const [name, spec] of Object.entries(question.params)) {
    const raw = given[name] ?? spec.default;
    if (raw === undefined || raw === null || raw === '') throw new Error(`${question.id}: parameter '${name}' has no value and no default`);
    if (spec.type === 'number') {
      const n = Number(raw);
      if (!Number.isFinite(n)) throw new Error(`${question.id}: parameter '${name}' must be a number, got '${raw}'`);
      out[name] = n;
    }
    else if (spec.type === 'date') out[name] = isoDate(String(raw), now, `${question.id}: parameter '${name}'`);
    else out[name] = String(raw);
  }
  for (const name of Object.keys(given))
    if (!(name in question.params)) throw new Error(`${question.id} takes no parameter '${name}'`);
  return out;
}

/** `-7d` is seven days before [now] at midnight UTC; anything else must be an ISO date or timestamp. */
export function isoDate(value: string, now: Date, what: string): string {
  const relative = RELATIVE.exec(value.trim());
  if (relative) {
    const d = new Date(Date.UTC(now.getUTCFullYear(), now.getUTCMonth(), now.getUTCDate() + Number(relative[1])));
    return d.toISOString();
  }
  const parsed = new Date(value);
  if (Number.isNaN(parsed.getTime())) throw new Error(`${what} must be an ISO date or a relative day count like -7d, got '${value}'`);
  return parsed.toISOString();
}

export async function ask(conn: KuzuConnection, question: Question, given: Record<string, unknown>, sources?: Record<string, string>): Promise<Answer> {
  const params = resolveParams(question, given);
  const started = Date.now();
  const rows = await run(conn, question.cypher, Object.keys(params).length ? params : undefined);
  const notes = sourceCaveats(sources).filter((n) => question.needs.some((s) => n.startsWith(`${s} `) || (s === 'dart' && n.startsWith('Dart'))));
  if (question.status === 'blocked') notes.unshift(`blocked: ${question.blocked_by}`);
  return {...rows, question, params, ms: Date.now() - started, notes};
}
