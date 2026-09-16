/// `grok kg serve`: the graph browser over one generation (core/docs/knowledge-graph/vis/). A loopback HTTP
/// server hands the page its render blob and answers, from the read-only index, what the page asks per
/// click: a node with its edge groups, an edge, a Cypher statement, one of the bounded operations.
import * as fs from 'fs';
import * as http from 'http';
import * as path from 'path';
import * as zlib from 'zlib';
import {TypeSystem} from './types';
import {KuzuConnection, KuzuDatabase, run, quote} from './kuzu';
import {impact, testsFor, explain, find, resolveTarget, DEFAULT_LIMIT} from './ops';
import {caveatNotes} from './answer';
import {Manifest} from './generation';
import {visDir, BLOB, INDEX, SCHEMA} from './vis';
import {Question, ask} from './questions';

export const PAGE_DIR = path.join(__dirname, 'vis');
/** Rows a Cypher answer may carry into the page. */
export const QUERY_CAP = 5000;
const OPS: Record<string, string> = {impact: 'impact', 'tests-for': 'tests-for', explain: 'explain', find: 'find'};
const MIME: Record<string, string> = {'.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8', '.json': 'application/json', '.bin': 'application/octet-stream', '.svg': 'image/svg+xml'};

export interface ServeOptions {
  genDir: string;
  /** Absolute, forward-slashed: the page builds editor links from it. */
  repoRoot: string;
  manifest: Manifest;
  questions: Question[];
  system: TypeSystem;
  db: KuzuDatabase;
  conn: KuzuConnection;
  port: number;
  host?: string;
}

export interface Served {
  server: http.Server;
  url: string;
  close(): Promise<void>;
}

class HttpError extends Error {
  constructor(readonly status: number, message: string) {
    super(message);
  }
}

/** Listens on the loopback address; resolves once the port is bound (0 picks a free one). */
export function serve(options: ServeOptions): Promise<Served> {
  const host = options.host ?? '127.0.0.1';
  const queue = new Queue();
  const server = http.createServer((req, res) => {
    handle(options, queue, req).then((answer) => send(req, res, answer)).catch((e) => {
      const status = e instanceof HttpError ? e.status : 500;
      send(req, res, json({error: e.message ?? String(e)}, status));
    });
  });
  return new Promise((resolve, reject) => {
    server.once('error', reject);
    server.listen(options.port, host, () => {
      const address = server.address() as {port: number};
      resolve({server, url: `http://${host}:${address.port}/`, close: () => new Promise((done) => server.close(() => done()))});
    });
  });
}

interface Answer {
  status: number;
  type: string;
  body: Buffer;
}

function json(value: unknown, status = 200): Answer {
  return {status, type: MIME['.json'], body: Buffer.from(JSON.stringify(value), 'utf8')};
}

async function handle(o: ServeOptions, queue: Queue, req: http.IncomingMessage): Promise<Answer> {
  const url = new URL(req.url ?? '/', 'http://localhost');
  const route = url.pathname;
  if (!route.startsWith('/api/')) return page(route);
  if (route === '/api/manifest')
    return json({...o.manifest, notes: caveatNotes(o.manifest), gen: path.basename(o.genDir), repo_root: o.repoRoot.replace(/\\/g, '/')});
  if (route === '/api/schema' || route === '/api/index' || route === '/api/graph.bin')
    return file(path.join(visDir(o.genDir), route === '/api/schema' ? SCHEMA : route === '/api/index' ? INDEX : BLOB));
  if (route === '/api/node') return queue.add(() => node(o, required(url, 'id')));
  if (route === '/api/edge') return queue.add(() => edge(o.conn, required(url, 'from'), required(url, 'to'), required(url, 'kind')));
  if (route === '/api/query') {
    if (req.method !== 'POST') throw new HttpError(405, 'POST a JSON body {cypher, params, limit}');
    const body = JSON.parse((await read(req)) || '{}');
    if (typeof body.cypher !== 'string' || !body.cypher.trim()) throw new HttpError(400, 'cypher is required');
    return queue.add(() => query(o.conn, body.cypher, limitOf(body.limit), body.params));
  }
  if (route === '/api/questions') return json(o.questions.map(({file, ...q}) => q));
  if (route.startsWith('/api/ask/')) {
    const id = route.slice('/api/ask/'.length);
    const question = o.questions.find((q) => q.id === id);
    if (!question) throw new HttpError(404, `no question ${id}`);
    const given = Object.fromEntries([...url.searchParams].filter(([k]) => k !== 'limit'));
    return queue.add(async () => {
      const answer = await ask(o.conn, question, given, o.manifest).catch((e) => { throw new HttpError(400, e.message); });
      const {question: q, ...rest} = answer;
      return json({...rest, id: q.id, highlight: q.highlight});
    });
  }
  const op = route.startsWith('/api/op/') ? OPS[route.slice('/api/op/'.length)] : undefined;
  if (op) return queue.add(() => operation(o, op, required(url, 'target'), limitOf(url.searchParams.get('limit'))));
  throw new HttpError(404, `no route ${route}`);
}

function required(url: URL, name: string): string {
  const value = url.searchParams.get(name);
  if (!value) throw new HttpError(400, `${name} is required`);
  return value;
}

function limitOf(value: unknown): number {
  const n = Number(value);
  return Number.isFinite(n) && n > 0 ? Math.min(Math.round(n), QUERY_CAP) : DEFAULT_LIMIT;
}

/** The page and its files: names only, from the one folder, so no path reaches the disk. */
function page(route: string): Answer {
  const name = route === '/' ? 'index.html' : path.basename(route);
  const file = path.join(PAGE_DIR, name);
  if (!/^[\w.-]+$/.test(name) || !fs.existsSync(file)) throw new HttpError(404, `no file ${name}`);
  return {status: 200, type: MIME[path.extname(name)] ?? 'application/octet-stream', body: fs.readFileSync(file)};
}

function file(p: string): Answer {
  if (!fs.existsSync(p)) throw new HttpError(404, `${path.basename(p)} is not built; restart grok kg serve`);
  return {status: 200, type: MIME[path.extname(p)] ?? 'application/octet-stream', body: fs.readFileSync(p)};
}

/** The node itself and its edge groups, the `edges` section `explain` answers with. */
async function node(o: ServeOptions, id: string): Promise<Answer> {
  const conn = o.conn;
  const target = await resolveTarget(conn, id);
  if (!target) throw new HttpError(404, `no node ${id}`);
  const found = await run(conn, `MATCH (n) WHERE n.${quote('id')} = $id RETURN n`, {id});
  const edges = (await explain(conn, target, {limit: QUERY_CAP, groups: o.manifest.edge_groups})).sections.find((s) => s.title === 'edges')!.rows;
  return json({node: found.rows[0].n, edges});
}

async function edge(conn: KuzuConnection, from: string, to: string, kind: string): Promise<Answer> {
  const {rows} = await run(conn, `MATCH (a)-[e]->(b) WHERE a.${quote('id')} = $from AND b.${quote('id')} = $to AND label(e) = $kind RETURN e`,
    {from, to, kind});
  if (!rows.length) throw new HttpError(404, `no ${kind} edge from ${from} to ${to}`);
  return json({edge: rows[0].e, from, to, kind});
}

/** A page never receives more than [limit] rows, whatever LIMIT the statement carries; one row more is asked for, so
 * `truncated` is a fact and not a guess. */
async function query(conn: KuzuConnection, cypher: string, limit: number, params?: Record<string, unknown>): Promise<Answer> {
  const statement = cypher.trim().replace(/;\s*$/, '');
  const bounded = /\bLIMIT\s+\d+\s*$/i.test(statement) ? statement : `${statement} LIMIT ${limit + 1}`;
  const started = Date.now();
  const bound = params && typeof params === 'object' && Object.keys(params).length ? params : undefined;
  const {columns, rows} = await run(conn, bounded, bound);
  return json({columns, rows: rows.slice(0, limit), ms: Date.now() - started, truncated: rows.length > limit});
}

async function operation(o: ServeOptions, op: string, text: string, limit: number): Promise<Answer> {
  const options = {limit, groups: o.manifest.edge_groups};
  let result;
  if (op === 'find') result = await find(o.conn, text, options);
  else {
    const target = await resolveTarget(o.conn, text);
    if (!target) throw new HttpError(404, `${text}: no such node`);
    result = op === 'impact' ? await impact(o.conn, target, options) : op === 'tests-for' ? await testsFor(o.conn, target, options) : await explain(o.conn, target, options);
  }
  result.notes = caveatNotes(o.manifest);
  return json(result);
}

function read(req: http.IncomingMessage): Promise<string> {
  return new Promise((resolve, reject) => {
    const chunks: Buffer[] = [];
    req.on('data', (c: Buffer) => chunks.push(c));
    req.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')));
    req.on('error', reject);
  });
}

function send(req: http.IncomingMessage, res: http.ServerResponse, answer: Answer): void {
  const headers: Record<string, string | number> = {'Content-Type': answer.type, 'Cache-Control': 'no-store'};
  let body = answer.body;
  if (body.length > 1024 && /\bgzip\b/.test(String(req.headers['accept-encoding'] ?? ''))) {
    body = zlib.gzipSync(body);
    headers['Content-Encoding'] = 'gzip';
  }
  headers['Content-Length'] = body.length;
  res.writeHead(answer.status, headers);
  res.end(body);
}

/** One statement at a time over the single connection, in arrival order. */
class Queue {
  private tail: Promise<unknown> = Promise.resolve();

  add<T>(work: () => Promise<T>): Promise<T> {
    const next = this.tail.then(work, work);
    this.tail = next.catch(() => undefined);
    return next;
  }
}
