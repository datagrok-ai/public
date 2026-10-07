/** HMMER worker: one engine instance (the compiled module is shared by the
 * pool), the ANARCI germline database and data loaded on first use, and
 * user model databases loaded by key. Requests are answered in order. */
import {HmmDatabase, HmmerEngine, type QueryResult} from '../hmmer/engine.ts';
import {align, numberFromAlignment, type AnarciData, type SequenceResult} from '../hmmer/anarci/anarci.ts';
import type {ChainGermlines} from '../hmmer/anarci/germline.ts';
import type {WorkerRequest, WorkerResponse} from './protocol';

const ctx = self as unknown as Worker;
let engine: Promise<HmmerEngine> | null = null;
let base = '';
let anarciData: Promise<AnarciData> | null = null;
const germlines = new Map<string, ChainGermlines | null>();
const databases = new Map<string, HmmDatabase>();

/** Response bytes, gunzipped when they are a gzip stream (servers may or may not decode .gz). */
async function bytesOf(response: Response): Promise<Uint8Array> {
  const bytes = new Uint8Array(await response.arrayBuffer());
  if (bytes.length < 2 || bytes[0] !== 0x1f || bytes[1] !== 0x8b) return bytes;
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function fetchFile(name: string): Promise<Response> {
  const response = await fetch(base + name);
  if (!response.ok) throw new Error(`Failed to load ${name} (HTTP ${response.status})`);
  return response;
}

async function loadAnarci(): Promise<AnarciData> {
  const [hmmer, db, lengths, species] = await Promise.all([
    engine!, fetchFile('ALL.hmm.h3m.gz').then(bytesOf),
    fetchFile('hmm-lengths.json').then((r) => r.json()), fetchFile('species.json').then((r) => r.json()),
  ]);
  return {
    database: hmmer.loadDatabase(db), hmmLengths: lengths, allSpecies: species,
    germlines: (chain) => germlines.get(chain) ?? null,
  };
}

async function anarci(request: Extract<WorkerRequest, {op: 'anarci'}>): Promise<SequenceResult[]> {
  const data = await (anarciData ??= loadAnarci());
  const alignments = align(data, request.sequences, request.options);
  if (request.options.assignGermline) {
    const chains = new Set(alignments.flatMap((a) => a.details.map((d) => d.chain_type)));
    await Promise.all([...chains].filter((c) => !germlines.has(c)).map(async (chain) => {
      try {
        germlines.set(chain, await (await fetchFile(`germlines-${chain}.json`)).json());
      } catch {
        germlines.set(chain, null);
      }
    }));
  }
  return request.sequences.map(([name, sequence], i) =>
    numberFromAlignment(data, name, sequence, alignments[i], request.options));
}

function database(key: string): HmmDatabase {
  const db = databases.get(key);
  if (!db) throw new Error(`Model database ${key} is not loaded in this worker`);
  return db;
}

async function handle(request: WorkerRequest): Promise<unknown> {
  switch (request.op) {
  case 'init':
    base = request.base;
    engine = HmmerEngine.instantiate(request.module);
    await engine;
    return null;
  case 'anarci':
    return anarci(request);
  case 'load': {
    const hmmer = await engine!;
    databases.get(request.key)?.free();
    const db = hmmer.loadDatabase(new Uint8Array(request.bytes));
    databases.set(request.key, db);
    return db.info;
  }
  case 'unload':
    databases.get(request.key)?.free();
    databases.delete(request.key);
    return null;
  case 'scan':
    return database(request.key).scan(request.sequences, request.options) satisfies QueryResult[];
  case 'search':
    // One part of a split search: a list, like every other batched request.
    return [{offset: request.offset,
      result: database(request.key).search(request.model, request.sequences, request.options) satisfies QueryResult}];
  }
}

let queue: Promise<void> = Promise.resolve();
ctx.addEventListener('message', (event: MessageEvent<{id: number; request: WorkerRequest}>) => {
  const {id, request} = event.data;
  queue = queue.then(async () => {
    let response: WorkerResponse;
    try {
      response = {id, ok: true, result: await handle(request)};
    } catch (e) {
      response = {id, ok: false, error: e instanceof Error ? e.message : String(e)};
    }
    ctx.postMessage(response);
  });
});
