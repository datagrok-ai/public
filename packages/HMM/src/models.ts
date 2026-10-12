/** Profile HMM sources: the built-in Pfam library shipped with the package,
 * Pfam families fetched from InterPro by accession, and user HMM files
 * (HMMER2/3 text or hmmpress `.h3m`, optionally gzipped). Each source is
 * loaded once into the worker pool under a content key. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type {DatabaseInfo} from './hmmer/engine.ts';
import {HmmerPool} from './pool';
import {cacheGet, cachePut} from './cache';

/** The built-in library (see files in wasm/pfam/: built by Rusty-HMMER web/pfam). */
export const BUILT_IN = 'pfam-biologics';

export interface LoadedModels {
  key: string;
  label: string;
  info: DatabaseInfo;
  /** Pfam clan accession of every model (null: no clan), or null when unknown
   * (user files); used for clan competition (hmmer/clans.ts). */
  clans: (string | null)[] | null;
  /** The model bytes as loaded (HMMER text or `.h3m`, gunzipped). */
  bytes: Uint8Array;
}

/** Bytes, gunzipped when they are a gzip stream. */
export async function gunzipIfNeeded(bytes: Uint8Array): Promise<Uint8Array> {
  if (bytes.length < 2 || bytes[0] !== 0x1f || bytes[1] !== 0x8b) return bytes;
  const stream = new Blob([bytes.slice()]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function contentKey(prefix: string, bytes: Uint8Array): Promise<string> {
  if (globalThis.crypto?.subtle) {
    const digest = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes.slice()));
    return prefix + ':' + Array.from(digest.slice(0, 12), (b) => b.toString(16).padStart(2, '0')).join('');
  }
  // Insecure (plain HTTP) origins have no SubtleCrypto: FNV-1a 32 plus the length.
  let hash = 0x811c9dc5;
  for (const b of bytes) hash = Math.imul(hash ^ b, 0x01000193) >>> 0;
  return `${prefix}:${hash.toString(16)}-${bytes.length}`;
}

/** Load model bytes into the pool (once per content). */
export async function loadModels(pool: HmmerPool, label: string, raw: Uint8Array,
  clans: (string | null)[] | null = null): Promise<LoadedModels> {
  const bytes = await gunzipIfNeeded(raw);
  const key = await contentKey('models', bytes);
  const info = await pool.loadDatabase(key, bytes);
  if (clans && clans.length !== info.models.length) clans = null;
  return {key, label, info, clans, bytes};
}

const builtIn = new Map<string, Promise<LoadedModels>>();

/** The built-in Pfam library (loaded on first use). */
export function loadBuiltIn(pool: HmmerPool, pkg: DG.Package): Promise<LoadedModels> {
  let loading = builtIn.get(pkg.webRoot);
  if (!loading) {
    loading = (async () => {
      const [models, index] = await Promise.all([`${BUILT_IN}.h3m.gz`, `${BUILT_IN}.json`].map(async (name) => {
        const response = await fetch(`${pkg.webRoot}dist/${name}`);
        if (!response.ok) throw new Error(`Failed to load the built-in Pfam library (HTTP ${response.status})`);
        return response;
      }));
      const manifest: {models: {clan: string | null}[]} = await index.json();
      return loadModels(pool, 'Pfam 37.0 biologics library', new Uint8Array(await models.arrayBuffer()),
        manifest.models.map((m) => m.clan ?? null));
    })();
    loading.catch(() => builtIn.delete(pkg.webRoot));
    builtIn.set(pkg.webRoot, loading);
  }
  return loading;
}

const germlines = new Map<string, Promise<LoadedModels>>();

/** ANARCI's germline models (all species and chains, IMGT 2026 build). */
export function loadGermlines(pool: HmmerPool, pkg: DG.Package): Promise<LoadedModels> {
  let loading = germlines.get(pkg.webRoot);
  if (!loading) {
    loading = (async () => {
      const response = await fetch(`${pkg.webRoot}dist/ALL.hmm.h3m.gz`);
      if (!response.ok) throw new Error(`Failed to load the ANARCI germline models (HTTP ${response.status})`);
      return loadModels(pool, 'ANARCI germline models', new Uint8Array(await response.arrayBuffer()));
    })();
    loading.catch(() => germlines.delete(pkg.webRoot));
    germlines.set(pkg.webRoot, loading);
  }
  return loading;
}

/** Pfam accessions (PF00001) in free text. */
export function parseAccessions(text: string): string[] {
  return [...new Set((text.toUpperCase().match(/PF\d{5}/g) ?? []))];
}

/** A Pfam family downloaded from InterPro, as cached in this browser. */
interface PfamEntry {
  /** HMMER3 text of the family's HMM. */
  hmm: Uint8Array;
  clan: string | null;
  fetched: number;
}

/** Cached families are refreshed after this long (InterPro follows Pfam releases). */
const CACHE_DAYS = 30;
const INTERPRO = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/';

async function download(url: string, accession: string): Promise<Response> {
  let response: Response;
  try {
    response = await grok.dapi.fetchProxy(url);
  } catch (e) {
    throw new Error(`Could not reach InterPro for Pfam ${accession}: ${e instanceof Error ? e.message : e}`);
  }
  if (response.status === 404 || response.status === 204)
    throw new Error(`Pfam ${accession} was not found in InterPro (mistyped or obsolete accession?)`);
  if (!response.ok) throw new Error(`InterPro returned HTTP ${response.status} for Pfam ${accession}`);
  return response;
}

/** One family's HMM and clan from InterPro (the current Pfam release). */
async function fetchEntry(accession: string): Promise<PfamEntry> {
  const [hmmResponse, metadata] = await Promise.all([
    download(`${INTERPRO}${accession}?annotation=hmm`, accession),
    download(`${INTERPRO}${accession}`, accession).then((r) => r.json()).catch(() => null),
  ]);
  const hmm = await gunzipIfNeeded(new Uint8Array(await hmmResponse.arrayBuffer()));
  if (!new TextDecoder().decode(hmm.subarray(0, 6)).startsWith('HMMER'))
    throw new Error(`InterPro did not return an HMM for Pfam ${accession}`);
  return {hmm, clan: metadata?.metadata?.set_info?.accession ?? null, fetched: metadata ? Date.now() : 0};
}

/** A family from the browser cache, refreshed from InterPro when stale or absent
 * (a stale copy is used when InterPro cannot be reached). */
async function pfamEntry(accession: string): Promise<PfamEntry> {
  const cached = await cacheGet<PfamEntry>(accession);
  if (cached && Date.now() - cached.fetched < CACHE_DAYS * 86_400_000) return cached;
  try {
    const entry = await fetchEntry(accession);
    // Without metadata (fetched = 0) the clan is unknown: keep it out of the cache.
    if (entry.fetched > 0) await cachePut(accession, entry);
    return entry;
  } catch (e) {
    if (cached) return cached;
    throw e;
  }
}

/** Pfam family HMMs from InterPro by accession (cached in the browser), with their clans. */
export async function fetchPfam(pool: HmmerPool, accessions: string[]): Promise<LoadedModels> {
  if (accessions.length === 0) throw new Error('No Pfam accessions (PF00000) given');
  const entries = await Promise.all(accessions.map(pfamEntry));
  const size = entries.reduce((n, e) => n + e.hmm.length, 0);
  const bytes = new Uint8Array(size);
  let at = 0;
  for (const {hmm} of entries) {
    bytes.set(hmm, at);
    at += hmm.length;
  }
  return loadModels(pool, `Pfam ${accessions.join(', ')}`, bytes, entries.map((e) => e.clan));
}

/** Load a user HMM file from a Datagrok file share. */
export async function loadFile(pool: HmmerPool, file: DG.FileInfo): Promise<LoadedModels> {
  return loadModels(pool, file.name, await file.readAsBytes());
}

/** Load HMM text given directly (e.g. by a script). */
export async function loadText(pool: HmmerPool, text: string): Promise<LoadedModels> {
  return loadModels(pool, 'HMM', new TextEncoder().encode(text));
}

/** Display label of a model: name, accession and description. */
export function modelLabel(info: DatabaseInfo, index: number): string {
  const m = info.models[index];
  const parts = [m.name, m.accession ? `(${m.accession})` : '', m.description ? `— ${m.description}` : ''];
  return parts.filter((s) => s).join(' ');
}
