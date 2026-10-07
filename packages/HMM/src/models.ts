/** Profile HMM sources: the built-in Pfam library shipped with the package,
 * Pfam families fetched from InterPro by accession, and user HMM files
 * (HMMER2/3 text or hmmpress `.h3m`, optionally gzipped). Each source is
 * loaded once into the worker pool under a content key. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type {DatabaseInfo} from './hmmer/engine.ts';
import {HmmerPool} from './pool';

/** The built-in library (see files in wasm/pfam/: built by Rusty-HMMER web/pfam). */
export const BUILT_IN = 'pfam-biologics';

export interface LoadedModels {
  key: string;
  label: string;
  info: DatabaseInfo;
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
export async function loadModels(pool: HmmerPool, label: string, raw: Uint8Array): Promise<LoadedModels> {
  const bytes = await gunzipIfNeeded(raw);
  const key = await contentKey('models', bytes);
  const info = await pool.loadDatabase(key, bytes);
  return {key, label, info};
}

const builtIn = new Map<string, Promise<LoadedModels>>();

/** The built-in Pfam library (loaded on first use). */
export function loadBuiltIn(pool: HmmerPool, pkg: DG.Package): Promise<LoadedModels> {
  let loading = builtIn.get(pkg.webRoot);
  if (!loading) {
    loading = (async () => {
      const response = await fetch(`${pkg.webRoot}dist/${BUILT_IN}.h3m.gz`);
      if (!response.ok) throw new Error(`Failed to load the built-in Pfam library (HTTP ${response.status})`);
      return loadModels(pool, 'Pfam 37.0 biologics library', new Uint8Array(await response.arrayBuffer()));
    })();
    loading.catch(() => builtIn.delete(pkg.webRoot));
    builtIn.set(pkg.webRoot, loading);
  }
  return loading;
}

/** Pfam accessions (PF00001) in free text. */
export function parseAccessions(text: string): string[] {
  return [...new Set((text.toUpperCase().match(/PF\d{5}/g) ?? []))];
}

/** Fetch Pfam family HMMs from InterPro (current Pfam release) through the server proxy. */
export async function fetchPfam(pool: HmmerPool, accessions: string[]): Promise<LoadedModels> {
  if (accessions.length === 0) throw new Error('No Pfam accessions (PF00000) given');
  const parts = await Promise.all(accessions.map(async (accession) => {
    const url = `https://www.ebi.ac.uk/interpro/api/entry/pfam/${accession}?annotation=hmm`;
    const response = await grok.dapi.fetchProxy(url);
    if (!response.ok) throw new Error(`Pfam ${accession}: HTTP ${response.status}`);
    return gunzipIfNeeded(new Uint8Array(await response.arrayBuffer()));
  }));
  const size = parts.reduce((n, p) => n + p.length, 0);
  const bytes = new Uint8Array(size);
  let at = 0;
  for (const part of parts) {
    bytes.set(part, at);
    at += part.length;
  }
  return loadModels(pool, `Pfam ${accessions.join(', ')}`, bytes);
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
