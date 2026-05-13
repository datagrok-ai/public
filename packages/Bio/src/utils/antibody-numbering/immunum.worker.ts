// @ts-ignore -- file-loader emits the .wasm file to the package's dist/ dir.
// The returned string is the path relative to the emitted bundle; we take only
// the basename because the worker already lives in the same dist/ directory.
import wasmEmittedPath from 'immunum/immunum_bg.wasm';
// @ts-ignore -- plain JS port of immunum.js
import {Annotator, initImmunum} from './immunum-glue.js';

import type {ImmunumNumberingRow, ImmunumWorkerRequest, ImmunumWorkerResponse} from './types';

const ctx: Worker = self as any;

const DEFAULT_CHAINS = ['H', 'K', 'L'];
let _annotatorCache: {key: string; annotator: any} | null = null;
let _initPromise: Promise<void> | null = null;

function resolveWasmUrl(): string {
  // Worker script URL (e.g. http://host/packages/Bio/dist/NNN.js) is the right
  // base — the wasm sits in the same dist/ folder. Stripping off any path
  // prefix from file-loader's output (which may include a publicPath like
  // 'dist/') keeps us at the basename regardless of how webpack resolved it.
  const path = wasmEmittedPath as unknown as string;
  const basename = path.substring(path.lastIndexOf('/') + 1);
  return new URL(basename, self.location.href).toString();
}

async function ensureInit(): Promise<void> {
  if (!_initPromise) {
    _initPromise = (async () => {
      const wasmUrl = resolveWasmUrl();
      const resp = await fetch(wasmUrl);
      if (!resp.ok) throw new Error(`Failed to fetch immunum WASM from ${wasmUrl} (HTTP ${resp.status})`);
      const bytes = await resp.arrayBuffer();
      await initImmunum(bytes);
    })();
  }
  return _initPromise;
}

function getAnnotator(chains: string[], scheme: string, minConfidence: number | null): any {
  const key = `${scheme}|${chains.join(',')}|${minConfidence ?? 'null'}`;
  if (_annotatorCache && _annotatorCache.key === key) return _annotatorCache.annotator;
  if (_annotatorCache) {
    try { _annotatorCache.annotator.free(); } catch { /* ignore */ }
  }
  const annotator = new Annotator(chains, scheme, minConfidence ?? undefined);
  _annotatorCache = {key, annotator};
  return annotator;
}

const VALID_AA = new Set('ACDEFGHIKLMNPQRSTVWY');

/** Normalizes a raw input cell to a single-letter AA sequence: strips gaps, whitespace,
 *  non-standard characters, and uppercases. Matches the Python `extract_sequence` helper
 *  used by the antpack reference script so downstream char-index maps line up. */
function extractSequence(raw: string | null | undefined): string {
  if (!raw || typeof raw !== 'string') return '';
  let out = '';
  const s = raw.trim();
  for (let i = 0; i < s.length; i++) {
    const ch = s[i].toUpperCase();
    if (VALID_AA.has(ch)) out += ch;
  }
  return out;
}

function chainToGroup(chain: string | null): string {
  if (!chain) return '';
  if (chain === 'H') return 'Heavy';
  if (chain === 'K' || chain === 'L') return 'Light';
  return '';
}

function numberOne(annotator: any, rawSeq: string): ImmunumNumberingRow {
  const seq = extractSequence(rawSeq);
  if (!seq || seq.length < 10) {
    return {
      positionNames: '',
      chainType: '',
      chainCode: '',
      numberingDetail: [],
      numberingMap: {},
      confidence: 0,
      error: seq ? 'sequence too short' : 'empty sequence',
    };
  }

  let result;
  try {
    result = annotator.number(seq);
  } catch (e: any) {
    return {
      positionNames: '', chainType: '', chainCode: '',
      numberingDetail: [], numberingMap: {}, confidence: 0,
      error: (e && e.message) || String(e),
    };
  }

  if (!result || result.error || !result.numbering) {
    return {
      positionNames: '', chainType: '', chainCode: '',
      numberingDetail: [], numberingMap: {}, confidence: result?.confidence ?? 0,
      error: result?.error ?? 'no alignment',
    };
  }

  // numbering is a Map<posCode, residue>; iterate in insertion (IMGT) order.
  // Each entry corresponds to a single residue at consecutive input indices
  // starting at query_start (immunum semantics).
  const posNames: string[] = [];
  const numberingDetail: {position: string; aa: string}[] = [];
  const numberingMap: Record<string, number> = {};

  let charIdx = (result.query_start ?? 0) | 0;
  const numbering: Map<string, string> = result.numbering;
  for (const [posCode, residue] of numbering.entries()) {
    posNames.push(posCode);
    numberingDetail.push({position: posCode, aa: residue});
    numberingMap[posCode] = charIdx;
    charIdx++;
  }

  return {
    positionNames: posNames.join(', '),
    chainType: chainToGroup(result.chain),
    chainCode: result.chain ?? '',
    numberingDetail,
    numberingMap,
    confidence: result.confidence ?? 0,
    error: '',
  };
}

ctx.addEventListener('message', async (e: MessageEvent) => {
  const req: ImmunumWorkerRequest = e.data.req;
  const port: MessagePort = e.ports[0];

  try {
    await ensureInit();

    if (req.op === 'init') {
      port.postMessage({ok: true} satisfies ImmunumWorkerResponse);
      return;
    }

    if (req.op === 'number') {
      const chains = req.chains && req.chains.length > 0 ? req.chains : DEFAULT_CHAINS;
      const annotator = getAnnotator(chains, req.scheme, req.minConfidence ?? null);
      const rows: ImmunumNumberingRow[] = new Array(req.sequences.length);
      for (let i = 0; i < req.sequences.length; i++)
        rows[i] = numberOne(annotator, req.sequences[i]);
      port.postMessage({ok: true, rows} satisfies ImmunumWorkerResponse);
      return;
    }

    port.postMessage({ok: false, error: `Unknown op: ${(req as any).op}`} satisfies ImmunumWorkerResponse);
  } catch (err: any) {
    port.postMessage({ok: false, error: (err && err.message) || String(err)} satisfies ImmunumWorkerResponse);
  }
});
