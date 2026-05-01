/**
 * Hover tooltip for OligoNucleotide cells.
 *
 * Shows monomer details synchronously, then async-loads RDKit structures
 * for whatever was hovered:
 *   - Nucleotide chip → sugar + base + 3'-phosphate (3 small structures)
 *   - PS / phosphate-mod linkage in the gap → just the phosphate structure
 *   - Conjugate pill → that conjugate's structure
 *
 * Structures come from the **central Bio monomer library** (not the
 * SequenceTranslator-internal one). Symbols are resolved via the alias map
 * so legacy / vendor shorthand still finds the canonical entry.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {getMonomerLibHelper, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';

import {HitResult} from './canvas-renderer';
import {
  canonicalPhosphateSymbol, canonicalSugarSymbol,
  ParsedNucleotide, resolveConjugate, resolvePhosphate, resolveSugar,
} from './types';

const STRUCT_W = 110;
const STRUCT_H = 90;

type StructureKind = 'sugar' | 'base' | 'phosphate' | 'conjugate';

/** Show a tooltip for the hovered hit. Idempotent — replaces any current tooltip. */
export function showMonomerTooltip(hit: HitResult, clientX: number, clientY: number): void {
  const root = buildTooltipRoot(hit);
  ui.tooltip.show(root, clientX + 14, clientY + 14);
  void renderStructuresAsync(root, hit);
}

function buildTooltipRoot(hit: HitResult): HTMLElement {
  const m = hit.monomer;
  const title = ui.divText('', {style: {fontWeight: '600', marginBottom: '4px'}});
  const detailRows: HTMLElement[] = [];

  // Linkage hover (PS marker in the gap) — shape the tooltip around the linkage
  if (hit.linkage) {
    const phos = resolvePhosphate(hit.linkage.phosphateSymbol);
    title.textContent = `${phos.meta.name} linkage`;
    detailRows.push(
      kv('Strand', hit.strand === 'sense' ? 'Sense' : 'Antisense'),
      kv('Between', `pos ${hit.position + 1} – ${hit.position + 2}`),
      kv('HELM symbol', `[${canonicalPhosphateSymbol(hit.linkage.phosphateSymbol)}]`, true),
    );
  } else if (m.kind === 'nucleotide') {
    const nt = m as ParsedNucleotide;
    const sugarRes = resolveSugar(nt.sugar, nt.base);
    const phosRes = resolvePhosphate(nt.phosphate);
    title.textContent = `${nt.base ?? '?'}${hit.position + 1} — ${sugarRes.meta.name}`;
    detailRows.push(
      kv('Strand', hit.strand === 'sense' ? 'Sense' : 'Antisense'),
      kv('Position', String(hit.position + 1)),
      kv('Base', nt.base ?? '—'),
      kv('Sugar', sugarRes.meta.name),
      kv('3\'-linkage', phosRes.meta.name),
      kv('HELM', nt.raw, true),
    );
  } else {
    const conj = resolveConjugate(m.symbol);
    title.textContent = conj.meta.name;
    detailRows.push(
      kv('Strand', hit.strand === 'sense' ? 'Sense' : 'Antisense'),
      kv('Position', String(hit.position + 1)),
      kv('Type', 'Terminal conjugate / linker'),
      kv('HELM', m.raw, true),
    );
  }

  // Reserved area for structures — populated async.
  // Holds 1 (linkage / conjugate) or 3 (full nucleotide) small canvases side by side.
  const structHost = ui.div([], {style: {
    marginTop: '8px',
    display: 'flex', flexDirection: 'row', gap: '6px',
    flexWrap: 'wrap',
  }});
  structHost.dataset.role = 'struct';
  const placeholder = ui.divText('Loading structures…', {style: {
    fontSize: '11px', color: 'var(--grey-4, #888)',
  }});
  structHost.appendChild(placeholder);

  return ui.divV([title, ...detailRows, structHost], {
    style: {fontSize: '12px', lineHeight: '1.45', minWidth: '240px', maxWidth: '380px'},
  });
}

function kv(k: string, v: string, mono = false): HTMLElement {
  return ui.divH([
    ui.divText(k, {style: {color: 'var(--grey-4, #888)', width: '90px', flexShrink: '0'}}),
    ui.divText(v, {style: mono ? {fontFamily: 'ui-monospace, Menlo, monospace', wordBreak: 'break-all'} : {}}),
  ], {style: {gap: '6px'}});
}

/** Cached merged Bio monomer library. */
let _bioLibPromise: Promise<IMonomerLib> | null = null;
function getBioLib(): Promise<IMonomerLib> {
  if (!_bioLibPromise) {
    _bioLibPromise = (async () => {
      const helper = await getMonomerLibHelper();
      return helper.getMonomerLib();
    })();
  }
  return _bioLibPromise;
}

/** Try the original symbol first, then the canonical alias, across all polymer types. */
function findMonomerMolfile(lib: IMonomerLib, rawSymbol: string, kind: StructureKind): string | null {
  const candidates: string[] = [rawSymbol];
  if (kind === 'sugar') candidates.push(canonicalSugarSymbol(rawSymbol));
  else if (kind === 'phosphate') candidates.push(canonicalPhosphateSymbol(rawSymbol));
  // Bases (A/C/G/U/T) and conjugates: try the symbol as-is.

  const polymerTypes = lib.getPolymerTypes();
  for (const sym of candidates) {
    for (const pt of polymerTypes) {
      const monomer = lib.getMonomer(pt, sym);
      if (monomer && monomer.molfile) return monomer.molfile;
    }
  }
  return null;
}

async function renderStructuresAsync(root: HTMLElement, hit: HitResult): Promise<void> {
  const host = root.querySelector('[data-role="struct"]') as HTMLDivElement | null;
  if (!host) return;

  // Determine which structures to show
  let requests: { label: string; symbol: string; kind: StructureKind }[];
  if (hit.linkage) {
    requests = [{label: 'Linkage', symbol: hit.linkage.phosphateSymbol, kind: 'phosphate'}];
  } else if (hit.monomer.kind === 'nucleotide') {
    const nt = hit.monomer as ParsedNucleotide;
    requests = [
      {label: 'Sugar', symbol: nt.sugar, kind: 'sugar'},
    ];
    if (nt.base) requests.push({label: 'Base', symbol: nt.base, kind: 'base'});
    if (nt.phosphate) requests.push({label: '3\'-linkage', symbol: nt.phosphate, kind: 'phosphate'});
  } else {
    requests = [{label: 'Conjugate', symbol: hit.monomer.symbol, kind: 'conjugate'}];
  }

  let lib: IMonomerLib;
  try {
    lib = await getBioLib();
  } catch {
    host.innerHTML = '';
    host.appendChild(ui.divText('Monomer library unavailable', {
      style: {fontSize: '11px', color: 'var(--grey-4, #888)'},
    }));
    return;
  }

  host.innerHTML = '';

  for (const req of requests) {
    const cell = makeStructureCell(req.label);
    host.appendChild(cell.root);
    const molfile = findMonomerMolfile(lib, req.symbol, req.kind);
    if (!molfile) {
      cell.body.textContent = `Not in HELMCore: ${req.symbol}`;
      continue;
    }
    void drawMolfileCached(cell.body, molfile, cacheKeyFor(req));
  }
}

/** Cache key for a structure render. */
function cacheKeyFor(req: {kind: StructureKind; symbol: string}): string {
  let canonical = req.symbol;
  if (req.kind === 'sugar') canonical = canonicalSugarSymbol(req.symbol);
  else if (req.kind === 'phosphate') canonical = canonicalPhosphateSymbol(req.symbol);
  return `${req.kind}:${canonical}`;
}

function makeStructureCell(label: string): { root: HTMLElement; body: HTMLDivElement } {
  const body = ui.div([], {style: {
    // width: `${STRUCT_W}px`, height: `${STRUCT_H}px`,
    // background: 'var(--grey-1, #f4f4f4)', border: '1px solid var(--grey-2, #ddd)',
    borderRadius: '4px',
    display: 'flex', alignItems: 'center', justifyContent: 'center',
    fontSize: '10px', color: 'var(--grey-4, #888)', textAlign: 'center', padding: '4px',
  }}) as HTMLDivElement;
  const cap = ui.divText(label, {style: {
    fontSize: '10px', color: 'var(--grey-5, #555)', textAlign: 'center', marginTop: '2px',
    letterSpacing: '0.04em', textTransform: 'uppercase', fontWeight: '600',
  }});
  const root = ui.divV([body, cap], {style: {display: 'inline-flex'}});
  return {root, body};
}

/** Cache of rendered structures by kind+canonical-symbol.
 *
 * We cache the actual element returned by grok.chem.drawMolecule and reuse
 * the same instance across tooltips. We deliberately do NOT clone — when
 * `drawMolecule` returns a canvas, its bitmap is painted asynchronously
 * (RDKit init etc.), so a clone taken right after the call ends up blank.
 * Since only one tooltip is visible at a time, moving the cached node
 * between tooltip hosts via `appendChild` is safe: the previous host is
 * being torn down. */
const _structCache = new Map<string, HTMLElement>();

function drawMolfileCached(host: HTMLDivElement, molfile: string, cacheKey: string): void {
  try {
    let cached = _structCache.get(cacheKey);
    if (!cached) {
      cached = grok.chem.drawMolecule(molfile, STRUCT_W, STRUCT_H, false);
      if (_structCache.size > 256) _structCache.clear();
      _structCache.set(cacheKey, cached);
    }
    host.innerHTML = '';
    // appendChild moves the node from any previous parent automatically.
    host.appendChild(cached);
  } catch {
    host.textContent = 'Render failed';
  }
}
