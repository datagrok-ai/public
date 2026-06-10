/**
 * Context-panel widget for an OligoNucleotide cell.
 *
 * - Per-cell summary (lengths, modification counts, conjugates).
 * - Color legend for the column (one-time visual key, not per-cell repeat).
 * - Quick actions: copy HELM, open in pattern designer (TODO).
 *
 * Registered via `tags: panel` + `input: semantic_value oligo {semType: OligoNucleotide}`
 * so the platform shows it automatically when an oligo cell is selected.
 */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {parseHelmDuplex} from './helm-parser';
import {
  BASE_COLORS, FALLBACK_COLOR,
  canonicalPhosphateSymbol, canonicalSugarSymbol,
  isCanonicalBase,
  ParsedNucleotide, resolveConjugate, resolvePhosphate, resolveSugar,
} from './types';
import {getNaturalAnalog} from './analog-cache';
import {getMonomerColors} from './monomer-colors';
import {describeOverhangs} from './alignment';

export function buildOligoPanel(value: DG.SemanticValue): DG.Widget {
  const helm: string = value.value ?? '';
  const model = parseHelmDuplex(helm);

  const sLen = model.sense.monomers.filter((m) => m.kind === 'nucleotide').length;
  const aLen = model.antisense ? model.antisense.monomers.filter((m) => m.kind === 'nucleotide').length : 0;

  // Aggregate modification usage in this oligo
  const sugarCounts = new Map<string, number>();
  const phosCounts = new Map<string, number>();
  const conjCounts = new Map<string, number>();
  /** Custom (non-canonical) base symbols, e.g. `cpm6A`, `5BrU`, `psiU`. */
  const baseCounts = new Map<string, number>();
  const allMonomers = [...model.sense.monomers, ...(model.antisense?.monomers ?? [])];
  for (const m of allMonomers) {
    if (m.kind === 'nucleotide') {
      const nt = m as ParsedNucleotide;
      bump(sugarCounts, nt.sugar);
      if (nt.phosphate) bump(phosCounts, nt.phosphate);
      if (nt.base && !isCanonicalBase(nt.base))
        bump(baseCounts, nt.base);
    } else {
      bump(conjCounts, m.symbol);
    }
  }

  const root = ui.divV([], {style: {fontSize: '12px'}});

  // Summary
  const summary: Record<string, string> = {
    'Sense length': `${sLen} nt`,
    'Antisense length': aLen ? `${aLen} nt` : 'single-strand',
    'Modifications used': humanizeModSet(sugarCounts, phosCounts, baseCounts),
    'Conjugates': conjCounts.size ?
      Array.from(conjCounts.entries()).map(([s, n]) => `${resolveConjugate(s).meta.name} ×${n}`).join(', ') :
      '—',
  };
  if (model.antisense) summary['Duplex'] = describeDuplexAlignment(model);
  root.appendChild(section('Summary', ui.tableFromMap(summary)));

  // Legend — only modifications actually present in this cell.
  // Note: there's no "Copy" section here — the platform already adds a
  // default "Actions | Copy value" entry for any cell.
  root.appendChild(section('Legend', buildCellLegend(sugarCounts, phosCounts, conjCounts, baseCounts)));

  return DG.Widget.fromRoot(root);
}

function bump(map: Map<string, number>, key: string): void {
  map.set(key, (map.get(key) ?? 0) + 1);
}

/** One-line human summary of how the two strands line up: blunt vs the
 * overhangs implied by the alignment, plus where the register came from
 * (explicit HELM pairs vs auto-aligned by complementarity). */
function describeDuplexAlignment(model: ReturnType<typeof parseHelmDuplex>): string {
  const o = describeOverhangs(model);
  if (!o) return '—';
  const ends: string[] = [];
  if (o.sense5) ends.push(`5' sense +${o.sense5}`);
  if (o.anti3) ends.push(`3' antisense +${o.anti3}`);
  if (o.sense3) ends.push(`3' sense +${o.sense3}`);
  if (o.anti5) ends.push(`5' antisense +${o.anti5}`);
  const shape = ends.length ? `overhangs: ${ends.join(', ')}` : 'blunt';
  const src = o.source === 'explicit' ? 'from HELM pairs' :
    o.source === 'auto' ? 'auto-aligned' : '';
  return `${o.paired} bp, ${shape}${src ? ` (${src})` : ''}`;
}

function humanizeModSet(
  sugars: Map<string, number>, phos: Map<string, number>, bases: Map<string, number>,
): string {
  const parts: string[] = [];
  // Aggregate by canonical name so legacy + canonical symbols collapse correctly
  const sugarByName = new Map<string, number>();
  for (const [sym, n] of sugars.entries()) {
    const meta = resolveSugar(sym, null).meta;
    if (meta.short === 'RNA' || meta.short === 'DNA') continue;
    sugarByName.set(meta.short, (sugarByName.get(meta.short) ?? 0) + n);
  }
  for (const [name, n] of sugarByName.entries()) parts.push(`${name} ×${n}`);

  const phosByName = new Map<string, number>();
  for (const [sym, n] of phos.entries()) {
    const meta = resolvePhosphate(sym).meta;
    if (meta.short === 'PO') continue;
    phosByName.set(meta.short, (phosByName.get(meta.short) ?? 0) + n);
  }
  for (const [name, n] of phosByName.entries()) parts.push(`${name} ×${n}`);

  // Custom bases — list each by its raw HELM symbol
  for (const [sym, n] of bases.entries()) parts.push(`${sym} ×${n}`);

  return parts.length ? parts.join(', ') : 'unmodified';
}

function section(title: string, body: HTMLElement): HTMLElement {
  return ui.divV([
    ui.divText(title, {style: {
      fontWeight: '600', marginTop: '8px', marginBottom: '4px',
      color: 'var(--grey-5, #555)', textTransform: 'uppercase', fontSize: '10px', letterSpacing: '0.04em',
    }}),
    body,
  ]);
}

/** Build a legend filtered to modifications actually present in this cell.
 * One row per unique mod, with the count to the right. Items collapse by
 * canonical symbol so legacy/canonical aliases share a row. Colors come from
 * the same `getMonomerColors()` path the canvas renderer uses — so what's
 * drawn on the chip and what's in the legend swatch are always in sync. The
 * local mod meta only supplies the fallback color and human-readable name. */
function buildCellLegend(
  sugars: Map<string, number>, phos: Map<string, number>,
  conjs: Map<string, number>, bases: Map<string, number>,
): HTMLElement {
  type Item = { label: string; color: string; count: number };
  const items: Item[] = [];

  // Sugars — collapse by canonical symbol so e.g. mR + m → one row. Canonical
  // ribose/deoxyribose are included too since the canvas now draws a stripe
  // for every sugar, not just the modified ones.
  const sugarByCanon = new Map<string, number>();
  for (const [sym, n] of sugars.entries()) {
    const c = canonicalSugarSymbol(sym);
    sugarByCanon.set(c, (sugarByCanon.get(c) ?? 0) + n);
  }
  for (const [c, n] of sugarByCanon.entries()) {
    const meta = resolveSugar(c, null).meta;
    const libColor = getMonomerColors('sugar', c).backgroundcolor;
    items.push({label: meta.name, color: libColor ?? meta.color, count: n});
  }

  // Phosphate / linkage mods — show ALL linkages used in the cell (including
  // canonical `p`) since the canvas now draws an apex for every linkage.
  const phosByCanon = new Map<string, number>();
  for (const [sym, n] of phos.entries()) {
    const c = canonicalPhosphateSymbol(sym);
    phosByCanon.set(c, (phosByCanon.get(c) ?? 0) + n);
  }
  for (const [c, n] of phosByCanon.entries()) {
    const meta = resolvePhosphate(c).meta;
    const libColor = getMonomerColors('linker', c).backgroundcolor;
    items.push({label: `${meta.name} (linkage)`, color: libColor ?? meta.color, count: n});
  }

  // Custom bases — prefer the library's base color; fall back to natural-
  // analog palette and finally to FALLBACK_COLOR.
  for (const [sym, n] of bases.entries()) {
    const libColor = getMonomerColors('base', sym).backgroundcolor;
    const analog = getNaturalAnalog(sym);
    const color = libColor ?? (analog && BASE_COLORS[analog]) ?? FALLBACK_COLOR;
    const label = analog ? `${sym} (base, analog ${analog})` : `${sym} (base)`;
    items.push({label, color, count: n});
  }

  // Conjugates
  for (const [sym, n] of conjs.entries()) {
    const meta = resolveConjugate(sym).meta;
    const libColor = getMonomerColors('chem', sym).backgroundcolor;
    items.push({label: meta.name, color: libColor ?? meta.color, count: n});
  }

  if (items.length === 0) {
    return ui.divText('No modifications. Chip color reflects the base (A/C/G/U).', {
      style: {fontSize: '12px', color: 'var(--grey-4, #888)', lineHeight: '1.4'},
    });
  }

  const rows = items.map(({label, color, count}) => ui.divH([
    ui.div([], {style: {
      width: '14px', height: '14px', borderRadius: '3px',
      background: color, border: '1px solid rgba(0,0,0,0.15)', flexShrink: '0',
    }}),
    ui.divText(label, {style: {fontSize: '12px', flex: '1'}}),
    ui.divText(`×${count}`, {style: {
      fontSize: '11px', color: 'var(--grey-4, #888)', fontVariantNumeric: 'tabular-nums',
    }}),
  ], {style: {gap: '6px', alignItems: 'center', marginBottom: '3px'}}));

  return ui.divV(rows);
}
