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
  canonicalPhosphateSymbol, canonicalSugarSymbol,
  ParsedNucleotide, resolveConjugate, resolvePhosphate, resolveSugar,
} from './types';

export function buildOligoPanel(value: DG.SemanticValue): DG.Widget {
  const helm: string = value.value ?? '';
  const model = parseHelmDuplex(helm);

  const sLen = model.sense.monomers.filter((m) => m.kind === 'nucleotide').length;
  const aLen = model.antisense ? model.antisense.monomers.filter((m) => m.kind === 'nucleotide').length : 0;

  // Aggregate modification usage in this oligo
  const sugarCounts = new Map<string, number>();
  const phosCounts = new Map<string, number>();
  const conjCounts = new Map<string, number>();
  const allMonomers = [...model.sense.monomers, ...(model.antisense?.monomers ?? [])];
  for (const m of allMonomers) {
    if (m.kind === 'nucleotide') {
      const nt = m as ParsedNucleotide;
      bump(sugarCounts, nt.sugar);
      if (nt.phosphate) bump(phosCounts, nt.phosphate);
    } else {
      bump(conjCounts, m.symbol);
    }
  }

  const root = ui.divV([], {style: {fontSize: '12px'}});

  // Summary
  root.appendChild(section('Summary', ui.tableFromMap({
    'Sense length': `${sLen} nt`,
    'Antisense length': aLen ? `${aLen} nt` : 'single-strand',
    'Modifications used': humanizeModSet(sugarCounts, phosCounts),
    'Conjugates': conjCounts.size ?
      Array.from(conjCounts.entries()).map(([s, n]) => `${resolveConjugate(s).meta.name} ×${n}`).join(', ') :
      '—',
  })));

  // Legend — only modifications actually present in this cell.
  // Note: there's no "Copy" section here — the platform already adds a
  // default "Actions | Copy value" entry for any cell.
  root.appendChild(section('Legend', buildCellLegend(sugarCounts, phosCounts, conjCounts)));

  return DG.Widget.fromRoot(root);
}

function bump(map: Map<string, number>, key: string): void {
  map.set(key, (map.get(key) ?? 0) + 1);
}

function humanizeModSet(sugars: Map<string, number>, phos: Map<string, number>): string {
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
 * canonical symbol so legacy/canonical aliases share a row. */
function buildCellLegend(
  sugars: Map<string, number>, phos: Map<string, number>, conjs: Map<string, number>,
): HTMLElement {
  type Item = { label: string; color: string; count: number };
  const items: Item[] = [];

  // Collapse by canonical symbol so e.g. mR + m → one row
  const sugarByCanon = new Map<string, number>();
  for (const [sym, n] of sugars.entries()) {
    const c = canonicalSugarSymbol(sym);
    if (c === 'r' || c === 'd') continue; // unmodified
    sugarByCanon.set(c, (sugarByCanon.get(c) ?? 0) + n);
  }
  for (const [c, n] of sugarByCanon.entries()) {
    const meta = resolveSugar(c, null).meta;
    items.push({label: meta.name, color: meta.color, count: n});
  }

  const phosByCanon = new Map<string, number>();
  for (const [sym, n] of phos.entries()) {
    const c = canonicalPhosphateSymbol(sym);
    if (c === 'p') continue;
    phosByCanon.set(c, (phosByCanon.get(c) ?? 0) + n);
  }
  for (const [c, n] of phosByCanon.entries()) {
    const meta = resolvePhosphate(c).meta;
    items.push({label: `${meta.name} (linkage)`, color: meta.color, count: n});
  }

  for (const [sym, n] of conjs.entries()) {
    const meta = resolveConjugate(sym).meta;
    items.push({label: meta.name, color: meta.color, count: n});
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
