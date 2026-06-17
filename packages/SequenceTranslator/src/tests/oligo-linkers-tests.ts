import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {parseHelmDuplex, canonicalizeHelm} from '../oligo-renderer/helm-parser';
import {computeLayout} from '../oligo-renderer/canvas-renderer';
import {isLinkerSymbol} from '../oligo-renderer/types';
import {isLinkerMonomer} from '../oligo-renderer/alignment';
import {LINKER_SS, LINKER_BULGE} from './oligo-fixtures';

/* ---------------------------------------------------------------- *
 * Standalone backbone linkers (`p`, `[sp]`, …): leading / consecutive / middle
 * / trailing runs are recognized as linkers (not conjugates), rendered as arcs
 * (no chips), and participate in alignment so an asymmetric run opens a gap
 * (with a widened arc on the strand that lacks it).
 * ---------------------------------------------------------------- */

category('OligoRenderer: standalone linkers — parsing', () => {
  test('isLinkerSymbol recognizes phosphate symbols', async () => {
    expect(isLinkerSymbol('p'), true);
    expect(isLinkerSymbol('sp'), true);
    expect(isLinkerSymbol('P'), true); // legacy alias → p
    expect(isLinkerSymbol('sP'), true); // legacy alias → sp
    expect(isLinkerSymbol('s2p'), true);
    expect(isLinkerSymbol('GalNAc'), false);
    expect(isLinkerSymbol('Chol'), false);
    expect(isLinkerSymbol(''), false);
  });

  test('leading standalone phosphates parse as linkers', async () => {
    const m = parseHelmDuplex(LINKER_SS);
    expect(m.sense.monomers[0].kind, 'linker');
    expect(m.sense.monomers[1].kind, 'linker');
    expect(m.sense.monomers[2].kind, 'linker');
    // first nucleotide follows the 3-phosphate cap
    expect(m.sense.monomers[3].kind, 'nucleotide');
  });

  test('consecutive mid-strand linkers all parse as linkers', async () => {
    const m = parseHelmDuplex(LINKER_SS);
    // monomers 12,13,14 are the three consecutive [sp] (3 lead + 9 nt before them)
    expect(m.sense.monomers[12].kind, 'linker');
    expect(m.sense.monomers[13].kind, 'linker');
    expect(m.sense.monomers[14].kind, 'linker');
    expect(m.sense.monomers.filter((x) => x.kind === 'linker').length, 6); // 3 cap + 3 mid
    expect(m.sense.monomers.filter((x) => x.kind === 'nucleotide').length, 19);
  });

  test('trailing standalone linker parses as a linker', async () => {
    const m = parseHelmDuplex('RNA1{r(C)p.r(A)p.r(G).[sp]}$$$$');
    const last = m.sense.monomers[m.sense.monomers.length - 1];
    expect(last.kind, 'linker');
    if (last.kind === 'linker') expect(last.symbol, 'sp');
  });

  test('a base-less conjugate (GalNAc) is NOT treated as a linker', async () => {
    const m = parseHelmDuplex('RNA1{r(C)p.r(A).[GalNAc]}$$$$');
    const last = m.sense.monomers[m.sense.monomers.length - 1];
    expect(last.kind, 'conjugate');
    expect(isLinkerMonomer(last), false);
  });

  test('isLinkerMonomer is true for parsed linker monomers', async () => {
    const m = parseHelmDuplex(LINKER_SS);
    expect(isLinkerMonomer(m.sense.monomers[0]), true); // a `p`
    expect(isLinkerMonomer(m.sense.monomers[12]), true); // an `sp`
    expect(isLinkerMonomer(m.sense.monomers[3]), false); // a nucleotide
  });

  test('canonicalizeHelm preserves standalone linkers as backbone monomers', async () => {
    const out = canonicalizeHelm(LINKER_SS);
    expect(out.includes('p.p.p'), true, `expected the 3-phosphate cap, got: ${out}`);
    expect(out.includes('[sp].[sp].[sp]'), true, `expected the consecutive sp run, got: ${out}`);
  });
});

category('OligoRenderer: standalone linkers — single-strand rendering', () => {
  test('standalone linkers render as arcs, never as chips', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_SS));
    expect(layout.textOnlyFallback, false);
    // every chip is a nucleotide — no linker (and no conjugate) chips
    expect(layout.senseChips.every((c) => c.monomer.kind === 'nucleotide'), true);
    expect(layout.senseChips.length, 19);
  });

  test('every leading / middle / trailing linker becomes an arc', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_SS));
    // 3 cap arcs + 18 inter-base gaps (one carrying the 3 extra sp) =
    // 3 + (17 gaps × 1 phosphate) + (1 gap × (phosphate + 3 sp)) = 3 + 17 + 4 = 24
    expect(layout.senseLinks.length, 24);
  });

  test('consecutive linkers sit in a row (monotonically increasing X)', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_SS));
    // the 3 leading cap arcs are the left-most links; their centers increase
    const leadXs = layout.senseLinks
      .map((l) => l.x + l.w / 2)
      .sort((a, b) => a - b)
      .slice(0, 3);
    expect(leadXs[0] < leadXs[1] && leadXs[1] < leadXs[2], true);
  });

  test('leading arcs flow into the first monomer (not too far left)', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_SS));
    const firstChip = layout.senseChips[0];
    const centers = layout.senseLinks.map((l) => l.x + l.w / 2).sort((a, b) => a - b);
    // the right-most leading cap arc sits just before the first chip, exactly
    // where a normal phosphate arc sits (center ≈ chipLeft − chipGap/2)
    expect(Math.abs(centers[2] - (firstChip.x - layout.chipGap / 2)) < 0.6, true,
      `rightmost lead arc=${centers[2]} expected≈${firstChip.x - layout.chipGap / 2}`);
  });
});

category('OligoRenderer: standalone linkers — duplex alignment & gaps', () => {
  test('a sense-only linker bulge keeps the duplex blunt and bases aligned', async () => {
    const m = parseHelmDuplex(LINKER_BULGE);
    expect(m.sense.monomers.filter((x) => x.kind === 'linker').length, 2);
    const layout = computeLayout(1100, 90, m);
    expect(layout.shift, 0);
    expect(layout.senseChips.length, 19);
    expect(layout.antiChips.length, 19);
    // paired bases past the bulge still share their column X exactly
    const s12 = layout.senseChips.find((c) => c.col === 12)!;
    const a12 = layout.antiChips.find((c) => c.col === 12)!;
    expect(Math.abs(s12.x - a12.x) < 0.01, true,
      `bulge broke alignment at col 12: sense=${s12.x}, anti=${a12.x}`);
  });

  test('the linker bulge widens its gap relative to a normal gap', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_BULGE));
    const sx = (c: number): {x: number; w: number} => layout.senseChips.find((k) => k.col === c)!;
    const gapAfter = (c: number): number => sx(c + 1).x - (sx(c).x + sx(c).w);
    // the 2 standalone linkers live between base 9 and base 10
    expect(gapAfter(9) > gapAfter(0) + 1, true,
      `linker gap not widened: gap@9=${gapAfter(9).toFixed(1)} gap@0=${gapAfter(0).toFixed(1)}`);
  });

  test('the antisense draws a single wider arc across the bulge', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_BULGE));
    // 19-nt antisense → 18 inter-base arcs, none split by standalone linkers
    expect(layout.antiLinks.length, 18);
    const widths = layout.antiLinks.map((l) => l.w);
    const max = Math.max(...widths); const min = Math.min(...widths);
    // the arc spanning the gap that sense fills with 2 linkers is much wider
    expect(max > min + 10, true, `expected a widened antisense arc: max=${max.toFixed(1)} min=${min.toFixed(1)}`);
  });

  test('sense draws an arc per linker monomer in the bulge', async () => {
    const layout = computeLayout(1100, 90, parseHelmDuplex(LINKER_BULGE));
    // 18 inter-base gaps, one of which carries 2 extra standalone-linker arcs → 20
    expect(layout.senseLinks.length, 20);
    // none of the standalone linkers became a chip
    expect(layout.senseChips.every((c) => c.monomer.kind === 'nucleotide'), true);
  });
});
