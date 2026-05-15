import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {parseHelmDuplex, looksLikeHelm, canonicalizeHelm} from '../oligo-renderer/helm-parser';
import {ParsedNucleotide} from '../oligo-renderer/types';
import {computeLayout, hitTest, drawDuplex} from '../oligo-renderer/canvas-renderer';
import {
  resolveSugar, resolvePhosphate, resolveConjugate,
  hashColor, SUGAR_MODS, PHOSPHATE_MODS, CONJUGATE_MODS,
} from '../oligo-renderer/types';

// HELMCore canonical: lowercase backbone (`r`, `m`, `d`, `lna`, `fl2r`, `p`, `sp`).
// Single-char sugars/phosphates are unbracketed; multi-char must be bracketed.
const SAMPLE_DUPLEX =
  'RNA1{m(G)[sp].m(A)[sp].[fl2r](C)p.[fl2r](U)p.r(G)p.r(A)p.r(A)p.r(U)p.r(A)p.r(U)p.r(A)p.r(A)p.r(A)p.' +
  'r(C)p.r(U)p.r(U)p.m(G)p.m(U)[sp].m(G)[sp].[L3]}|' +
  'RNA2{[fl2r](C)[sp].m(A)[sp].m(A)p.m(G)p.r(U)p.r(U)p.r(U)p.r(A)p.r(U)p.r(A)p.r(U)p.r(U)p.r(C)p.' +
  'r(A)p.r(U)p.m(C)p.m(A)[sp].m(G)[sp].r(U)}$$$$';

const SAMPLE_SS = 'RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(C)[sp].d(A)[sp]}$$$$';

// Legacy / Pistoia-style symbols that should still render via the alias map.
const LEGACY_DUPLEX =
  'RNA1{[mR](G)[sP].[mR](A)[sP].[fR](C)P.R(U)P}|RNA2{[mR](C)[sP].R(A)P.R(U)P.R(G)P}$$$$';

category('OligoRenderer: parser', () => {
  test('looksLikeHelm — positive', async () => {
    expect(looksLikeHelm(SAMPLE_DUPLEX), true);
    expect(looksLikeHelm('RNA1{R(A)P}'), true);
    expect(looksLikeHelm('PEPTIDE1{A.C.G}'), true);
  });

  test('looksLikeHelm — negative', async () => {
    expect(looksLikeHelm(''), false);
    expect(looksLikeHelm('AUCGUACGUAGCUAGCAU'), false);
    expect(looksLikeHelm('Af Cf Gf Uf'), false);
  });

  test('parses two-strand duplex', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    expect(m.sense.type, 'RNA');
    expect(m.sense.monomers.length, 20); // 19 nucleotides + L3 conjugate
    expect(m.antisense !== null, true);
    expect(m.antisense!.monomers.length, 19);
  });

  test('parses single-strand', async () => {
    const m = parseHelmDuplex(SAMPLE_SS);
    expect(m.sense.monomers.length, 6);
    expect(m.antisense, null);
    // First nucleotide has LNA sugar, base C, PS phosphate
    const first = m.sense.monomers[0];
    expect(first.kind, 'nucleotide');
    if (first.kind === 'nucleotide') {
      expect(first.sugar, 'lna');
      expect(first.base, 'C');
      expect(first.phosphate, 'sp');
    }
  });

  test('single-char unbracketed deoxyribose `d` parses correctly', async () => {
    // Regression: previously `dR(T)` was emitted (invalid HELM); the parser
    // would split `d` as the sugar and choke on the rest. With canonical `d(T)`,
    // sugar=d, base=T, phosphate empty.
    const m = parseHelmDuplex('RNA1{d(T)p.d(C)p.d(G)}$$$$');
    expect(m.sense.monomers.length, 3);
    const first = m.sense.monomers[0];
    if (first.kind === 'nucleotide') {
      expect(first.sugar, 'd');
      expect(first.base, 'T');
      expect(first.phosphate, 'p');
    }
  });

  test('parses bracketed conjugate', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const last = m.sense.monomers[m.sense.monomers.length - 1];
    expect(last.kind, 'conjugate');
    if (last.kind === 'conjugate') expect(last.symbol, 'L3');
  });

  test('handles empty / malformed without throwing', async () => {
    const m1 = parseHelmDuplex('');
    expect(m1.sense.monomers.length, 0);
    expect(m1.antisense, null);
    const m2 = parseHelmDuplex('not a helm string');
    expect(m2.sense.monomers.length, 0);
  });

  test('passes through unknown sugar/phosphate symbols', async () => {
    const helm = 'RNA1{[xr](A)[zp].[yr](C)p}$$$$';
    const m = parseHelmDuplex(helm);
    const m0 = m.sense.monomers[0];
    if (m0.kind === 'nucleotide') {
      expect(m0.sugar, 'xr');
      expect(m0.phosphate, 'zp');
    }
  });

  test('parser strips brackets from multi-char base in parens', async () => {
    // Bug: pre-fix, cleanupHelmSymbol was not called on the base, so base === '[5Br-dC]'.
    // Post-fix: cleanupHelmSymbol('[5Br-dC]') → '5Br-dC'.
    const dup = parseHelmDuplex('RNA1{r([5Br-dC])p}$$$$');
    const m = dup.sense.monomers[0] as ParsedNucleotide;
    expect(m.kind, 'nucleotide');
    expect(m.sugar, 'r');
    // pre-fix: '[5Br-dC]'  post-fix: '5Br-dC'
    expect(m.base, '5Br-dC');
    expect(m.phosphate, 'p');
  });

  test('parser strips brackets from multi-char base on both strands', async () => {
    // Covers both strands and different modification combinations.
    // pre-fix: base === '[5meC]', '[5Br-dC]', '[5fU]'  post-fix: '5meC', '5Br-dC', '5fU'
    const dup = parseHelmDuplex('RNA1{r([5meC])p.[fl2r]([5Br-dC])[sp]}|RNA2{r([5fU])p}$$$$');
    const senseFirst = dup.sense.monomers[0] as ParsedNucleotide;
    const senseSecond = dup.sense.monomers[1] as ParsedNucleotide;
    const antiFirst = dup.antisense!.monomers[0] as ParsedNucleotide;
    // pre-fix: '[5meC]'  post-fix: '5meC'
    expect(senseFirst.base, '5meC');
    expect(senseSecond.sugar, 'fl2r');
    // pre-fix: '[5Br-dC]'  post-fix: '5Br-dC'
    expect(senseSecond.base, '5Br-dC');
    expect(senseSecond.phosphate, 'sp');
    // pre-fix: '[5fU]'  post-fix: '5fU'
    expect(antiFirst.base, '5fU');
  });

  test('parser leaves bare-letter base unchanged', async () => {
    // Regression guard: single-letter bases must NOT be transformed.
    // A, G, C, T, U are all single-char and do not have brackets → unchanged by cleanupHelmSymbol.
    const dup = parseHelmDuplex('RNA1{r(A)p.r(G)p.r(C)p.r(T)p.r(U)p}$$$$');
    const bases = dup.sense.monomers.map((m) => (m as ParsedNucleotide).base);
    // Expected: ['A', 'G', 'C', 'T', 'U'] — unchanged for all five bases.
    expect(bases.join(','), 'A,G,C,T,U');
  });

  test('canonicalizeHelm re-brackets multi-char base on output', async () => {
    // serializeCanonicalMonomer emits `([${base}])` when base.length > 1.
    // pre-fix: emitted '(5Br-dC)' (invalid HELM, missing brackets).
    // post-fix: emits '([5Br-dC])' — round-trip is valid HELM.
    const out = canonicalizeHelm('RNA1{r([5Br-dC])p}$$$$');
    // Must contain the bracketed base form.
    expect(out.includes('([5Br-dC])'), true,
      `expected canonicalized HELM to contain '([5Br-dC])', got: ${out}`);
    // Must NOT contain the bare (unbracketed) form.
    expect(out.includes('(5Br-dC)') && !out.includes('([5Br-dC])'), false,
      `must not emit unbracketed (5Br-dC) without wrapping brackets, got: ${out}`);
  });

  test('canonicalizeHelm keeps single-letter base unbracketed', async () => {
    // serializeCanonicalMonomer: base.length === 1 → `(${base})`, not `([${base}])`.
    // Expected output contains 'r(A)p' — single letter stays bare.
    const out = canonicalizeHelm('RNA1{r(A)p}$$$$');
    expect(out.includes('r(A)p'), true,
      `expected single-letter base to stay unbracketed, got: ${out}`);
    // Single-letter base must NOT get double-bracketed.
    expect(out.includes('([A])'), false,
      `single-letter base must NOT be wrapped in brackets, got: ${out}`);
  });
});

category('OligoRenderer: modification dictionary', () => {
  test('resolves HELMCore canonical sugars', async () => {
    expect(resolveSugar('m', 'A').meta.short, '2\'-OMe');
    expect(resolveSugar('fl2r', 'C').meta.short, '2\'-F');
    expect(resolveSugar('lna', 'G').meta.short, 'LNA');
    expect(resolveSugar('moe', 'A').meta.short, '2\'-MOE');
  });

  test('aliases legacy uppercase symbols to canonical', async () => {
    // `mR` (legacy) and `m` (canonical) must resolve identically
    expect(resolveSugar('mR', 'A').meta.name, resolveSugar('m', 'A').meta.name);
    expect(resolveSugar('fR', 'C').meta.name, resolveSugar('fl2r', 'C').meta.name);
    expect(resolveSugar('LR', 'G').meta.name, resolveSugar('lna', 'G').meta.name);
    expect(resolveSugar('dR', 'T').meta.name, resolveSugar('d', 'T').meta.name);
  });

  test('uses base-canonical color for unmodified ribose', async () => {
    const a = resolveSugar('r', 'A');
    const c = resolveSugar('r', 'C');
    expect(a.color !== c.color, true,
      `Unmodified A and C must use distinct canonical colors, got ${a.color} for both`);
  });

  test('falls back to hash color for unknown sugar', async () => {
    const x = resolveSugar('zzz', 'A');
    expect(x.color.startsWith('hsl('), true);
    expect(x.meta.name, 'zzz');
  });

  test('hashColor is deterministic', async () => {
    expect(hashColor('foo'), hashColor('foo'));
    expect(hashColor('foo') !== hashColor('bar'), true);
  });

  test('PS phosphate (canonical and legacy) resolve to same color', async () => {
    expect(resolvePhosphate('sp').color, PHOSPHATE_MODS['sp'].color);
    expect(resolvePhosphate('sP').color, PHOSPHATE_MODS['sp'].color);
  });

  test('GalNAc / L3 / Chol resolve as conjugates', async () => {
    expect(resolveConjugate('GalNAc').meta.category, 'conjugate');
    expect(resolveConjugate('L3').meta.category, 'conjugate');
    expect(resolveConjugate('Chol').meta.color, CONJUGATE_MODS['Chol'].color);
  });
});

category('OligoRenderer: layout', () => {
  test('comfortable cell — chips fit width and height', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(500, 70, m);
    expect(layout.textOnlyFallback, false);
    expect(layout.chipW > 5, true, 'chip width should be larger than minimum');
    // chipW is bounded above by the cell's height budget. V_PAD=10 is the
    // top/bottom breathing room; the rest is divided among 2 strands +
    // 1 strand-gap + 2 apex zones with the fixed chip aspect ratio.
    const heightBudget = (70 - 2 * 10) / (2 + 0.5 + 0.9) / 1.25;
    expect(layout.chipW <= heightBudget + 0.01, true,
      `chipW should be bounded by height (got ${layout.chipW}, cap ${heightBudget})`);
    expect(layout.antiY > layout.senseY, true);
    expect(layout.senseChips.length > 0, true);
    expect(layout.antiChips.length > 0, true);
    // Top apex must sit at least V_PAD=10 from the cell's top edge.
    expect(layout.senseY - layout.apexH >= 10 - 0.01, true,
      `top apex zone violates 10px top margin: senseY=${layout.senseY}, apexH=${layout.apexH}`);
  });

  test('chipW grows with cell height (no hard max cap)', async () => {
    // Same duplex, taller cells. Without a hardcoded MAX_CHIP_W cap, chipW
    // must grow proportionally with cellH (height budget governs the cap).
    // Smallest cell here is 60px to stay above the V_PAD=10 floor where
    // the height budget would otherwise dip under MIN_CHIP_W.
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layoutSmall = computeLayout(1000, 60, m);
    const layoutMid   = computeLayout(1000, 100, m);
    const layoutBig   = computeLayout(1000, 180, m);
    expect(layoutMid.chipW > layoutSmall.chipW, true,
      `expected chipW to grow with height: ${layoutSmall.chipW} → ${layoutMid.chipW}`);
    expect(layoutBig.chipW > layoutMid.chipW, true,
      `expected chipW to keep growing: ${layoutMid.chipW} → ${layoutBig.chipW}`);
    // Previously this was capped at MAX_CHIP_W=17. With the cap removed and
    // a 180px-tall cell, chipW comfortably exceeds 30 even after V_PAD=10
    // takes 20px out of the vertical budget.
    expect(layoutBig.chipW > 30, true,
      `chipW should exceed previous 17px cap on tall cells: got ${layoutBig.chipW}`);
  });

  test('chipW grows with cell width when not height-limited', async () => {
    // Wide+tall cell with a small duplex (4 chips per strand). Width budget
    // is generous, so chipW should be governed by the height cap.
    const helm = 'RNA1{m(A)p.m(C)p.m(G)p.m(U)}|RNA2{m(U)p.m(G)p.m(C)p.m(A)}$$$$';
    const m = parseHelmDuplex(helm);
    const tallWide = computeLayout(2000, 140, m);
    // Should fit easily; chipW capped by height budget. (cellH - 2*V_PAD)
    // / heightFactor / aspect.
    const heightCap = (140 - 2 * 10) / (2 + 0.5 + 0.9) / 1.25;
    expect(Math.abs(tallWide.chipW - heightCap) < 0.5, true,
      `tall+wide cell: chipW should ride the height cap (got ${tallWide.chipW}, expected ~${heightCap.toFixed(1)})`);
  });

  test('top apex and bottom apex sit at least 10px from cell edges', async () => {
    // V_PAD=10 guarantees vertical breathing room so the apex tips never
    // touch the cell border.
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    for (const [w, h] of [[600, 70], [1000, 120], [1500, 200]]) {
      const layout = computeLayout(w, h, m);
      if (layout.textOnlyFallback) continue;
      const topApexEdge = layout.senseY - layout.apexH;
      const bottomApexEdge = layout.antiY >= 0 ? layout.antiY + layout.chipH + layout.apexH : -1;
      expect(topApexEdge >= 10 - 0.01, true,
        `top apex too close to top at ${w}×${h}: ${topApexEdge}`);
      if (bottomApexEdge > 0)
        expect(h - bottomApexEdge >= 10 - 0.01, true,
          `bottom apex too close to bottom at ${w}×${h}: cellH-edge=${h - bottomApexEdge}`);
    }
  });

  test('multi-char base ellipsizes — chipW does NOT shrink for everyone else', async () => {
    // Tight cell with one long custom base on sense. Show-everything priority
    // means: keep chipW as big as uniform-fit allows; ellipsize the long
    // multi-char chip rather than shrinking every chip on the row.
    const canonical = parseHelmDuplex(
      'RNA1{m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)}|' +
      'RNA2{m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)}$$$$');
    const multichar = parseHelmDuplex(
      'RNA1{m([cpm6A])p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)}|' +
      'RNA2{m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)}$$$$');
    const cellW = 300; const cellH = 50;
    const lCanon = computeLayout(cellW, cellH, canonical);
    const lMulti = computeLayout(cellW, cellH, multichar);
    // Sizing is governed by uniform-fit, which depends only on chip COUNT
    // and conjugate widths — not on whether a base is single- or multi-char.
    // So a multi-char base must not pull chipW down.
    expect(Math.abs(lMulti.chipW - lCanon.chipW) < 0.5, true,
      `multi-char chipW=${lMulti.chipW.toFixed(2)} drifted from canonical chipW=${lCanon.chipW.toFixed(2)}; ` +
      `sizing should be based on uniform fit, not widened`);
    // Every monomer still draws (just the multi-char one shows an ellipsized label).
    expect(lMulti.senseChips.length, multichar.sense.monomers.length);
    expect(lMulti.antiChips.length, multichar.antisense!.monomers.length);
  });

  test('leading conjugate does NOT cause chipW to stall on roomy cells', async () => {
    // Regression: `fitsInBudget` used to double-count the leading conjugate
    // width (added `leadW` to the running total AND `widths[0]` already
    // contained the pill width), so the binary search bailed early on cells
    // where the conjugate strand had enough room. Net effect was the layout
    // chose a chipW much smaller than what the cell could actually hold.
    // This is the exact HELM the user reported, on a wide+tall cell.
    const helm =
      'RNA1{[Chol].m(G)[sp].[fl2r](A)[sp].m(C)p.[fl2r](U)p.m(G)p.[fl2r](A)p.m(A)p.[fl2r](U)p.m(A)p.' +
      '[fl2r](U)p.m(A)p.[fl2r](A)p.m(A)p.[fl2r](C)p.m(U)p.[fl2r](U)p.m(G)[sp].[fl2r](U)[sp].m(G)[L3]}|' +
      'RNA2{m(C)[sp].[fl2r]([br8A])[sp].m(C)p.[fl2r](A)p.m(A)p.[fl2r](G)p.m(U)p.[fl2r](U)p.m(U)p.' +
      '[fl2r](A)p.m(U)p.[fl2r](A)p.m(U)p.[fl2r]([cneT])p.m(C)p.[fl2r](A)p.m(G)[sp].' +
      '[fl2r]([c7io7A])[sp].m(C)}$$$$';
    const m = parseHelmDuplex(helm);
    const layoutNoConj = computeLayout(1000, 120, parseHelmDuplex(
      'RNA1{m(G)p.m(A)p.m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)p.m(U)p.m(G)}|' +
      'RNA2{m(C)p.m(A)p.m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)p.m(U)p.m(C)}$$$$'));
    const layoutWithConj = computeLayout(1000, 120, m);
    // With a 1000×120 cell, the height cap allows chipW ~28. The version
    // without conjugates rides that cap; the version with [Chol]/[L3] +
    // multi-char bases should pick a chipW that's a meaningful fraction
    // (≥ 60%) of that — not collapse back to single digits.
    expect(layoutWithConj.chipW >= layoutNoConj.chipW * 0.6, true,
      `with-conjugate chipW=${layoutWithConj.chipW.toFixed(2)} dropped too far ` +
      `below no-conjugate chipW=${layoutNoConj.chipW.toFixed(2)} on the same roomy cell`);
    // And every monomer must place.
    expect(layoutWithConj.senseChips.length, m.sense.monomers.length);
    expect(layoutWithConj.antiChips.length, m.antisense!.monomers.length);
  });

  test('leading conjugate reduces chipW but no monomer is dropped', async () => {
    // With a Cholesterol conjugate at sense 5', the conjugate pill consumes
    // a chunk of horizontal budget. chipW shrinks compared to the same data
    // without the conjugate — but all chips still place.
    const noConj = parseHelmDuplex(
      'RNA1{m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)}|' +
      'RNA2{m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)}$$$$');
    const withConj = parseHelmDuplex(
      'RNA1{[Chol].m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)}|' +
      'RNA2{m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)}$$$$');
    const cellW = 350; const cellH = 70;
    const lNo = computeLayout(cellW, cellH, noConj);
    const lW = computeLayout(cellW, cellH, withConj);
    expect(lW.chipW <= lNo.chipW + 0.01, true,
      `with-conjugate chipW=${lW.chipW} must be ≤ no-conjugate chipW=${lNo.chipW}`);
    // Every monomer (incl. the conjugate pill) is placed.
    expect(lW.senseChips.length, withConj.sense.monomers.length);
    expect(lW.antiChips.length, withConj.antisense!.monomers.length);
  });

  test('tight cell falls back to uniform widths with ellipsis (no text-only)', async () => {
    // A duplex with a long base name on a moderately tight cell: widened
    // would overflow even at MIN_CHIP_W, but uniform widths (with the
    // multi-char chip ellipsized) still fit.
    const m = parseHelmDuplex(
      'RNA1{m([cpm6A])p.m([5Br-dC])p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)p.m(U)p.m(A)p.m(C)p.m(G)}|' +
      'RNA2{m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)p.m(A)p.m(U)p.m(G)p.m(C)}$$$$');
    const layout = computeLayout(220, 50, m);
    // Not text-only — everything still draws as chips, just with truncated labels.
    expect(layout.textOnlyFallback, false,
      'tight-but-survivable cell should not collapse to text-only');
    // All chips placed (no chip dropped).
    expect(layout.senseChips.length, m.sense.monomers.length);
    expect(layout.antiChips.length, m.antisense!.monomers.length);
  });

  test('compact cell — both strands still fit in 28px row', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(300, 28, m);
    expect(layout.antiY > 0, true);
    expect(layout.chipH * 2 + layout.strandGap <= 28, true,
      'two strands must fit within 28px row');
  });

  test('very narrow cell triggers text-only fallback', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(60, 30, m);
    expect(layout.textOnlyFallback, true);
    expect(layout.senseChips.length, 0);
  });

  test('aspect ratio is preserved across sizes', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    for (const [w, h] of [[300, 50], [500, 70], [800, 90]]) {
      const layout = computeLayout(w, h, m);
      const ratio = layout.chipH / layout.chipW;
      // ASPECT_H_OVER_W = 1.25 in canvas-renderer.ts
      expect(Math.abs(ratio - 1.25) < 0.01, true,
        `aspect ratio drifted at ${w}x${h}: got ${ratio.toFixed(3)}`);
    }
  });

  test('single-strand layout', async () => {
    const m = parseHelmDuplex(SAMPLE_SS);
    const layout = computeLayout(400, 50, m);
    expect(layout.antiY, -1);
    expect(layout.antiChips.length, 0);
  });

  test('antisense reversed by default for pair-alignment', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(800, 70, m);
    expect(layout.antiReversed, true);
    // First antisense chip in display = last antisense monomer in data
    const antiLast = m.antisense!.monomers[m.antisense!.monomers.length - 1];
    expect(layout.antiChips[0].monomer === antiLast, true);
  });

  test('pairAlign: false leaves antisense in data order', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(800, 70, m, {pairAlign: false});
    expect(layout.antiReversed, false);
    expect(layout.antiChips[0].monomer === m.antisense!.monomers[0], true);
  });

  test('conjugate occupies a wider position; following chip does not overlap', async () => {
    // KRAS-like duplex with [L3] conjugate at sense 3'
    const helm = 'RNA1{m(G)p.m(A)p.m(C)p.[L3]}|RNA2{m(C)p.m(A)p.m(A)p.r(G)}$$$$';
    const m = parseHelmDuplex(helm);
    const layout = computeLayout(600, 70, m);
    const conj = layout.senseChips[3];
    expect(conj.monomer.kind, 'conjugate');
    expect(conj.w >= layout.chipW, true, 'conjugate must be at least chip-wide');
  });

  test('reversed-antisense linkage owner is the lower-indexed pair member', async () => {
    // The `phosphate` field on a nucleotide always means "linkage immediately
    // AFTER this monomer in 5'→3' data order" — so it lives on the lower-indexed
    // end of the bond. When antisense is displayed reversed, the gap to the
    // right of display index i pairs data positions (N-1-i, N-2-i); the owner
    // must be the lower one, i.e. `monomers[i+1]` in the reversed array, not `m`.
    //
    // Small fixture: antisense `r(A)[sp].r(C)p.r(G)` →
    //   data 0 → A, phos=sp   (linkage 0↔1, sp)
    //   data 1 → C, phos=p    (linkage 1↔2, p)
    //   data 2 → G, phos=''   (terminal)
    // Reversed display: G, C, A.
    //   antiLinks[0] = gap right of display 0 → between data 2 and 1 → p, owner 1
    //   antiLinks[1] = gap right of display 1 → between data 1 and 0 → sp, owner 0
    const helm = 'RNA1{r(A)p.r(C)p.r(G)}|RNA2{r(A)[sp].r(C)p.r(G)}$$$$';
    const m = parseHelmDuplex(helm);
    const layout = computeLayout(600, 70, m);
    expect(layout.antiReversed, true);
    expect(layout.antiLinks.length, 2);
    expect(layout.antiLinks[0].ownerOrigIdx, 1);
    expect(layout.antiLinks[0].phosphateSymbol, 'p');
    expect(layout.antiLinks[1].ownerOrigIdx, 0);
    expect(layout.antiLinks[1].phosphateSymbol, 'sp');
  });

  test('reversed-antisense draws ALL linkages including the leftmost-data sp', async () => {
    // Regression: pre-fix, reversed-strand placement read `m.phosphate` for
    // the gap to the right of display i, but that field belongs to the bond
    // on the OTHER side of `m`. The net effect was a one-index shift across
    // the row plus a dropped link at the terminal display position — so the
    // 4 sp linkages on this antisense ended up as 3, in the wrong gaps.
    const helm =
      'RNA1{m(C)[sp].m(A)[sp].m(U)p.m(G)p.m(G)p.m(U)p.m(U)p.m(G)p.m(A)p.m(A)p.' +
      'm(C)p.m(A)p.m(U)p.m(G)p.m(A)p.m(G)p.m(C)[sp].m(A)[sp].m(A)[L3]}|' +
      'RNA2{m(U)[sp].m(U)[sp].m(G)p.m(C)p.m(U)p.m(C)p.m(A)p.m(U)p.m(G)p.m(U)p.' +
      'm(U)p.m(C)p.m(A)p.m(A)p.m(C)p.m(C)p.m(A)[sp].m(U)[sp].m(G)}$$$$';
    const m = parseHelmDuplex(helm);
    const layout = computeLayout(1200, 90, m);
    expect(layout.antiReversed, true);
    // Antisense has 19 nucleotides → 18 inter-nucleotide gaps, all drawn.
    expect(layout.antiLinks.length, 18,
      `expected 18 antisense linkages, got ${layout.antiLinks.length}`);
    // 4 of them must be `sp`. With reversal, the sp at data 0↔1 maps to the
    // RIGHTMOST display gap and the sp at data 17↔18 maps to the LEFTMOST —
    // so sp owners, in display order, are 17, 16, 1, 0.
    const spOwners = layout.antiLinks
      .filter((l) => l.phosphateSymbol === 'sp')
      .map((l) => l.ownerOrigIdx);
    expect(spOwners.length, 4,
      `expected 4 sp linkages on antisense, got ${spOwners.length}`);
    expect(spOwners.join(','), '17,16,1,0',
      `sp owners (display order) must be 17,16,1,0; got ${spOwners.join(',')}`);

    // And sense — 19 nucleotides + L3 conjugate. 18 nucleotide-to-nucleotide
    // gaps; the bond into the L3 conjugate is not a phosphate, so it isn't
    // pushed as a link. 4 of those 18 are sp (data owners 0, 1, 16, 17).
    expect(layout.senseLinks.length, 18,
      `expected 18 sense linkages, got ${layout.senseLinks.length}`);
    const senseSpOwners = layout.senseLinks
      .filter((l) => l.phosphateSymbol === 'sp')
      .map((l) => l.ownerOrigIdx)
      .sort((a, b) => a - b);
    expect(senseSpOwners.join(','), '0,1,16,17');
  });
});

category('OligoRenderer: hit testing', () => {
  test('hits each chip via stored positions', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(800, 70, m);

    const first = layout.senseChips[0];
    const hit = hitTest(first.x + first.w / 2, layout.senseY + layout.chipH / 2, m, layout);
    expect(hit !== null, true);
    expect(hit!.strand, 'sense');
    expect(hit!.position, 0);

    // Gap between chips at chip-row-Y should always miss now — PS linkage
    // markers live in the apex zone above (sense) / below (antisense) the
    // chip row, not in the inter-chip gap at chip-row-Y.
    const gapX = first.x + first.w + layout.chipGap * 0.1;
    const miss = hitTest(gapX, layout.senseY + layout.chipH / 2, m, layout);
    expect(miss, null);
  });

  test('antisense hit returns original position despite reversed display', async () => {
    const m = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = computeLayout(800, 70, m);
    // Click on the leftmost antisense chip in display — that's the LAST monomer in data
    const leftmost = layout.antiChips[0];
    const hit = hitTest(leftmost.x + leftmost.w / 2, layout.antiY + layout.chipH / 2, m, layout);
    expect(hit !== null, true);
    expect(hit!.strand, 'antisense');
    expect(hit!.position, m.antisense!.monomers.length - 1,
      'leftmost AS chip in pair-aligned display = last monomer in data');
  });

  test('hovering the apex above a PS linkage returns the linkage', async () => {
    // PS linkage between sense pos 0 and 1. The apex sits above the sense
    // chip row (since sense's decoration side is "top"), peak at senseY-apexH.
    const helm = 'RNA1{m(G)[sp].m(A)p.m(C)p}$$$$';
    const m = parseHelmDuplex(helm);
    const layout = computeLayout(600, 70, m);
    const c0 = layout.senseChips[0];
    const c1 = layout.senseChips[1];
    const midX = (c0.x + c0.w + c1.x) / 2;
    // Apex zone is [senseY - apexH, senseY]; aim mid-zone.
    const apexY = layout.senseY - layout.apexH / 2;
    const hit = hitTest(midX, apexY, m, layout);
    expect(hit !== null, true);
    expect(hit!.linkage !== undefined, true);
    expect(hit!.linkage!.phosphateSymbol, 'sp');
  });
});

category('OligoRenderer: drawing smoke', () => {
  test('draws without error on offscreen canvas', async () => {
    const canvas = document.createElement('canvas');
    canvas.width = 500; canvas.height = 70;
    const g = canvas.getContext('2d')!;
    const model = parseHelmDuplex(SAMPLE_DUPLEX);
    const layout = drawDuplex(g, 0, 0, 500, 70, model);
    expect(layout.textOnlyFallback, false);
  });

  test('handles unknown modifications without throwing', async () => {
    const helm = 'RNA1{[xr](A)[zp].[yr](C)p}$$$$';
    const canvas = document.createElement('canvas');
    canvas.width = 200; canvas.height = 50;
    const g = canvas.getContext('2d')!;
    const model = parseHelmDuplex(helm);
    drawDuplex(g, 0, 0, 200, 50, model);
    // Verifying no throw is the assertion.
    expect(true, true);
  });

  test('renders legacy uppercase symbols via alias map', async () => {
    const canvas = document.createElement('canvas');
    canvas.width = 400; canvas.height = 50;
    const g = canvas.getContext('2d')!;
    const model = parseHelmDuplex(LEGACY_DUPLEX);
    drawDuplex(g, 0, 0, 400, 50, model);
    expect(model.sense.monomers.length, 4);
    expect(model.antisense !== null, true);
  });
});
