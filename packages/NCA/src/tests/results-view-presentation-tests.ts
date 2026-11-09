/**
 * Presentation tests for the per-subject results grid (drift guard) and the
 * aggregated ("By group") grid (GAP-10).
 *
 * The drift guard exists because `C0_EXTRAP` shipped with a raw `C0_EXTRAP`
 * header and DG's generic min/avg/max/std stats tooltip: it was added to the
 * writer schema and to `PRESENTATION_COLUMN_ORDER` but to none of the
 * presentation maps, and nothing failed. The writer drift-guards its own side
 * (`FIXED_RESULT_COLUMNS` vs the built `cols` array); these tests are the
 * presentation-side counterpart.
 *
 * `aggDisplayNames` is the presentation side of a cross-module string contract:
 * the aggregator writes a CV/Mean/SD `label` tag as `"${statWord} ${sourceLabel}"`
 * (e.g. "%CV CL/F"), and `aggDisplayNames` strips the `"${statWord} "` prefix back
 * off to recover the bare parameter ("CL/F") — preserving the apparent-clearance
 * (CL/F) vs true-clearance (CL) distinction in the grouped header. These tests pin
 * that helper directly (the aggregator side is pinned in
 * adapter-results-aggregator-tests.ts) and pin the end-to-end coupling by feeding a
 * real aggregator-produced label through the helper.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {
  aggDisplayNames, missingPresentationEntries,
  colourBandProjections, DEFAULT_RESULTS_THRESHOLDS,
  buildHeaderTooltips, buildC0TooltipSpec, buildBwVariesTooltipSpec,
} from '../views/results-view';
import {richTipEl} from '../views/rich-tooltip';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, indomethContract, indomethScalars, multicompoundContract, multicompoundScalars,
  runPipelineOrThrow, makeSyntheticObs, singleSubjectContract, singleSubjectScalars,
} from './test-utils';

category('Views: results grid presentation completeness', () => {
  test('every presentation column has a friendly header (no raw SCHEMA_NAME in a grid header)', async () => {
    const {noHeader} = missingPresentationEntries();
    expect(noHeader.length, 0,
      `columns missing a DISPLAY_HEADERS entry: [${noHeader.join(', ')}] — ` +
      'each would render its raw UPPERCASE_UNDERSCORE schema name in the grid header ' +
      'next to C_max / t_max / AUC_0–t. Add the friendly name, or add it to ' +
      'HEADER_EXEMPT if the writer sets a route-dependent `label` tag instead.');
  });

  test('every presentation column has a header tooltip (no DG stats-tooltip fallback)', async () => {
    const {noTooltip} = missingPresentationEntries();
    expect(noTooltip.length, 0,
      `columns missing a buildHeaderTooltips entry: [${noTooltip.join(', ')}] — ` +
      'onCellTooltip returns false for these, so DG falls through to its generic ' +
      'min/avg/max/std tooltip, which names no parameter and explains nothing.');
  });

  test('every numeric presentation column has an explicit number format', async () => {
    const {noFormat} = missingPresentationEntries();
    expect(noFormat.length, 0,
      `columns missing a NUMBER_FORMATS entry: [${noFormat.join(', ')}] — ` +
      'each renders at DG default precision beside 3-decimal siblings (this was ' +
      'the third symptom of BUG-12). Add the format, or add it to FORMAT_EXEMPT ' +
      'if the column is non-numeric or deliberately unformatted (TMAX / ' +
      'LAMBDAZ_NPOINTS). Note what "unformatted" actually buys: DG defaults a ' +
      'Float64 column to TEN decimals, which is how t_lag shipped as ' +
      '0.0000000000 — exemption is a real cost, not a neutral omission.');
  });

  test('every COLOURED column explains its colour in the header tooltip', async () => {
    const {noLegend} = missingPresentationEntries();
    expect(noLegend.length, 0,
      `coloured columns whose header tooltip carries no legend: [${noLegend.join(', ')}] — ` +
      'the scientist sees a tinted cell and hovering never says what the tint means. ' +
      'This is how MRT and Vss shipped (user QA 2026-07-30): both are painted by the ' +
      'cross-column fragility hook in configureResultsGrid, and both had a plain ' +
      'tooltip. Attach the band via withLegend() in buildHeaderTooltips.');
  });
});

category('Views: colour bands (one source for the rule and the legend)', () => {
  /** Every number a band's text mentions, in the order it appears. */
  const nums = (s: string): string[] => s.match(/[0-9]+(?:\.[0-9]+)?/g) ?? [];

  test('every band paints and explains the SAME thresholds', async () => {
    // The invariant the single source exists to hold: each threshold a band
    // asserts must appear in BOTH projections — the DG rule that paints the cell
    // and the legend that explains it. Before this, the two were separate string
    // literals that merely happened to agree.
    const bad: string[] = [];
    for (const {column, ruleKey, legend, thresholds} of colourBandProjections()) {
      if (ruleKey === null || legend === null) continue; // cross-column / hidden
      const inRule = nums(ruleKey);
      const inLegend = nums(legend);
      for (const t of thresholds) {
        const s = String(t);
        if (!inRule.includes(s))
          bad.push(`${column}: threshold ${s} missing from rule "${ruleKey}"`);
        if (!inLegend.includes(s))
          bad.push(`${column}: rule "${ruleKey}" paints at ${s}, legend "${legend}" never says so`);
      }
    }
    expect(bad.length, 0, bad.join(' | '));
  });

  test('the legend never states a number the rule does not paint at', async () => {
    // The converse direction. `AUC_DF` reads "≥ 3 — adequate" while DG needs a
    // bounded range (`3-100`) — the 100 is a rendering artifact and is correctly
    // absent from the legend, where it would read as a real ceiling on the
    // Satterthwaite df. This test pins the direction that artifact could break.
    for (const {column, ruleKey, legend} of colourBandProjections()) {
      if (ruleKey === null || legend === null) continue;
      const inRule = nums(ruleKey);
      for (const n of nums(legend)) {
        expect(inRule.includes(n), true,
          `${column}: legend "${legend}" states ${n}, but rule "${ruleKey}" does not paint at it`);
      }
    }
  });

  test('all three grids with colour are covered, sparse included', async () => {
    const grids = new Set(colourBandProjections().map((r) => r.grid));
    for (const g of ['results', 'sparseNca', 'sparseComposite'])
      expect(grids.has(g as never), true, `${g} bands must be under the guard`);
  });

  test('float noise never reaches the legend or the DG rule (adj-R² acceptance band edge)', async () => {
    // Measured, not assumed: a threshold reached by arithmetic (the former
    // `minRSquared + 0.05`) was `0.9500000000000001` in IEEE-754 at some inputs, and that
    // value reached BOTH the conditional range key and the tooltip. The adj-R² band edge is
    // now the ACCEPTANCE threshold R_acc (lambdaAcceptanceRSquared), entered directly, but a
    // value pasted / typed with trailing precision must still be cleaned — fmtThreshold does
    // it. Iterating R_acc across the spec range 0.70–0.95 (plus a noisy value) proves the
    // ruleKey and legend never leak 7+ decimal digits.
    for (const lambdaAcceptanceRSquared of [0.70, 0.80, 0.90, 0.95, 0.9500000000000001]) {
      const rows = colourBandProjections({
        ...DEFAULT_RESULTS_THRESHOLDS, lambdaAcceptanceRSquared,
      }).filter((r) => r.column === 'LAMBDAZ_ADJ_R2');
      expect(rows.length > 0, true, 'adj-R² bands present');
      for (const r of rows) {
        for (const text of [r.ruleKey, r.legend]) {
          if (text === null) continue;
          expect(/[0-9]\.[0-9]{7,}/.test(text), false,
            `R_acc=${lambdaAcceptanceRSquared} leaked IEEE-754 noise into "${text}"`);
        }
      }
    }
  });

  test('the cross-column MRT / Vss band is declared and explained', async () => {
    // The reported defect: both are painted from sibling columns, so no DG rule
    // exists to notice them missing. They must still carry a legend.
    for (const column of ['MRT', 'VSS']) {
      const rows = colourBandProjections().filter((r) => r.column === column);
      expect(rows.length, 1, `${column} declares exactly one (fragile) band`);
      expect(rows[0].ruleKey, null, `${column} is painted cross-column, not by a DG rule`);
      const legend = rows[0].legend ?? '';
      for (const token of ['%AUMCextrap', 'adj. R²', 'N points']) {
        expect(legend.includes(token), true,
          `${column} legend must name the sibling condition "${token}" — it reads "${legend}"`);
      }
    }
  });
});

category('Views: aggregated grid presentation (aggDisplayNames)', () => {
  test('label preferred: %CV on a CL/F group keeps the apparent-clearance distinction', async () => {
    const r = aggDisplayNames('%CV', 'CL', '%CV CL/F');
    expect(r.bare, 'CL/F', 'bare recovered from label, not static "CL"');
    expect(r.header, '%CV CL/F', 'header = label (CL/F preserved)');
  });

  test('label preferred: Mean / SD of a CL/F group', async () => {
    expect(aggDisplayNames('Mean', 'CL', 'Mean CL/F').header, 'Mean CL/F', 'Mean CL/F');
    expect(aggDisplayNames('SD', 'CL', 'SD CL/F').header, 'SD CL/F', 'SD CL/F');
  });

  test('no label -> static friendly name (null and undefined both falsy)', async () => {
    const n = aggDisplayNames('Mean', 'CMAX', null);
    expect(n.bare, 'C_max', 'bare = static friendly');
    expect(n.header, 'Mean C_max', 'header composed from static friendly');
    expect(aggDisplayNames('%CV', 'CMAX', undefined).header, '%CV C_max', 'undefined label -> static');
  });

  test('prefix mismatch degrades VISIBLY to static (never silent corruption)', async () => {
    // If the aggregator prefix and statWord ever diverge, the strip fails and the
    // header falls back to the static param name — a visible "%CV CL" mislabel a
    // reviewer/test catches, not a silent CL/F->CL data confusion.
    const r = aggDisplayNames('%CV', 'CL', 'CV% CL/F'); // legacy / mismatched prefix
    expect(r.bare, 'CL', 'mismatch -> static friendly');
    expect(r.header, '%CV CL', 'mismatch -> composed from static');
  });

  test('unknown param falls back to its raw name', async () => {
    const r = aggDisplayNames('Mean', 'FOO', null);
    expect(r.bare, 'FOO', 'unknown param -> raw name');
    expect(r.header, 'Mean FOO', 'header from raw name');
  });

  test('end-to-end coupling: helper recovers CL/F from a real aggregator label', async () => {
    // Build a real aggregate so the label is produced by the aggregator, not a
    // literal — if the aggregator prefix and the "%CV" statWord ever drift apart,
    // the startsWith strip fails and bare degrades to "CL", failing this test.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'COMPOUND', ['A', 'A', 'B']),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 1.5, 0.8])),
    ]);
    df.col('CL')!.setTag('label', 'CL/F');
    const agg = aggregateNcaResults(df, 'COMPOUND', {params: ['CL']});

    const cvLabel = agg.col('CV_PCT_CL')!.tags['label'];
    expect(cvLabel, '%CV CL/F', 'aggregator writes the %CV prefix');
    expect(aggDisplayNames('%CV', 'CL', cvLabel).header, '%CV CL/F',
      'helper recovers CL/F from the real aggregator label (prefix in lockstep)');
  });
});

// =====================================================================
// GROK-20960 AC-N.3.4 — the C0_EXTRAP CELL tooltip's back-extrapolation line:
// ONE extra line, read from PROVENANCE_JSON.c0,
// shown only when the dose-time value was ESTIMATED. Asserted on the REAL
// tooltip element text over real pipeline output (indometh-1 = 20.6 % of
// AUC₀–∞), and asserted ABSENT on an observed row and on a partial row.
// =====================================================================
category('Views: C0_EXTRAP cell tooltip — back-extrapolated share (AC-N.3.4)', () => {
  const header = () => buildHeaderTooltips(DEFAULT_RESULTS_THRESHOLDS).C0_EXTRAP;
  const textOf = (spec: ReturnType<typeof buildC0TooltipSpec>): string =>
    spec === null ? '' : richTipEl(spec).textContent ?? '';

  test('indometh-1 (logslope c0): the line reads "back-extrapolated: 20.6 % of AUC₀–∞"', async () => {
    const df = await loadFixture('indometh.csv');
    const {results} = runPipelineOrThrow(df, indomethContract(df), indomethScalars(), DEFAULT_RULES);
    let row = -1;
    for (let i = 0; i < results.rowCount; i++)
      if (String(results.col('SUBJECT')!.get(i)) === '1') {row = i; break;}
    expect(row >= 0, true, 'indometh subject 1 present');
    const spec = buildC0TooltipSpec(results, row, header());
    expect(spec !== null, true, 'an estimated c0 row gets the cell tooltip');
    const text = textOf(spec);
    expect(text.includes('back-extrapolated: 20.6 % of AUC₀–∞'), true, `tooltip text was: ${text}`);
    // Exactly ONE extra line — the header spec plus the share, nothing else.
    expect(spec!.extra!.length, 1, 'one line, no column, no band');
    expect(text.startsWith('C_0'), true, 'the header explanation is kept above the line');
    // Every indometh subject is estimated → every row carries the line.
    for (let i = 0; i < results.rowCount; i++)
      expect(buildC0TooltipSpec(results, i, header()) !== null, true, `indometh row ${i} has the line`);
  }, {timeout: 30000});

  test('observed dose-time value (05 multicompound): NO line — nothing was extrapolated', async () => {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const {results} = runPipelineOrThrow(df, multicompoundContract(df), multicompoundScalars(), DEFAULT_RULES);
    for (let i = 0; i < results.rowCount; i++)
      expect(buildC0TooltipSpec(results, i, header()), null, `05 row ${i}: observed c0 → DG default tooltip`);
  }, {timeout: 30000});

  test('partial IV-bolus row (no AUC₀–∞): NO line — there is no share to report', async () => {
    // Two post-dose points: c0 inserts (logslope) but λz cannot fit (minPoints 3).
    const df = makeSyntheticObs({
      subjects: ['P', 'P'], times: Float64Array.from([0.5, 1]), concs: Float64Array.from([4, 3]),
    });
    const {results, outcome} = runPipelineOrThrow(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    expect(results.col('STATUS')!.get(0), 'partial', 'the row is partial');
    expect(outcome.perProfile[0].result.provenance.c0?.method, 'logslope', 'c0 was estimated…');
    expect(buildC0TooltipSpec(results, 0, header()), null, '…but with no AUC₀–∞ the line is omitted');
    // And the JSON carries null for the share (not NaN, which JSON cannot hold).
    const json = JSON.parse(results.col('PROVENANCE_JSON')!.get(0) as string) as
      {c0: {pctAucBackExtrap: number | null}};
    expect(json.c0.pctAucBackExtrap, null, 'PROVENANCE_JSON.c0.pctAucBackExtrap is null on partial');
  }, {timeout: 30000});
});

// nca-cl-bw-normalised-units AD-2: a profile whose body weight varies under a per-kg dose carries
// BW_VARIES_WITHIN_PROFILE; the note is shown on the absolute companions + BW_KG it qualifies,
// and on no other cell.
category('Views: companion cell tooltip — body weight varies (AD-2)', () => {
  const NOTE = 'Body weight varies within the profile (0.2–0.3 kg). CL, Vz and Vss per kg are exact ' +
    '(the dose is per kg); the absolute values use the first row\'s weight (0.2 kg).';
  const frameWith = (warnings: string): DG.DataFrame => DG.DataFrame.fromColumns([
    DG.Column.fromStrings('SUBJECT', ['A', 'B']),
    DG.Column.fromFloat64Array('CL_ABS', Float64Array.from([0.06, 0.07])),
    DG.Column.fromStrings('WARNINGS_JSON', [warnings, '[]']),
  ]);
  const tips = () => buildHeaderTooltips(DEFAULT_RESULTS_THRESHOLDS);

  test('the note is appended under the header on CL_ABS / VZ_ABS / VSS_ABS / BW_KG', async () => {
    const df = frameWith(JSON.stringify([{code: 'BW_VARIES_WITHIN_PROFILE', severity: 'info', message: NOTE}]));
    for (const c of ['CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG']) {
      const spec = buildBwVariesTooltipSpec(df, c, 0, tips()[c]);
      expect(spec !== null, true, `${c}: the note is shown`);
      expect(spec!.extra!.length, 1, `${c}: one line`);
      expect((richTipEl(spec!).textContent ?? '').includes(NOTE), true, `${c}: the profile's own note`);
    }
  });

  test('no note on a profile without it, nor on a per-kg or unrelated column', async () => {
    const df = frameWith(JSON.stringify([{code: 'BW_VARIES_WITHIN_PROFILE', severity: 'info', message: NOTE}]));
    expect(buildBwVariesTooltipSpec(df, 'CL_ABS', 1, tips().CL_ABS), null, 'row without the note → DG default');
    expect(buildBwVariesTooltipSpec(df, 'CL', 0, tips().CL), null, 'per-kg CL is exact → no note');
    const other = frameWith(JSON.stringify([{code: 'PER_KG_NO_BW', severity: 'info', message: 'x'}]));
    expect(buildBwVariesTooltipSpec(other, 'CL_ABS', 0, tips().CL_ABS), null, 'another code → no note');
    expect(buildBwVariesTooltipSpec(frameWith('not json'), 'CL_ABS', 0, tips().CL_ABS), null, 'bad JSON → no note');
  });
});
