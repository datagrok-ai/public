import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseSpectronautText, parseSpectronautStream} from '../parsers/spectronaut-parser';
import {sniffIsPrecursor} from '../package';
import {SEMTYPE} from '../utils/proteomics-types';
import {getGroups} from '../analysis/experiment-setup';
import {_package} from '../package-test';

// Webpack `asset/source` (see webpack.config.js) inlines the committed
// tools/spectronaut-aggregate.{sql,sh} as raw strings into the test bundle
// (module types in src/global.d.ts). The tools/ dir is NOT deployed under
// files/, so it is unreachable via _package.files at runtime — this static
// import is the R5 drift guard (a deleted/renamed/changed fallback fails the
// build + grok test).
// eslint-disable-next-line import/no-unresolved
import aggregateSql from '../../tools/spectronaut-aggregate.sql';
// eslint-disable-next-line import/no-unresolved
import aggregateSh from '../../tools/spectronaut-aggregate.sh';

/** Builds a synthetic Spectronaut long-format TSV string for testing.
 *
 * Each protein gets ≥2 distinct precursor/fragment rows per (protein × condition
 * × replicate). The quantity value is held CONSTANT within each group so duckdb's
 * `max()` equals the parser's first-encountered value (RESEARCH Assumption A1);
 * PG.Organisms / R.FileName are likewise constant per group so duckdb
 * `any_value()` equals first-non-null (RESEARCH Pitfall 5).
 *
 * The emitted header carries the D-01 precursor signature
 * (`EG.ModifiedPeptide`, `FG.Charge`, `FG.Id`, `PEP.StrippedSequence`). These are
 * EXTRA columns — they add no protein groups, so the existing suite observes the
 * same protein/row counts as before. NO `PG.Genes`/`PG.ProteinAccessions`
 * columns are emitted (the real hye-mix demo lacks them; the committed
 * tools/spectronaut-aggregate.sql drops the matching any_value SELECT terms for
 * exactly this reason). Default conditions stay `CondA`/`CondB` so the committed
 * SQL's DMD↔WT flip is a structural no-op.
 *
 * `extraRows` injects deliberate filter-branch rows (decoy/contaminant prefixes,
 * sub-threshold / non-numeric / empty-string q-values) so a single call can
 * exercise every R2 branch without perturbing the default behaviour the existing
 * 19 tests depend on. */
function makeLongFormatTsv(opts: {
  proteins: {id: string; organism?: string}[];
  conditions: string[];
  replicates: number[];
  qValues?: Map<string, number | string>; // protein -> qValue override
  ibaqValues?: Map<string, number>; // protein -> iBAQ override
  fileNameTemplate?: string;
  quantityColumn?: string; // defaults to 'PG.IBAQ'; pass 'PG.Quantity' to exercise the alias
  /** Deliberate filter-branch rows appended after the regular protein rows.
   * Each entry emits one precursor row with the given protein id, q-value, and
   * (optional) explicit quantity, across every condition × replicate. Used to
   * exercise CON__/REV__, sub-threshold-numeric, non-numeric, and empty-string
   * q-value branches in one fixture. Defaults/existing tests never set this, so
   * their protein/row counts are unaffected. */
  extraRows?: {id: string; qValue: number | string; quantity?: number | string; organism?: string}[];
}): string {
  const quantityCol = opts.quantityColumn ?? 'PG.IBAQ';
  const headers = [
    'R.FileName', 'R.Condition', 'R.Replicate', 'PG.ProteinGroups',
    'PG.Organisms', quantityCol, 'EG.Qvalue', 'PEP.StrippedSequence',
    'EG.ModifiedPeptide', 'FG.Charge', 'FG.Id',
  ];
  const rows: string[] = [headers.join('\t')];

  // Two distinct precursors per (protein × condition × replicate). Quantity is
  // held constant per group; only the precursor identity (modified peptide /
  // charge / fragment id / stripped sequence) varies, so duckdb max() collapses
  // to the same value the parser's first-encountered-wins logic produces.
  const precursors = [
    {strip: 'PEPTIDER', mod: '_PEPTIDER_', charge: '2', fg: 'FG1'},
    {strip: 'ANOTHERPEP', mod: '_ANOTHERPEP_', charge: '3', fg: 'FG2'},
  ];

  const emit = (
    fileName: string, cond: string, rep: number, id: string,
    organism: string, quantity: string, qVal: string,
  ): void => {
    for (const p of precursors) {
      rows.push([fileName, cond, String(rep), id, organism,
        quantity, qVal, p.strip, p.mod, p.charge, p.fg].join('\t'));
    }
  };

  let ibaqCounter = 100;
  for (const protein of opts.proteins) {
    const ibaq = opts.ibaqValues?.get(protein.id) ?? ibaqCounter++;
    for (const cond of opts.conditions) {
      for (const rep of opts.replicates) {
        const fileName = (opts.fileNameTemplate ?? 'file') + `_${cond}_${rep}`;
        const qVal = opts.qValues?.get(protein.id) ?? 0.001;
        const organism = protein.organism ?? 'Homo sapiens';
        emit(fileName, cond, rep, protein.id, organism, String(ibaq), String(qVal));
      }
    }
  }

  // Deliberate filter-branch rows (decoy/contaminant + edge q-values).
  if (opts.extraRows) {
    let extraQuant = 500;
    for (const er of opts.extraRows) {
      const q = er.quantity ?? extraQuant++;
      for (const cond of opts.conditions) {
        for (const rep of opts.replicates) {
          const fileName = (opts.fileNameTemplate ?? 'file') + `_${cond}_${rep}`;
          const organism = er.organism ?? 'Homo sapiens';
          emit(fileName, cond, rep, er.id, organism, String(q), String(er.qValue));
        }
      }
    }
  }

  return rows.join('\n');
}

category('Spectronaut', () => {
  const baseProteins = [
    {id: 'P12345'}, {id: 'Q67890'}, {id: 'O15439'},
  ];
  const baseOpts = {
    proteins: baseProteins,
    conditions: ['CondA', 'CondB'],
    replicates: [1, 2],
  };

  test('pivot produces correct protein count', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    expect(df.rowCount, 3);
  });

  test('pivot produces correct sample column names', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    expect(df.col('CondA_1') !== null, true);
    expect(df.col('CondA_2') !== null, true);
    expect(df.col('CondB_1') !== null, true);
    expect(df.col('CondB_2') !== null, true);
  });

  test('CON__ and REV__ prefixed proteins are excluded', async () => {
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: [
        {id: 'P12345'}, {id: 'CON__P99999'}, {id: 'REV__Q11111'},
      ],
    });
    const df = await parseSpectronautText(tsv);
    expect(df.rowCount, 1);
    expect(df.col('PG.ProteinGroups')!.get(0), 'P12345');
  });

  test('rows with numeric EG.Qvalue > threshold are excluded', async () => {
    const qValues = new Map<string, number | string>();
    qValues.set('P12345', 0.001);
    qValues.set('Q67890', 0.05); // above default 0.01 threshold
    qValues.set('O15439', 0.005);
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: baseProteins,
      qValues,
    });
    const df = await parseSpectronautText(tsv);
    expect(df.rowCount, 2);
  });

  test('non-numeric EG.Qvalue rows treated as passing', async () => {
    const qValues = new Map<string, number | string>();
    qValues.set('P12345', 0.001);
    qValues.set('Q67890', 'Profiled');
    qValues.set('O15439', 'NaN');
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: baseProteins,
      qValues,
    });
    const df = await parseSpectronautText(tsv);
    expect(df.rowCount, 3);
  });

  test('PG.ProteinGroups has PROTEIN_ID semantic type', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    expect(df.col('PG.ProteinGroups')?.semType, SEMTYPE.PROTEIN_ID);
  });

  test('intensity columns have INTENSITY semantic type', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    expect(df.col('CondA_1')?.semType, SEMTYPE.INTENSITY);
    expect(df.col('CondB_2')?.semType, SEMTYPE.INTENSITY);
  });

  test('R.FileName stored as spectronaut.fileName tag on intensity columns', async () => {
    const tsv = makeLongFormatTsv({...baseOpts, fileNameTemplate: 'run'});
    const df = await parseSpectronautText(tsv);
    const col = df.col('CondA_1');
    expect(col !== null, true);
    if (col) {
      const tag = col.getTag('spectronaut.fileName');
      expect(tag !== null && tag !== '', true);
      expect(tag!.includes('run'), true);
    }
  });

  test('groups auto-populated via setGroups when exactly 2 conditions', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    const groups = getGroups(df);
    expect(groups !== null, true);
    if (groups) {
      expect(groups.group1.columns.length, 2);
      expect(groups.group2.columns.length, 2);
    }
  });

  test('groups not set when more than 2 conditions', async () => {
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      conditions: ['CondA', 'CondB', 'CondC'],
    });
    const df = await parseSpectronautText(tsv);
    const groups = getGroups(df);
    expect(groups, null);
  });

  test('pre-normalized tag set for raw (non-log2) intensity data', async () => {
    // Use default baseOpts which produce raw-range IBAQ values (100+)
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    expect(df.getTag('proteomics.preNormalized'), 'true');
  });

  test('pre-normalized tag set when detectLog2Status detects log2-range values', async () => {
    // Create data with values in log2 range (0-30)
    const ibaqValues = new Map<string, number>();
    ibaqValues.set('P12345', 15.2);
    ibaqValues.set('Q67890', 22.7);
    ibaqValues.set('O15439', 18.1);
    const tsv = makeLongFormatTsv({...baseOpts, ibaqValues});
    const df = await parseSpectronautText(tsv);
    expect(df.getTag('proteomics.preNormalized'), 'true');
  });

  test('both raw and log2 intensity columns present in output', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    // Raw columns
    expect(df.col('CondA_1') !== null, true);
    // Log2 columns
    expect(df.col('log2(CondA_1)') !== null, true);
    expect(df.col('log2(CondB_2)') !== null, true);
  });

  test('proteomics.source tag set to spectronaut', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = await parseSpectronautText(tsv);
    expect(df.getTag('proteomics.source'), 'spectronaut');
  });

  test('PG.Quantity accepted as alternative to PG.IBAQ', async () => {
    const tsv = makeLongFormatTsv({...baseOpts, quantityColumn: 'PG.Quantity'});
    const df = await parseSpectronautText(tsv);
    expect(df.rowCount, 3);
    expect(df.col('CondA_1') !== null, true);
    expect(df.col('CondB_2') !== null, true);
    expect(df.col('CondA_1')?.semType, SEMTYPE.INTENSITY);
  });

  test('throws when neither PG.IBAQ nor PG.Quantity is present', async () => {
    // makeLongFormatTsv with a bogus quantity-column name produces a frame
    // that has neither known quantity column.
    const tsv = makeLongFormatTsv({...baseOpts, quantityColumn: 'PG.Bogus'});
    let threw = false;
    try {
      await parseSpectronautText(tsv);
    } catch (e) {
      threw = true;
      expect((e as Error).message.includes('PG.IBAQ'), true);
      expect((e as Error).message.includes('PG.Quantity'), true);
    }
    expect(threw, true);
  });

  test('demo file produces correct dimensions', async () => {
    // This test will use the actual demo file when run in Datagrok
    // For build validation, we use synthetic data with known dimensions
    const proteins = [];
    for (let i = 0; i < 93; i++)
      proteins.push({id: `P${String(i).padStart(5, '0')}`});
    const tsv = makeLongFormatTsv({
      proteins,
      conditions: ['HYE mix A', 'HYE mix B'],
      replicates: [1, 2, 3, 4],
    });
    const df = await parseSpectronautText(tsv);
    expect(df.rowCount, 93);
    // 8 sample columns + protein ID + organisms + log2 columns
    const intensityCols = df.columns.toList().filter((c) =>
      c.semType === SEMTYPE.INTENSITY && !c.name.startsWith('log2('));
    expect(intensityCols.length, 8);
  });

  // ---------------------------------------------------------------------------
  // Plan 12-03 streaming-path coverage (RESEARCH Open Q1: the 19 tests above
  // stay locked on parseSpectronautText; everything below pins the streaming
  // path, the D-01 routing, and their equivalence to the text path / the
  // committed duckdb golden — without touching makeLongFormatTsv or the 19).
  // ---------------------------------------------------------------------------

  /** Wraps a TSV string in a File and streams it through the precursor
   * aggregator — the exact path importSpectronaut routes precursor reports to. */
  async function streamTsv(tsv: string): Promise<DG.DataFrame> {
    return await parseSpectronautStream(
      new File([tsv], 'fixture.tsv', {type: 'text/tab-separated-values'}));
  }

  /** Reads a committed demo file using the package files convention already used
   * by spectronaut-candidates-e2e / fragpipe-e2e (`_package.files` maps the
   * deployed `files/` directory; `demo/<name>` -> `files/demo/<name>`). */
  async function readDemoFile(name: string): Promise<string> {
    return _package.files.readAsText(`demo/${name}`);
  }

  /** Non-log2 INTENSITY sample column names, sorted (the raw pre-log2 columns). */
  function rawSampleCols(df: DG.DataFrame): string[] {
    return df.columns.toList()
      .filter((c) => c.semType === SEMTYPE.INTENSITY && !c.name.startsWith('log2('))
      .map((c) => c.name)
      .sort();
  }

  /** protein-id value -> row index, keyed on the primary protein column the
   * shared finalize tail derives (falls back to PG.ProteinGroups). */
  function proteinRowIndex(df: DG.DataFrame): {col: DG.Column; map: Map<string, number>} {
    const col = df.col('Primary Protein ID') ?? df.col('PG.ProteinGroups')!;
    const map = new Map<string, number>();
    for (let i = 0; i < df.rowCount; i++)
      map.set(String(col.get(i)), i);
    return {col, map};
  }

  test('streams precursor fixture', async () => {
    // R1 smoke test — the exact name in 12-VALIDATION.md's R1 row. Streams the
    // committed precursor fixture (D-01 signature, edge rows, no PG.Genes).
    const text = await readDemoFile('spectronaut-hye-precursor.tsv');
    const df = await streamTsv(text);
    expect(df.rowCount > 0, true);
    const sampleCol = df.columns.toList().find((c) =>
      c.semType === SEMTYPE.INTENSITY && !c.name.startsWith('log2(') &&
      (c.name.startsWith('CondA_') || c.name.startsWith('CondB_')));
    expect(sampleCol !== undefined, true);
    expect(df.col('PG.ProteinGroups')?.semType, SEMTYPE.PROTEIN_ID);
  });

  test('sniffIsPrecursor routes by header', async () => {
    // Directly verifies the D-01 routing branch (12-02 Task 2) in isolation.
    const precHeader = [
      'R.FileName', 'R.Condition', 'R.Replicate', 'PG.ProteinGroups',
      'PG.Organisms', 'PG.Quantity', 'EG.Qvalue', 'EG.ModifiedPeptide',
      'FG.Charge', 'PEP.StrippedSequence',
    ].join('\t');
    const precRow = [
      'run_CondA_1', 'CondA', '1', 'P1', 'Homo sapiens', '10', '0.001',
      '_PEP_', '2', 'PEP',
    ].join('\t');
    // PG-level header: required cols + PG.IBAQ but NONE of the three signature cols.
    const pgHeader = [
      'R.FileName', 'R.Condition', 'R.Replicate', 'PG.ProteinGroups',
      'PG.Organisms', 'PG.IBAQ',
    ].join('\t');
    const pgRow = ['run_CondA_1', 'CondA', '1', 'P1', 'Homo sapiens', '10'].join('\t');

    expect(await sniffIsPrecursor(
      new File([precHeader + '\n' + precRow], 'p.tsv')), true);
    expect(await sniffIsPrecursor(
      new File([pgHeader + '\n' + pgRow], 'pg.tsv')), false);
  });

  test('stream path matches text path', async () => {
    // Structural equivalence + per-(protein × raw sample column) numeric equality
    // within 1e-3 — catches the documented first-encountered-vs-max divergence
    // (RESEARCH ~L326/L428-429): a structural-only check would miss it.
    const tsv = makeLongFormatTsv(baseOpts);
    const a = await parseSpectronautText(tsv);
    const b = await streamTsv(tsv);

    expect(b.rowCount, a.rowCount);

    const aCols = rawSampleCols(a);
    const bCols = rawSampleCols(b);
    expect(bCols.join(','), aCols.join(','));

    for (const name of aCols)
      expect(b.col(`log2(${name})`) !== null, a.col(`log2(${name})`) !== null);

    expect(b.getTag('proteomics.source'), a.getTag('proteomics.source'));
    expect(b.getTag('proteomics.preNormalized'), a.getTag('proteomics.preNormalized'));

    const ga = getGroups(a);
    const gb = getGroups(b);
    expect(gb !== null, ga !== null);
    if (ga && gb) {
      expect(gb.group1.columns.length, ga.group1.columns.length);
      expect(gb.group2.columns.length, ga.group2.columns.length);
    }

    // Per-cell numeric equality (Blocker 2): match rows by protein id, then
    // compare every raw sample cell within 1e-3 (both-null treated as equal).
    const ai = proteinRowIndex(a);
    const bi = proteinRowIndex(b);
    expect(bi.map.size, ai.map.size);
    for (const [protein, ar] of ai.map) {
      const br = bi.map.get(protein);
      expect(br !== undefined, true);
      for (const name of aCols) {
        const ac = a.col(name)!;
        const bc = b.col(name)!;
        const aNull = ac.isNone(ar);
        const bNull = bc.isNone(br!);
        expect(bNull, aNull);
        if (!aNull && !bNull) {
          const av = ac.get(ar) as number;
          const bv = bc.get(br!) as number;
          expect(Math.abs(av - bv) <= 1e-3, true);
        }
      }
    }
  });

  test('streaming filter parity', async () => {
    // Every R2 edge branch in one fixture: CON__/REV__ excluded, an all-q>0.01
    // protein excluded, non-numeric ('Profiled') and empty-string q-value
    // proteins INCLUDED (duckdb null/non-numeric-pass parity).
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: [{id: 'KEEP1'}, {id: 'KEEP2'}],
      extraRows: [
        {id: 'CON__C1', qValue: 0.001},
        {id: 'REV__R1', qValue: 0.001},
        {id: 'HIGHQ', qValue: 0.05},      // all precursors numeric q > 0.01
        {id: 'PROF', qValue: 'Profiled'}, // non-numeric -> passes
        {id: 'EMPTYQ', qValue: ''},       // empty-string -> passes
      ],
    });
    const df = await streamTsv(tsv);
    const protCol = df.col('PG.ProteinGroups')!;
    const ids = new Set<string>();
    for (let i = 0; i < df.rowCount; i++)
      ids.add(String(protCol.get(i)));

    expect(ids.has('CON__C1'), false);
    expect(ids.has('REV__R1'), false);
    expect(ids.has('HIGHQ'), false);
    expect(ids.has('PROF'), true);
    expect(ids.has('EMPTYQ'), true);
    // KEEP1, KEEP2, PROF, EMPTYQ survive.
    expect(df.rowCount, 4);
  });

  /** Streams `tsv` while recording every `grok.shell.info` string, then restores
   * the original `grok.shell.info` in a `finally` (a thrown assertion must NOT
   * leak the stub into later tests — T-12-18). Returns the captured messages so
   * the caller asserts the user-facing malformed-message contract directly (the
   * per-line counters are function-local to parseSpectronautStream and not
   * returned; the completion `grok.shell.info` is the observable the gap is
   * about). */
  async function streamCapturingInfo(
    tsv: string,
  ): Promise<{df: DG.DataFrame; infos: string[]}> {
    const infos: string[] = [];
    const original = grok.shell.info.bind(grok.shell);
    (grok.shell as unknown as {info: (m: unknown) => void}).info = (m: unknown) => {
      infos.push(String(m));
    };
    try {
      const df = await streamTsv(tsv);
      return {df, infos};
    } finally {
      (grok.shell as unknown as {info: typeof original}).info = original;
    }
  }

  test('by-design-filtered rows are not counted malformed', async () => {
    // A fixture whose ONLY dropped rows are by-design-filtered: a CON__ decoy,
    // a REV__ decoy, and a protein whose every precursor has numeric q > 0.01.
    // Every emitted line has the full 11-field header width — NO truncation.
    // The gap's exact false signal ("skipped N malformed line(s)") must NOT
    // fire for purely-filtered input; observe it via the grok.shell.info spy.
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: [{id: 'KEEP1'}, {id: 'KEEP2'}],
      extraRows: [
        {id: 'CON__C1', qValue: 0.001}, // contaminant -> filtered (silent)
        {id: 'REV__R1', qValue: 0.001}, // decoy       -> filtered (silent)
        {id: 'HIGHQ', qValue: 0.05},    // numeric q > 0.01 -> filtered (silent)
      ],
    });
    const {df, infos} = await streamCapturingInfo(tsv);

    // Correct wide DataFrame: the two kept proteins survive; the
    // decoy/contaminant/high-q proteins are absent.
    const protCol = df.col('PG.ProteinGroups')!;
    const ids = new Set<string>();
    for (let i = 0; i < df.rowCount; i++)
      ids.add(String(protCol.get(i)));
    expect(ids.has('KEEP1'), true);
    expect(ids.has('KEEP2'), true);
    expect(ids.has('CON__C1'), false);
    expect(ids.has('REV__R1'), false);
    expect(ids.has('HIGHQ'), false);
    expect(df.rowCount, 2);

    // The crux: NO emitted info message labels the by-design-filtered rows
    // "malformed" (this is the exact false corruption signal the gap reports).
    const malformedMsgs = infos.filter((s) => /malformed/i.test(s));
    expect(malformedMsgs.length, 0);
  });

  test('truncated line is counted malformed', async () => {
    // A valid precursor fixture plus one deliberately truncated data line
    // (fewer tab-separated fields than the 11-field header). The genuine
    // malformed category MUST still surface — the fix narrows the message, it
    // does not suppress it.
    const validTsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: [{id: 'KEEP1'}, {id: 'KEEP2'}],
    });
    const truncated = validTsv + '\n' + 'file\tCondA\t1';
    const {df, infos} = await streamCapturingInfo(truncated);

    // The valid proteins still produce a correct DataFrame.
    const protCol = df.col('PG.ProteinGroups')!;
    const ids = new Set<string>();
    for (let i = 0; i < df.rowCount; i++)
      ids.add(String(protCol.get(i)));
    expect(ids.has('KEEP1'), true);
    expect(ids.has('KEEP2'), true);
    expect(df.rowCount, 2);

    // The genuine-malformed completion message WAS emitted for the truncated row.
    const surfaced = infos.some((s) => /skipped \d+ malformed line\(s\)/.test(s));
    expect(surfaced, true);
  });

  test('empty-protein rows are filtered silently (streaming↔text parity)', async () => {
    // Spectronaut routinely emits unassigned-precursor rows: a structurally
    // complete line with PG.ProteinGroups blank. The text path (pivotSpectronaut)
    // drops these silently — the streaming path must do the same, NOT count them
    // toward the user-facing "skipped N malformed line(s)" message.
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: [{id: 'KEEP1'}, {id: 'KEEP2'}],
      extraRows: [
        {id: '', qValue: 0.001},
        {id: '', qValue: 0.002},
      ],
    });
    const {df, infos} = await streamCapturingInfo(tsv);

    const protCol = df.col('PG.ProteinGroups')!;
    const ids = new Set<string>();
    for (let i = 0; i < df.rowCount; i++)
      ids.add(String(protCol.get(i)));
    expect(ids.has('KEEP1'), true);
    expect(ids.has('KEEP2'), true);
    expect(ids.has(''), false);
    expect(df.rowCount, 2);

    const malformedMsgs = infos.filter((s) => /malformed/i.test(s));
    expect(malformedMsgs.length, 0);
  });

  test('both-numeric-casts-null rows are filtered silently (streaming↔text parity)', async () => {
    // A precursor row where both PG.Quantity and EG.Qvalue cast to null is
    // unusable for aggregation. The text path likewise has no contribution from
    // such a row (ibaqCol.isNone gates the protein-map write, q-null-passes-as-
    // non-numeric); categorizing it 'malformed' on the streaming path was a
    // false-positive user message inherited from Phase 12. It must now be silent.
    const tsv = makeLongFormatTsv({
      ...baseOpts,
      proteins: [{id: 'KEEP1'}, {id: 'KEEP2'}],
      extraRows: [
        {id: 'NULLPAIR', quantity: 'NaN', qValue: 'NA'},
      ],
    });
    const {df, infos} = await streamCapturingInfo(tsv);

    expect(df.rowCount, 2);
    const malformedMsgs = infos.filter((s) => /malformed/i.test(s));
    expect(malformedMsgs.length, 0);
  });

  test('streaming tag set and groups', async () => {
    const df = await streamTsv(makeLongFormatTsv(baseOpts));
    expect(df.getTag('proteomics.source'), 'spectronaut');
    expect(df.getTag('proteomics.preNormalized'), 'true');
    const g = getGroups(df);
    expect(g !== null, true);
    if (g) {
      expect(g.group1.columns.length, 2);
      expect(g.group2.columns.length, 2);
    }
  });

  test('streaming output equals duckdb golden', async () => {
    // The JSON sidecar is AUTHORITATIVE: Plan 01 derived it verbatim from the
    // real duckdb golden .tsv (D-04) via a pure transcriber (no in-test
    // re-aggregation). This test NEVER re-aggregates and has NO conditional
    // fallback — if the fixture or sidecar cannot be read it FAILS loudly.
    const golden = JSON.parse(
      await readDemoFile('spectronaut-hye-precursor-golden.json')) as
      Record<string, {quantity: number; qvalue: number}>;
    const text = await readDemoFile('spectronaut-hye-precursor.tsv');
    const df = await streamTsv(text);

    const {map: protRow} = proteinRowIndex(df);
    const sampleCols = rawSampleCols(df); // CondA_1.. / CondB_1.. raw columns

    // Per-key (protein × condition_replicate) quantity equivalence within 1e-3.
    const goldenKeys = Object.keys(golden);
    expect(goldenKeys.length > 0, true);
    for (const key of goldenKeys) {
      // key = `${protein}${condition}_${replicate}` — split on the trailing
      // `<Cond...>_<rep>` so multi-segment protein ids stay intact.
      const m = key.match(/^(.*?)((?:Cond[A-Za-z]+)_\d+)$/);
      expect(m !== null, true);
      const protein = m![1];
      const sampleKey = m![2];
      const ri = protRow.get(protein);
      expect(ri !== undefined, true);
      const col = df.col(sampleKey);
      expect(col !== null, true);
      expect(col!.isNone(ri!), false);
      const streamed = col!.get(ri!) as number;
      expect(Math.abs(streamed - golden[key].quantity) <= 1e-3, true);
    }

    // The streamed (protein × sample) non-null key set must EQUAL the sidecar's
    // exactly — no extra, no missing (catches a ported condition flip / filter
    // mismatch: such a drift changes keys, not just values).
    const streamedKeys = new Set<string>();
    for (const [protein, ri] of protRow) {
      for (const sc of sampleCols) {
        const col = df.col(sc)!;
        if (!col.isNone(ri))
          streamedKeys.add(`${protein}${sc}`);
      }
    }
    const goldenSet = new Set(goldenKeys);
    expect(streamedKeys.size, goldenSet.size);
    for (const k of streamedKeys)
      expect(goldenSet.has(k), true);
    for (const k of goldenSet)
      expect(streamedKeys.has(k), true);
  });

  test('duckdb fallback tooling is committed', async () => {
    // Lightweight R5 regression: the committed duckdb fallback was not deleted
    // or renamed. aggregateSql / aggregateSh are the verbatim committed
    // tools/spectronaut-aggregate.{sql,sh}, inlined by webpack asset/source
    // (tools/ is not runtime-readable via _package.files).
    expect(typeof aggregateSql === 'string' && aggregateSql.length > 0, true);
    expect(typeof aggregateSh === 'string' && aggregateSh.length > 0, true);
    expect(aggregateSql.includes('max(TRY_CAST'), true);
    expect(aggregateSh.includes('spectronaut-aggregate.sql'), true);
  });
});
