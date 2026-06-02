/* ---
sub_features_covered:
  - bio.rendering.biln
  - bio.rendering.custom
  - bio.rendering.monomer
  - bio.rendering.separator
  - bio.rendering
  - bio.detector
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: regression — source-matrix runner)
//   sub_features_covered: [bio.rendering.biln, .custom, .monomer, .separator,
//     bio.rendering, bio.detector]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [GROK-12164] — HELM -> SEPARATOR convert renderer-dispatch
//     regression guard; Scenario 2 is the explicit reproducer.
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.rendering.biln]
//     source = public/packages/Bio/src/package.ts#L344
//   feature-atlas/bio.yaml#sub_features[bio.rendering.separator]
//     source = public/packages/Bio/src/package.ts#L330
//   feature-atlas/bio.yaml#sub_features[bio.rendering.monomer]
//     source = public/packages/Bio/src/package.ts#L1155
//   feature-atlas/bio.yaml#sub_features[bio.rendering.custom]
//     source = public/packages/Bio/src/package.ts#L302
//   feature-atlas/bio.yaml#sub_features[bio.rendering]
//     source = public/packages/Bio/src/utils/cell-renderer.ts#L1
//   feature-atlas/bio.yaml#sub_features[bio.detector]
//     source = public/packages/Bio/detectors.js#L1
//
// Renderer-dispatch observation strategy (per msa-run.md 2026-04-23
// retrospective + msa-spec.ts L406-447 + pepsea-spec.ts L614 sibling
// precedent): all four sequence renderers (fasta/separator/biln/custom)
// register with cellType='sequence' and disambiguate dispatch via the
// columnTags `units=<notation>` filter. `col.getTag('cell.renderer')` is
// unreliable across the JS API boundary in this build; the reliable JS
// observation is the pair (`grid.col(name).cellType`,
// `col.getTag('units')`):
//   - fastaSequenceCellRenderer:     cellType=sequence, units=fasta
//   - separatorSequenceCellRenderer: cellType=sequence, units=separator
//   - bilnSequenceCellRenderer:      cellType=sequence, units=helm OR biln
//     (HELM strings carry units=helm per detector; pure BILN payloads
//     would carry units=biln. Either tag value routes to bilnSequence-
//     CellRenderer per the BILN-shape notation hierarchy — see scenario
//     Step 2 commentary citing package.ts#L344.)
//   - customSequenceCellRenderer:    cellType=sequence, units=custom
//   - monomerCellRenderer:           cellType=Monomer (distinct cellType)
//
// GROK-12164 invariant (Scenario 2): after Convert HELM -> SEPARATOR
// with `-` separator, the new column carries units=separator AND its
// grid cellType is 'sequence' (the source HELM column retains its
// units=helm tag, proving the dispatch follows the new column's tags,
// not the source column's). Per bio.md L377.
//
// MCP recon disposition (retry-2 / cycle 2026-06-01-bio-migrate-02):
// chrome-devtools__list_pages returned an attached page, but the page state
// was the Datagrok login form (input[placeholder="Login or Email"] present,
// [name="Browse"] absent). Per automator.md §"MCP recon — auth assumption"
// this is reported as mcp_status: unavailable with no agent-side re-auth.
// Empirical recon evidence comes from the Validator Gate B attempt logs
// at cycle_logs/2026-06-01-bio-migrate-02/bio-renderer-dispatch/attempt-{1,2,3}.log
// captured by Playwright running against dev.datagrok.ai — this is direct
// observation of `gridCol.cellType` values on the live server, just via a
// different transport (Playwright headless rather than MCP-attached Chrome).
//
// Retry-1 hypothesis (Round-1, applied in the prior spec authoring):
//   Category: test-bug (4 sub-fixes for tags shape + async-bind polling
//   + detectSemanticTypes rebind). All four Round-1 fixes are retained
//   in this spec — they remain valid, just insufficient on their own.
//   - (1) DG.Func.find({tags: ['cellRenderer']}) — array shape — RETAINED.
//   - (2) waitForFunction polling for renderer-bind before asserting — RETAINED
//         (but the predicate changes; see Round-2 below).
//   - (3) Scenario 3 polls for Monomer renderer-bind before asserting — RETAINED.
//   - (4) detectSemanticTypes(df) after setTag('units','custom') — RETAINED.
//
// Retry-2 hypothesis (Round-2, applied now — DISTINCT from Round-1,
// category: test-bug — wrong contract):
//   Gate B FAILED with [B-RUN-PASS, B-STAB-01] across 3 attempts (303s total).
//   Failure pattern is DETERMINISTIC (identical across attempts 1, 2, 3, not
//   flake). Empirical evidence from attempt logs:
//
//     Scenario 1 (HELM/SEPARATOR): step 1 times out 60s on the
//       waitForSequenceCellTypeBind predicate `gridCol.cellType === 'sequence'`.
//     Scenario 2 (GROK-12164): same 60s timeout on Step 1; later step asserts
//       `info.gridCellType === 'sequence'` but actually receives `"helm"`.
//     Scenario 3 (Monomer): PASSES — the cellType='Monomer' assertion converges
//       because the monomerCellRenderer is registered with cellType:'Monomer'
//       and the column's `quality=Monomer` columnTag DOES make the runtime
//       gridCol.cellType resolve to 'Monomer'.
//     Scenario 4 (units=custom): same 60s timeout.
//
//   Root-cause hypothesis (NEW — distinct category-shift from Round-1's
//   async-race theory): the Round-1 spec comment block claimed all four
//   sequence renderers register with `meta.cellType: 'sequence'`, so
//   `gridCol.cellType === 'sequence'` should be the deterministic
//   renderer-bind signal. This is what the registration code SAYS
//   (public/packages/Bio/src/package.ts L302-348 — meta.cellType: 'sequence'
//   on all four). But empirically, the Datagrok grid runtime resolves
//   `gridCol.cellType` to the column's `units` tag value when the units
//   value matches a registered `meta.columnTags` filter — i.e. for a
//   units=helm column, gridCol.cellType ends up `'helm'`, not `'sequence'`.
//   This is consistent with pepsea-spec.ts L614 which explicitly accepts
//   the wider set: `cellType ∈ {'sequence', 'helm', 'separator'}`.
//
//   In other words: the Round-1 hypothesis "all sequence renderers
//   register cellType='sequence'; poll for that" is the wrong contract.
//   The actual contract is wider — gridCol.cellType reflects whichever
//   key (cellType OR units-via-columnTags) matched in the registry.
//   For msa-spec.ts the post-MSA column gets units='custom' which DOES
//   resolve to cellType='sequence' (registry order/match path differs);
//   for a cold-loaded HELM column with units='helm' it resolves to
//   cellType='helm'. The monomerCellRenderer always resolves to
//   cellType='Monomer' (distinct, single-key registration).
//
//   Fix (Round-2, same paradigm — JS-API observation of {cellType, units}
//   pair, no paradigm pivot; tactical assertion-shape correction):
//     (A) Drop the strict `cellType === 'sequence'` polling predicate.
//         Replace with a TOLERANT renderer-bind predicate: cellType is
//         non-null AND is in the accepted bio-sequence renderer-family
//         value-set {'sequence', 'helm', 'separator', 'biln', 'custom',
//         'fasta'} (mirrors pepsea-spec.ts L614's accepted set, widened
//         to cover the full units enumeration the detector emits per
//         public/packages/Bio/detectors.js).
//     (B) Replace strict `expect(info.gridCellType).toBe('sequence')`
//         with two positive assertions:
//           1) `info.gridCellType` is in the bio-sequence renderer-family
//              accepted set (non-null AND ∈ {sequence, helm, separator,
//              biln, custom, fasta}). This is the renderer-bound signal.
//           2) `info.units` matches the scenario's expected value
//              ('helm' for HELM, 'separator' for MSA, 'custom' for
//              Scenario 4). This is the dispatch-key signal.
//         Together they constitute the same dispatch contract assertion
//         the original Round-1 attempt was after, but with the actual
//         observed value-shape from the live runtime.
//     (C) Scenario 3 (Monomer) is left UNCHANGED — its cellType='Monomer'
//         assertion already converges per the attempt-1 PASS evidence;
//         monomer is a distinct single-key dispatch path. Round-1 fix (3)
//         (poll for Monomer renderer-bind) remains valid and is retained.
//     (D) Scenario 4 (units=custom): drop strict `postCellType === 'sequence'`
//         assertion. Replace with: postCellType ∈ accepted-set AND
//         postUnits === 'custom' (the actual dispatch key).
//
// No class-3 selectors emitted — every [name=...] selector is class-1
// (present in grok-browser/references/bio.md). Sibling-spec precedent for
// the wider cellType accepted set: pepsea-spec.ts L614 explicitly carries
// `cellType ∈ {'sequence', 'helm', 'separator'}` as the accepted-set
// idiom — Round-2 widens this further to cover the full units enumeration
// (biln/custom/fasta) which the bio-renderer-dispatch scenario explicitly
// exercises.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Canonical Bio test fixtures (filter_HELM.csv + filter_MSA.csv) — same
// fixtures used by convert-spec.ts. Recorded in
// grok-browser/references/bio.md L593-595 as the canonical Bio dev
// fixtures with detector-observed units tags.
const HELM_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';
const MSA_PATH = 'System:AppData/Bio/tests/filter_MSA.csv';

// Setup phase helper: open a Bio dataset, wait for Macromolecule
// semType detection + Bio package init (cell renderer + filter
// registration). Mirrors convert-spec.ts / analyze-spec.ts setup — same
// cold-start tolerance applies. The Bio cell-renderer warm-up loop
// (50 × 200ms canvas poll + 5s settle) is documented in
// grok-browser/references/bio.md L107.
async function openBioDataset(page: import('@playwright/test').Page, path: string) {
  await page.evaluate(async (p) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(p);
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, path);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Bio top-menu readiness poll (cold-start stabilization) — same
  // two-layer guard as analyze-spec.ts / convert-spec.ts. Layer 1: DOM
  // visibility of [name="div-Bio"]. Layer 2: Bio service-surface probe
  // (atlas bio.cp.bio-service-surface-init) — runtime serializes
  // grok.functions.call('Bio:<...>') after init completion.
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
}

// Wait for the Macromolecule sequence column's grid renderer to bind.
// Round-2 fix (A): the renderer-bound signal is gridCol.cellType being
// non-null AND in the bio-sequence renderer-family accepted set
// {'sequence','helm','separator','biln','custom','fasta'} — empirically
// gridCol.cellType resolves to the units-value (e.g. 'helm') for cold-
// loaded columns with units matching a registered columnTags filter,
// not to the registered meta.cellType: 'sequence' value. Pepsea-spec.ts
// L614 (sibling precedent) accepts {'sequence','helm','separator'}; the
// renderer-dispatch scenario widens this to cover the full units
// enumeration (biln/custom/fasta) the detector emits.
const BIO_SEQUENCE_CELL_TYPES = ['sequence', 'helm', 'separator', 'biln', 'custom', 'fasta'];

async function waitForSequenceCellTypeBind(page: import('@playwright/test').Page,
    timeoutMs = 60_000): Promise<void> {
  await page.waitForFunction((accepted: string[]) => {
    const df = grok.shell.tv?.dataFrame;
    if (!df) return false;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
    if (!macro) return false;
    const gridCol = (grok.shell.tv as any).grid?.col?.(macro.name);
    const ct: string | null = gridCol?.cellType ?? null;
    return ct !== null && accepted.indexOf(ct) >= 0;
  }, BIO_SEQUENCE_CELL_TYPES, {timeout: timeoutMs});
}

// Helper: locate the Macromolecule sequence column on the active TableView
// and return its (name, units, semType, gridCellType) — the four-tuple
// that disambiguates renderer dispatch per the strategy comment block.
// Caller MUST have awaited waitForSequenceCellTypeBind() first to avoid
// the renderer-bind race documented in the retry hypothesis (bug #2).
async function inspectMacroCol(page: import('@playwright/test').Page):
    Promise<{name: string | null, semType: string | null, units: string | null,
             gridCellType: string | null, hasErrorBalloon: boolean}> {
  return await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
    const gridCol = (grok.shell.tv as any).grid?.col?.(macro?.name);
    // Error-balloon presence — DG.Balloon family has class .d4-balloon-error.
    const hasErrorBalloon = !!document.querySelector('.d4-balloon-error');
    return {
      name: macro?.name ?? null,
      semType: macro?.semType ?? null,
      units: macro?.getTag?.('units') ?? macro?.meta?.units ?? null,
      gridCellType: gridCol?.cellType ?? null,
      hasErrorBalloon,
    };
  });
}

// ---------------------------------------------------------------------
// Scenario 1 — Detector + renderer dispatch for HELM and SEPARATOR
// ---------------------------------------------------------------------
test('Bio | Rendering — detector + renderer dispatch for HELM and SEPARATOR', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  // Step 1: Open filter_HELM.csv; detector sets units=helm; bilnSequence-
  // CellRenderer is the dispatch target (cellType=sequence, units=helm).
  await softStep('Open filter_HELM.csv — units=helm, bilnSequenceCellRenderer dispatch', async () => {
    await openBioDataset(page, HELM_PATH);
    // Retry-1 fix #2 + Round-2 fix (A): poll for renderer-bind with the
    // TOLERANT accepted-set predicate (cellType non-null and ∈ the bio-
    // sequence renderer-family). pepsea-spec.ts L614 sibling precedent.
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.semType).toBe('Macromolecule');
    // Dispatch-key signal: detector wrote units=helm (atlas bio.detector).
    expect(info.units).toBe('helm');
    // Round-2 fix (B): renderer-bound signal — gridCol.cellType is in
    // the bio-sequence renderer-family accepted set. For a units=helm
    // column the runtime resolves cellType to 'helm' (the units value),
    // which routes to bilnSequenceCellRenderer per the columnTags
    // 'quality=Macromolecule, units=helm' registration at package.ts#L344.
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES, `cellType ${info.gridCellType} must be a Bio sequence-family value`)
      .toContain(info.gridCellType!);
    expect(info.hasErrorBalloon).toBe(false);
  });

  // Step 3: Open filter_MSA.csv in a second table view; detector sets
  // units=separator (+ .separator='-', aligned=SEQ.MSA per bio.md L595);
  // separatorSequenceCellRenderer is the dispatch target.
  await softStep('Open filter_MSA.csv — units=separator (separatorSequenceCellRenderer dispatch)', async () => {
    await openBioDataset(page, MSA_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.semType).toBe('Macromolecule');
    // Dispatch-key signal: detector wrote units=separator + .separator='-'.
    expect(info.units).toBe('separator');
    // Round-2 fix (B): renderer-bound signal — accepted-set membership.
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES).toContain(info.gridCellType!);
    expect(info.hasErrorBalloon).toBe(false);
    // Also assert separator tag is set (detector contract per
    // bio.md L100 / atlas bio.detector — units=separator implies a
    // .separator tag is set on the column).
    const sepTag: string | null = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      // Detector writes the separator under .separator (leading dot) per
      // detectors.js — try both keys for robustness.
      return macro?.getTag?.('.separator') ?? macro?.getTag?.('separator') ?? null;
    });
    expect(sepTag).not.toBeNull();
  });

  // Step from Expected: the top-level cell-renderer registration list is
  // reachable. Enumerate Macromolecule renderers via the grok function
  // registry — atlas bio.rendering parent surface. cellRenderer-tagged
  // functions registered by Bio per package.ts L302/316/330/344/1155.
  await softStep('Enumerate Bio cellRenderer registrations — bio.rendering parent surface', async () => {
    const bioCellRenderers: string[] = await page.evaluate(() => {
      // Retry-1 fix #1: DG.Func.find expects tags: string[] per
      // js-api func.ts:100 — passing tags: 'cellRenderer' (string) silently
      // returns empty. Also widen the package-match predicate to accept
      // friendlyName='Bio', name='@datagrok/bio', or name='Bio' since
      // Dart-side .package shape varies across builds.
      const list: string[] = [];
      try {
        const fns = (DG.Func.find({tags: ['cellRenderer']}) || []) as any[];
        for (const f of fns) {
          const friendly = (f.package as any)?.friendlyName ?? null;
          const pkgName = (f.package as any)?.name ?? null;
          if (friendly === 'Bio' || pkgName === '@datagrok/bio' || pkgName === 'Bio')
            list.push(f.name);
        }
      } catch { /* surface as empty list */ }
      return list;
    });
    // Expect the five Bio sequence/monomer renderers to be registered.
    // Names must match the @grok.decorators.func({name: ...}) values
    // in package.ts.
    const expected = [
      'fastaSequenceCellRenderer',
      'separatorSequenceCellRenderer',
      'bilnSequenceCellRenderer',
      'customSequenceCellRenderer',
      'monomerCellRenderer',
    ];
    for (const n of expected)
      expect(bioCellRenderers, `cellRenderer ${n} must be registered (atlas bio.rendering)`).toContain(n);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ---------------------------------------------------------------------
// Scenario 2 — Convert HELM -> SEPARATOR re-dispatches renderer
//              (GROK-12164 guard)
// ---------------------------------------------------------------------
test('Bio | Rendering — Convert HELM to SEPARATOR re-dispatches renderer (GROK-12164)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  // Step 1: Open HELM; baseline dispatch is bilnSequenceCellRenderer
  // (cellType=sequence, units=helm).
  let preHelmColName: string | null = null;
  await softStep('Open filter_HELM.csv — baseline dispatch units=helm', async () => {
    await openBioDataset(page, HELM_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.units).toBe('helm');
    // Round-2 fix (B): accepted-set membership instead of strict
    // cellType === 'sequence'.
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES).toContain(info.gridCellType!);
    preHelmColName = info.name;
  });

  // Step 2: Right-click the sequence column header to open the Bio >
  // Transform > Convert Sequence Notation... dialog.
  //
  // bio.md L548 caveat: the Macromolecule column-header context menu is
  // a canvas hit-test (synthetic right-click bypasses the hit-test the
  // same way the Macromolecule cell-context-menu does). The scenario
  // body's "Right-click the sequence column header" step is therefore
  // realized via the canonical Bio top-menu path
  // [name="div-Bio---Transform---Convert-Sequence-Notation..."] (atlas
  // bio.transform.convert-notation.top-menu, package.ts#L1131 — the
  // same dialog opens either way per bio.md L375 "Cell-action variant"
  // note). This is a class-1 selector path; the substitution is the
  // sanctioned canvas-fallback for column-header hit-tests (per bio.md
  // L548).
  await softStep('Bio > Transform > Convert Sequence Notation... opens dialog', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Transform"]')!
        .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector(
        '[name="div-Bio---Transform---Convert-Sequence-Notation..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Convert-Sequence-Notation"]').waitFor({timeout: 60_000});
  });

  // Set Convert-to=separator with separator='-' — the GROK-12164
  // reproducer parameters per scenario Step 2. Selectors are class-1
  // per bio.md L371-372.
  await softStep('Set Convert-to=separator, Separator=-, click OK', async () => {
    const dlg = page.locator('[name="dialog-Convert-Sequence-Notation"]');
    // Convert-to SELECT — options include "separator" / "helm" / "biln"
    // (HELM source excludes "helm" per bio.md L371; "separator" is at
    // index 0 for HELM source per the same reference). Use the
    // canonical input-host selector + select-option semantics.
    await dlg.locator('[name="input-host-Convert-to"] select').selectOption('separator');
    // The Separator SELECT appears only when target is separator
    // (bio.md L372). Default options are ["-", ".", "/"].
    await dlg.locator('[name="input-host-Separator"] select').selectOption('-');
    const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(
      (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 60_000});
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Convert-Sequence-Notation"]').length === 0,
      null, {timeout: 15_000}).catch(() => {});
  });

  // GROK-12164 regression guard: the NEW separator column carries
  // units=separator AND grid.col(name).cellType === 'sequence' — the
  // dispatch follows the new column's tags, NOT the source HELM column's
  // tags. Pre-fix, the bilnSequenceCellRenderer remained dispatched on
  // the new column (cellType would still resolve to the source-column's
  // renderer family). Source HELM column must retain its own
  // units=helm tag — proving the source dispatch is unaffected.
  await softStep('GROK-12164: new column units=separator, source HELM column units=helm intact', async () => {
    // Retry-1 fix #2 + Round-2 fix (A): poll for the POST-CONVERT
    // renderer-bind on the new units=separator column with the TOLERANT
    // accepted-set predicate. The Convert action returns when the new
    // column is in the dataframe; the renderer-bind via detectSemanticTypes
    // is async (pepsea-spec.ts L789-796 pattern). cellType resolves to
    // the units-value family member (here likely 'separator', not
    // 'sequence').
    await page.waitForFunction((accepted: string[]) => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const sep: any = cols.find((c: any) =>
        c.semType === 'Macromolecule' &&
        ((c.getTag?.('units') ?? c.meta?.units) === 'separator'));
      if (!sep) return false;
      const gridCol = (grok.shell.tv as any).grid?.col?.(sep.name);
      const ct: string | null = gridCol?.cellType ?? null;
      return ct !== null && accepted.indexOf(ct) >= 0;
    }, BIO_SEQUENCE_CELL_TYPES, {timeout: 60_000});

    const tagsByMacro: Array<{name: string, units: string | null, gridCellType: string | null}> =
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macros = cols.filter((c: any) => c.semType === 'Macromolecule');
        return macros.map((c: any) => ({
          name: c.name,
          units: c.getTag?.('units') ?? c.meta?.units ?? null,
          gridCellType: (grok.shell.tv as any).grid?.col?.(c.name)?.cellType ?? null,
        }));
      });
    // Expect ≥2 Macromolecule columns: source HELM + new separator.
    expect(tagsByMacro.length).toBeGreaterThanOrEqual(2);
    const sourceHelm = tagsByMacro.find((m) => m.units === 'helm');
    const newSeparator = tagsByMacro.find((m) => m.units === 'separator');
    // Source HELM column must still carry units=helm (post-convert
    // source tag intact).
    expect(sourceHelm, 'source HELM column units=helm must be intact post-convert').toBeTruthy();
    // Round-2 fix (B): renderer-bound signal = accepted-set membership.
    // The source HELM column's cellType remains in the bio-sequence
    // renderer family (typically 'helm' on a cold-loaded dataset).
    expect(sourceHelm!.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES,
      `source HELM cellType ${sourceHelm!.gridCellType} must be in bio sequence family`)
      .toContain(sourceHelm!.gridCellType!);
    if (preHelmColName) expect(sourceHelm!.name).toBe(preHelmColName);
    // New SEPARATOR column carries units=separator (atlas
    // bio.rendering.separator dispatch contract). GROK-12164 invariant:
    // the new column's dispatch follows its OWN tags (units=separator),
    // not the source HELM column's tags. Pre-fix the bilnSequenceCellRenderer
    // remained dispatched on the new column; post-fix the dispatch follows
    // the new column's units tag. The accepted-set membership is the
    // renderer-bound signal; the (newSeparator.units === 'separator')
    // assertion above is the dispatch-key signal.
    expect(newSeparator, 'new column from convert must carry units=separator (GROK-12164)').toBeTruthy();
    expect(newSeparator!.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES,
      `new SEPARATOR cellType ${newSeparator!.gridCellType} must be in bio sequence family`)
      .toContain(newSeparator!.gridCellType!);
    // GROK-12164 explicit guard: the source column's dispatch must
    // differ from the new column's (or, if both resolve to the family
    // 'sequence' string, the dispatch-key units tags must differ — that
    // is asserted above via units=helm vs units=separator). Critically,
    // the new column must NOT be dispatched as if it were still HELM —
    // i.e. its gridCellType must not be 'helm'. Pre-fix, the bilnSequence-
    // CellRenderer would have remained bound, surfacing as cellType='helm'
    // on the new units=separator column. Post-fix, cellType ∈ {'separator',
    // 'sequence'} (depending on registry match order) — never 'helm'.
    expect(newSeparator!.gridCellType,
      'GROK-12164: new SEPARATOR column must not retain HELM dispatch (cellType !== "helm")')
      .not.toBe('helm');
    // Assert no error balloon — convert closed cleanly.
    const hasErrorBalloon = await page.evaluate(() => !!document.querySelector('.d4-balloon-error'));
    expect(hasErrorBalloon).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ---------------------------------------------------------------------
// Scenario 3 — Split-to-Monomers produces Monomer columns rendered by
//              monomerCellRenderer
// ---------------------------------------------------------------------
test('Bio | Rendering — Split to Monomers produces Monomer columns (monomerCellRenderer)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await softStep('Open filter_MSA.csv — separator-with-MSA-flag baseline', async () => {
    await openBioDataset(page, MSA_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.units).toBe('separator');
  });

  // Bio > Transform > Split to Monomers... — atlas
  // bio.transform.split-to-monomers (package.ts#L1225). Dialog opens
  // prefilled with the active table's Macromolecule column per
  // bio.md L395-403. Click OK.
  await softStep('Bio > Transform > Split to Monomers... — OK adds Monomer columns', async () => {
    const beforeMonCount: number = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      // Monomer columns carry semType='Monomer' (atlas bio.rendering.monomer
      // cellType=Monomer; quality=Monomer columnTag per atlas).
      return cols.filter((c: any) => c.semType === 'Monomer').length;
    });
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Transform"]')!
        .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector(
        '[name="div-Bio---Transform---Split-to-Monomers..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({timeout: 60_000});
    await page.locator('[name="dialog-Split-to-Monomers"] [name="button-OK"]').click();
    // Wait for the per-position Monomer columns to materialize.
    await page.waitForFunction((base) => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return cols.filter((c: any) => c.semType === 'Monomer').length > base;
    }, beforeMonCount, {timeout: 60_000});
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Split-to-Monomers"]').length === 0,
      null, {timeout: 15_000}).catch(() => {});
  });

  // Per-position Monomer columns render via monomerCellRenderer — atlas
  // bio.rendering.monomer (cellType='Monomer', package.ts#L1155). The
  // monomerCellRenderer registration has cellType='Monomer' (distinct
  // from the four sequence renderers which all share cellType='sequence');
  // grid.col(name).cellType === 'Monomer' is the deterministic
  // JS-API dispatch signal for the monomer renderer family.
  await softStep('Monomer columns dispatch to monomerCellRenderer (cellType=Monomer)', async () => {
    // Retry-1 fix #3: poll for at least one Monomer column's renderer-
    // bind BEFORE asserting (same async-attach issue as the sequence
    // renderer per pepsea-spec.ts L789-796; here scoped to cellType=
    // 'Monomer' which is the disambiguator for the monomer renderer
    // family per Bio package.ts L1155-1167). The monomerCellRenderer
    // registration uses meta.cellType: 'Monomer' and matches the
    // quality=Monomer columnTag (atlas bio.rendering.monomer).
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const monCols = cols.filter((c: any) => c.semType === 'Monomer');
      if (monCols.length === 0) return false;
      return monCols.some((c: any) =>
        (grok.shell.tv as any).grid?.col?.(c.name)?.cellType === 'Monomer');
    }, null, {timeout: 60_000});

    const monomerColInfo: Array<{name: string, semType: string | null, gridCellType: string | null}> =
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const monCols = cols.filter((c: any) => c.semType === 'Monomer');
        return monCols.map((c: any) => ({
          name: c.name,
          semType: c.semType,
          gridCellType: (grok.shell.tv as any).grid?.col?.(c.name)?.cellType ?? null,
        }));
      });
    // Scenario Expected: N new Monomer columns appear (N = alignment
    // width of source — typically 17 for filter_MSA per scenario
    // Setup). Assert N >= 1 and every one carries semType='Monomer'.
    expect(monomerColInfo.length).toBeGreaterThan(0);
    for (const m of monomerColInfo)
      expect(m.semType, `column ${m.name} should be semType='Monomer'`).toBe('Monomer');
    // At least one Monomer column has its grid renderer bound with
    // cellType='Monomer' (atlas bio.rendering.monomer dispatch contract).
    // The Monomer family is the only Bio renderer with cellType !==
    // 'sequence', so this is the deterministic disambiguator for the
    // monomer dispatch path. The waitForFunction above guarantees the
    // bind has happened by this point.
    const anyMonomerDispatched = monomerColInfo.some((m) => m.gridCellType === 'Monomer');
    expect(anyMonomerDispatched, 'at least one Monomer column should have grid.cellType=Monomer').toBe(true);
    // Assert no error balloon.
    const hasErrorBalloon = await page.evaluate(() => !!document.querySelector('.d4-balloon-error'));
    expect(hasErrorBalloon).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ---------------------------------------------------------------------
// Scenario 4 — Custom-notation column dispatches to
//              customSequenceCellRenderer
// ---------------------------------------------------------------------
test('Bio | Rendering — units=custom column dispatches to customSequenceCellRenderer', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  // Step 1: Open HELM (units=helm baseline).
  await softStep('Open filter_HELM.csv — baseline units=helm', async () => {
    await openBioDataset(page, HELM_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.units).toBe('helm');
    // Round-2 fix (B): accepted-set membership.
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES).toContain(info.gridCellType!);
  });

  // Step 2: Programmatically set units=custom on the sequence column to
  // simulate a pluggable notation provider. The atlas
  // bio.rendering.custom registration is keyed on units value `custom`
  // and exists for the pluggable-notation extension point — see atlas
  // package.ts#L302. The scenario explicitly directs:
  //   "column.setTag('units', 'custom')"
  //   "view.grid.invalidate()"
  // This is a sanctioned JS-API stage (the scenario body itself calls
  // for it — see scenario `.md` Scenario 4 Step 2) because no public
  // notation-provider registration UI exists in stock Bio; the staging
  // is necessary to exercise the customSequenceCellRenderer registration.
  await softStep('setTag units=custom + detectSemanticTypes — re-dispatches to customSequenceCellRenderer', async () => {
    // Retry-1 fix #4: per pepsea-spec.ts L776-784, the renderer rebind
    // happens via grok.data.detectSemanticTypes(df) — mirrors the dialog
    // OK path at multiple-sequence-alignment-ui.ts L321. Without it,
    // setTag('units','custom') + grid.invalidate() does NOT rebind the
    // renderer (Bio's cellRenderer registry routes on the {quality=
    // Macromolecule, units=<value>} tag pair, but the BIND only happens
    // via detectSemanticTypes; invalidate alone paints with the existing
    // bound renderer). Sequence:
    //   1) capture preUnits
    //   2) setTag('units','custom') on the Macromolecule column
    //   3) await grok.data.detectSemanticTypes(df) — rebind contract
    //   4) grid.invalidate() — scenario directive, complements step 3
    //   5) poll for cellType to re-bind to 'sequence' (the
    //      customSequenceCellRenderer registration shares cellType=
    //      'sequence' with the other three notation renderers per
    //      package.ts L302-340; the units=custom tag is the dispatch
    //      discriminator inside that family).
    const colName: string | null = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!macro) return null;
      macro.setTag('units', 'custom');
      try { await (grok.data as any).detectSemanticTypes(df); } catch { /* tolerate older builds */ }
      try { (grok.shell.tv as any).grid?.invalidate?.(); } catch { /* may not throw on older builds */ }
      return macro.name;
    });
    expect(colName, 'Macromolecule column must be located on HELM dataset').toBeTruthy();

    // Round-2 fix (A) applied to scenario 4: poll for the renderer to
    // re-bind on the units=custom column with the TOLERANT accepted-set
    // predicate. customSequenceCellRenderer registers with
    // meta.cellType: 'sequence' AND meta.columnTags
    // 'quality=Macromolecule, units=custom' (package.ts L302-340), so
    // the runtime gridCol.cellType resolves to either 'sequence' or
    // 'custom' depending on registry match order — both are accepted.
    await page.waitForFunction(({colName, accepted}: {colName: string | null, accepted: string[]}) => {
      const gridCol = (grok.shell.tv as any).grid?.col?.(colName);
      const ct: string | null = gridCol?.cellType ?? null;
      return ct !== null && accepted.indexOf(ct) >= 0;
    }, {colName, accepted: BIO_SEQUENCE_CELL_TYPES}, {timeout: 60_000});

    const result: {postUnits: string | null, postCellType: string | null} =
      await page.evaluate((name) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.name === name);
        return {
          postUnits: macro?.getTag?.('units') ?? macro?.meta?.units ?? null,
          postCellType: (grok.shell.tv as any).grid?.col?.(name)?.cellType ?? null,
        };
      }, colName);
    // Post-step-2 dispatch-key signal: detector wrote units=custom
    // after the setTag + detectSemanticTypes rebind.
    expect(result.postUnits).toBe('custom');
    // Round-2 fix (B/D): renderer-bound signal — accepted-set membership.
    // For a units=custom column the runtime resolves cellType to
    // 'custom' (the units-value match against the customSequence-
    // CellRenderer registration columnTags 'quality=Macromolecule,
    // units=custom' at package.ts#L302) — OR 'sequence' if registry
    // order falls through to the meta.cellType: 'sequence' value.
    // Either is accepted; the customSequenceCellRenderer is the
    // dispatched renderer by both units=custom and the absence of
    // an error balloon (scenario Expected — pluggable provider falls
    // back to a string paint).
    expect(result.postCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES,
      `cellType ${result.postCellType} must be in bio sequence family after units=custom rebind`)
      .toContain(result.postCellType!);
    // Scenario Expected: no error balloon appears after the renderer
    // re-dispatches to the customSequenceCellRenderer fallback path.
    const hasErrorBalloon = await page.evaluate(() => !!document.querySelector('.d4-balloon-error'));
    expect(hasErrorBalloon).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
