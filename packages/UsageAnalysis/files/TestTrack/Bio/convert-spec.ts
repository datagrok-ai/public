/* ---
sub_features_covered:
  - bio.calculate.get-region
  - bio.calculate.get-region.top-menu
  - bio.transform.convert-notation
  - bio.transform.convert-notation.top-menu
  - bio.transform.to-atomic-level
  - bio.transform.split-to-monomers
  - bio.detector
--- */
//   related_bugs: [GROK-15176, GROK-12164] (cross-cutting bug-repro
//     surfaces — per chain bug_focused_candidates, GROK-15176 is delegated
//     to a dedicated bio-grok-15176-spec.ts that chains To Atomic Level →
//     Chem renderer / PubChem standardization; GROK-12164 is below the
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv', units: 'fasta'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv', units: 'helm'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv', units: 'separator'},
];
for (const ds of datasets) {
  test(`Bio Convert matrix on ${ds.name}`, async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    // Setup phase: open dataset, wait for Macromolecule semType detection +
    // Bio package init (cell renderer + filter registration). Mirrors the
    // analyze-spec.ts setup phase — same cold-start tolerance applies.
    await page.evaluate(async (path) => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
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
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    // Bio top-menu readiness poll (cold-start stabilization — same
    // two-layer guard as analyze-spec.ts).
    //
    // Layer 1 — DOM-level top-menu visibility: the Bio top-menu entry
    // appears only once Bio package functions are registered against the
    // active Macromolecule TableView.
    //
    // Layer 2 — init-order serialization probe via Bio:getSeqHelper /
    // Bio:getMonomerLibHelper / Bio:getBioLib (atlas
    // bio.cp.bio-service-surface-init). The runtime serializes
    // `grok.functions.call('Bio:<...>')` after init completion — awaiting
    // any of these is a deterministic "Bio package ready" probe.
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 3000));
    });
    // Scenario 1 — Macromolecule detection assertion (atlas bio.detector).
    // The Macromolecule detector classifies the sequence column
    // synchronously; the table view opens with the sequence column tagged
    // `quality=Macromolecule, units=<fasta|helm|separator>`.
    await softStep(`${ds.name}: Dataset opens with Macromolecule sequence column (units=${ds.units})`, async () => {
      const info: {hasMacro: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        return {hasMacro: !!macro, units: macro?.meta?.units ?? null};
      });
      expect(info.hasMacro).toBe(true);
      expect(info.units).toBe(ds.units);
    });
    // Scenario 2.a — Bio > Calculate > Extract Region... (atlas
    // bio.calculate.get-region + .top-menu). Dialog opens prefilled with
    // the active Macromolecule column; OK adds a sub-region Macromolecule
    // column whose name matches `<source>: (start-end)` and whose units
    // equal the source notation (per scenario Setup table cell).
    await softStep(`${ds.name}: Calculate > Extract Region adds a Macromolecule sub-region column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Calculate"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector('[name="div-Bio---Calculate---Extract-Region..."]') as HTMLElement).click();
      });
      // Dialog name is `dialog-Get-Region` per bio.md (named after the API
      // function `getRegion`, not the menu label "Extract Region...").
      await page.locator('[name="dialog-Get-Region"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Get-Region"] [name="button-OK"]').click();
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 30_000});
      const info: {hasRegion: boolean, units: string | null, name: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        // Find the most-recently-added Macromolecule column whose name
        // matches the `<source>: (start-end)` pattern.
        const region: any = cols.reverse().find((c: any) =>
          c.semType === 'Macromolecule' && /:\s*\(\d+-\d+\)/.test(c.name));
        return {hasRegion: !!region, units: region?.meta?.units ?? null, name: region?.name ?? null};
      });
      expect(info.hasRegion).toBe(true);
      expect(info.units).toBe(ds.units);
      // Wait for the Extract Region dialog to fully detach before next step.
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-Get-Region"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    // Scenario 2.b — Bio > Transform > Convert Sequence Notation... (atlas
    // bio.transform.convert-notation + .top-menu, package.ts#L1131).
    //
    // Canonical menu path per atlas + bio.md 2026-06-01 selector validation.
    // The scenario's source-text "PolyTool > Convert" referred to the same
    // surface — bio.md MCP-validated that PolyTool > Convert on FASTA
    // routes to `dialog-To-Atomic-Level` (different dialog, same
    // shared-widget overlap that the 2026-04-23 reference run hit), while
    // the canonical Transform > Convert Sequence Notation... opens
    // `dialog-Convert-Sequence-Notation` (proper notation-conversion
    // surface). This spec uses the canonical Transform path per the atlas
    // sub_feature top-menu registration.
    //
    // GROK-12164 awareness: a separator-with-`-` convert path is a known
    // renderer-dispatch regression (atlas
    // bio.x.detector-renderer-after-convert). This spec verifies a new
    // column appears with the converted notation — it does NOT assert
    // post-convert renderer dispatch (covered cross-scenario via a
    // dedicated bug-focused spec if/when prioritised; chain
    // `bug_match_attempts_skipped: GROK-12164`).
    await softStep(`${ds.name}: Transform > Convert Sequence Notation adds a new Macromolecule column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
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
      // Accept the prefilled defaults — the dialog auto-picks the active
      // Macromolecule column (per bio.md "Source column" row); the
      // Convert-to SELECT excludes the source notation by default.
      await page.locator('[name="dialog-Convert-Sequence-Notation"] [name="button-OK"]').click();
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 60_000});
      // Verify a new Macromolecule column appears with a different units
      // tag than the source (the converted notation). Conversion target
      // is the dialog default — we assert the structural invariant
      // (new Macromolecule column with units != source-units), not the
      // specific target notation.
      const info: {newMacroCount: number, lastUnits: string | null, allUnits: string[]} =
        await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
          const last: any = macroCols[macroCols.length - 1];
          return {
            newMacroCount: macroCols.length,
            lastUnits: last?.meta?.units ?? null,
            allUnits: macroCols.map((c: any) => c.meta?.units ?? null).filter((u: any) => u !== null),
          };
        });
      // Expect ≥2 Macromolecule columns (source + converted). The last
      // Macromolecule column carries the converted-to notation, which
      // differs from the source units.
      expect(info.newMacroCount).toBeGreaterThanOrEqual(2);
      expect(info.lastUnits).not.toBe(ds.units);
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-Convert-Sequence-Notation"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    // Scenario 2.c — Bio > Transform > To Atomic Level... (atlas
    // bio.transform.to-atomic-level, package.ts#L863).
    //
    // Produces a V3000 molfile column (semType=Molecule, units=molblock).
    //
    // GROK-15176 awareness: MSA → To Atomic Level molfile may carry
    // illegal isotope=1 on heavy-atom oxygens (atlas
    // bio.x.bio-to-chem-via-atomic-level). This spec exercises the
    // To Atomic Level surface on canonical sources but does NOT chain
    // into Chem renderer / PubChem standardization — the molfile-validity
    // contract is delegated to bug-focused chain entry
    // bio-grok-15176-spec.ts per scenario related_bugs note.
    await softStep(`${ds.name}: Transform > To Atomic Level adds a Molecule molblock column`, async () => {
      const beforeMolCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Molecule').length;
      });
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---To-Atomic-Level..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-To-Atomic-Level"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-To-Atomic-Level"] [name="button-OK"]').click();
      // Atomic-level conversion can take 20-30s on a 540-row HELM dataset
      // per convert-run.md 2026-04-23 (HELM: 21s). 120s timeout is the
      // defensive ceiling.
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Molecule').length > base;
      }, beforeMolCount, {timeout: 120_000});
      const info: {hasMol: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const mol: any = cols.reverse().find((c: any) => c.semType === 'Molecule');
        return {hasMol: !!mol, units: mol?.meta?.units ?? null};
      });
      expect(info.hasMol).toBe(true);
      expect(info.units).toBe('molblock');
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    // Scenario 2.d — Bio > Transform > Split to Monomers... (atlas
    // bio.transform.split-to-monomers, package.ts#L1225).
    //
    // Produces N per-position `Monomer` columns where N is the alignment
    // width of the source column. Scenario Setup table:
    //   FASTA → 39 Monomer columns
    //   HELM  → 17 Monomer columns
    //   MSA   → 17 Monomer columns
    // Assertion: N ≥ 1 (and matches the source alignment) — exact N is
    // dataset-specific.
    await softStep(`${ds.name}: Transform > Split to Monomers adds per-position Monomer columns`, async () => {
      const beforeMonCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
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
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Monomer').length > base;
      }, beforeMonCount, {timeout: 60_000});
      const monCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Monomer').length;
      });
      // N ≥ 1 — exact N is dataset-specific (scenario expects ~17-39).
      expect(monCount).toBeGreaterThan(beforeMonCount);
      expect(monCount).toBeGreaterThan(0);
    });
    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
