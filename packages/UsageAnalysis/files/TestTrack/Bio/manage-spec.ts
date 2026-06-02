/* ---
sub_features_covered:
  - bio.manage.libraries-view
  - bio.detector
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: smoke — single-action ui-smoke)
//   sub_features_covered: [bio.manage.libraries-view, bio.detector]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [] (chain analyzer skipped — manage.* surfaces have no
//     curated bug-library entries)
//   produced_from: migrated
//   coverage_type: smoke
//   unresolved_ambiguities: [ui-smoke-rule-1-trigger-fit-poor] (carried
//     forward; out-of-scope for this spec — chain-analyzer concern)
//
// Atlas provenance (consulted):
//   feature-atlas/bio.yaml#sub_features[bio.manage.libraries-view]
//     (package.ts#L1359 → manageLibrariesView top-menu function).
//   feature-atlas/bio.yaml#sub_features[bio.detector] (touched in setup
//     by the canonical HELM fixture's Macromolecule semType classification).
//
// Reference template: sibling bio/manage-spec.ts (prior cycle on HELM.csv)
// supplies the canonical body shape: setup → top-menu dispatch →
// view-name assertion → checkbox enumeration → toggle-and-restore.
// Selectors come from .claude/skills/grok-browser/references/bio.md#L459
// (MCP-validated 2026-06-01 on dev.datagrok.ai), specifically the
// `.monomer-lib-controls-form .ui-input-bool input[type="checkbox"]`
// per-library checkbox structure. The scenario .md's hint of
// `[name="input-host-<libfile>.json"]` is not in the bio.md reference
// — we follow the canonical reference per §"Scenario authority"
// (scenario itself cites bio.md as the selector source).
//
// Server-state independence: the per-library list depends on FileShare
// permissions on System:AppData/Bio/monomer-libraries; bio.md observed
// 4 libraries on the recon account but explicitly notes "varies per-user".
// Per scenario .md "library set is server-state dependent — assert
// structure, not a fixed library count": we assert ≥1 per-library
// checkbox host + per-row state isolation (toggle two distinct rows
// and check independence), NOT a fixed count.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Bio Manage Monomer Libraries (filter_HELM.csv)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Setup phase: open canonical HELM fixture, wait for Macromolecule
  // detector to classify the sequence column, wait for Bio package
  // cell-renderer / filter registration to settle. Mirrors the canonical
  // setup pattern from sibling specs (analyze-spec.ts, sequence-space-spec.ts).
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_HELM.csv');
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
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Bio top-menu readiness poll (cold-start stabilization — Bio package
  // sometimes needs an additional beat after grid canvas paints before
  // the top-menu Bio leaf is reliably enumerable on a cold init).
  await page.waitForFunction(() => {
    return !!document.querySelector('[name="div-Bio"]');
  }, null, {timeout: 30_000});

  await softStep('Open filter_HELM.csv and detect Macromolecule semType (bio.detector)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macromolCol = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        rows: df.rowCount,
        hasMacromolecule: !!macromolCol,
        macromolColName: macromolCol?.name || null,
      };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.hasMacromolecule).toBe(true);
    expect(info.macromolColName).not.toBeNull();
  });

  await softStep('Dispatch Bio > Manage > Monomer Libraries top-menu (bio.manage.libraries-view)', async () => {
    // Top-menu navigation pattern: click Bio root, hover Manage sub-menu
    // to populate its children, then click the leaf. Mirrors the canonical
    // pattern from sibling manage-spec.ts (prior cycle) and bio.md
    // reference selector `[name="div-Bio---Manage---Monomer-Libraries"]`.
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 500));
      const manage = document.querySelector('[name="div-Bio---Manage"]')!;
      manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 500));
      (document.querySelector('[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
      null, {timeout: 30_000});

    const view = await page.evaluate(() => {
      const v = grok.shell.v;
      const root = v?.root;
      // Section headers: scenario asserts both `Manage Monomer Libraries`
      // and `Manage Duplicate Monomer Symbols` headers are present in
      // the view. Empirical DOM shape (MCP-evidenced via Playwright
      // trace on prior attempt's error-context.md, 2026-06-01):
      //   - `Manage Duplicate Monomer Symbols` IS a real <h1>: a heading
      //     element whose textContent equals the title verbatim.
      //   - `Manage Monomer Libraries` is the View tab/panel TITLE — NOT
      //     an <hN> heading. It surfaces as `grok.shell.v.name` and as
      //     a `[cursor=pointer]` generic in the view header band. There
      //     is no exact-text-equal-to-"Manage Monomer Libraries"
      //     <hN>/<div>; the prior matcher used `textContent === ...` on
      //     a `div` set and missed because containing elements wrap the
      //     title alongside a counter badge.
      //
      // Detection contract revised: match `Manage Monomer Libraries`
      // by view.name (the canonical surface — same path Datagrok uses
      // to render the tab title), and match `Manage Duplicate Monomer
      // Symbols` by its <hN> heading element (durable structural anchor).
      let hasDuplicateHeading = false;
      if (root) {
        const headings = root.querySelectorAll('h1, h2, h3, h4');
        for (const el of Array.from(headings)) {
          const t = (el.textContent || '').trim();
          if (t === 'Manage Duplicate Monomer Symbols') {
            hasDuplicateHeading = true;
            break;
          }
        }
      }
      return {
        name: v?.name || null,
        type: v?.type || null,
        rootClass: root?.className || null,
        hasGrokView: !!root?.closest('.grok-view') || (root?.className || '').includes('grok-view'),
        hasDuplicateHeading,
      };
    });
    // `Manage Monomer Libraries` (the view's section/title) is asserted
    // via `grok.shell.v.name` — the canonical Datagrok handle. This
    // bypasses fragile DOM-text matching and aligns with the view-root
    // identification used in bio.md (L468: "where v.name === 'Manage
    // Monomer Libraries'").
    expect(view.name).toBe('Manage Monomer Libraries');
    expect(view.type).toBe('view');
    // `Manage Duplicate Monomer Symbols` section header — asserted as
    // a literal heading element (durable structural anchor: an actual
    // <h1> in the right-side panel).
    expect(view.hasDuplicateHeading).toBe(true);
  });

  // SCOPE_REDUCTION (within-slice): scenario step 3 asks to assert a
  // library Search box (`[name="input-host-Search"]`). That selector is
  // not in bio.md (the MCP-validated 2026-06-01 reference for the
  // Monomer Libraries view) and could not be live-MCP-observed this
  // cycle (profile auth stale despite prewarm — see dispatch yaml
  // mcp_observations). Per §"Selector provenance (3-class model)" we
  // do NOT emit class-3 (inferred / scenario-hinted) selectors. The
  // search-box presence assertion is deferred; the per-library
  // checkbox listing + per-row state isolation (the load-bearing
  // ui-smoke surface for `bio.manage.libraries-view`) IS asserted.
  await softStep('View renders per-library checkbox listing (>=1 row) + checkbox element per row', async () => {
    const listing = await page.evaluate(() => {
      const root = grok.shell.v.root;
      // Per-library rows (bio.md reference L470-471):
      //   `.monomer-lib-controls-form .ui-input-bool` (host)
      //   `.monomer-lib-controls-form .ui-input-bool input[type="checkbox"]` (control)
      const form = root.querySelector('.monomer-lib-controls-form');
      const libRows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
      const checkboxes: HTMLInputElement[] = libRows
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const labels = libRows.map((r) => {
        const span = r.querySelector('span');
        return span ? (span.textContent || '').trim() : '';
      });
      return {
        formPresent: !!form,
        rowCount: libRows.length,
        checkboxCount: checkboxes.length,
        labels,
        eachRowHasCheckbox: libRows.length > 0 &&
          libRows.every((r) => !!r.querySelector('input[type="checkbox"]')),
      };
    });
    // Class-1 assertions (per bio.md reference): per-library checkbox host
    // structure under `.monomer-lib-controls-form`.
    expect(listing.formPresent).toBe(true);
    expect(listing.rowCount).toBeGreaterThanOrEqual(1);
    expect(listing.checkboxCount).toBeGreaterThanOrEqual(1);
    // Per-row state isolation precondition: every visible row carries
    // an `input[type="checkbox"]` (the Datagrok ui-input-bool checkbox).
    // Scenario asserts this as the "checkbox element present inside
    // each library host" check.
    expect(listing.eachRowHasCheckbox).toBe(true);
  });

  await softStep('Toggle two libraries independently — per-row state isolation holds', async () => {
    // Scenario 2 assertion: toggling one library's checkbox transitions
    // it cleanly (checked→unchecked or vice versa) AND a second library's
    // checkbox remains independent (its state is unaffected by toggling
    // the first). After the test, restore the original states so the
    // test is idempotent against shared FileShare settings.
    const result = await page.evaluate(async () => {
      const root = grok.shell.v.root;
      const form = root.querySelector('.monomer-lib-controls-form');
      const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
      const cbs = rows
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      if (cbs.length < 2) {
        return {
          haveTwoRows: false,
          rowCount: cbs.length,
        };
      }
      const idxA = 0;
      const idxB = 1;
      const a0 = cbs[idxA].checked;
      const b0 = cbs[idxB].checked;
      // Toggle A
      cbs[idxA].click();
      await new Promise((r) => setTimeout(r, 2500));
      const cbs2 = (Array.from(form!.querySelectorAll('.ui-input-bool')) as HTMLElement[])
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const a1 = cbs2[idxA].checked;
      const b1 = cbs2[idxB].checked;
      // Toggle B (independently)
      cbs2[idxB].click();
      await new Promise((r) => setTimeout(r, 2500));
      const cbs3 = (Array.from(form!.querySelectorAll('.ui-input-bool')) as HTMLElement[])
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const a2 = cbs3[idxA].checked;
      const b2 = cbs3[idxB].checked;
      // Restore: toggle A back, toggle B back. Order chosen so both end
      // at original state regardless of platform debounce.
      cbs3[idxA].click();
      await new Promise((r) => setTimeout(r, 2000));
      cbs3[idxB].click();
      await new Promise((r) => setTimeout(r, 2500));
      const cbs4 = (Array.from(form!.querySelectorAll('.ui-input-bool')) as HTMLElement[])
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const aR = cbs4[idxA].checked;
      const bR = cbs4[idxB].checked;
      return {
        haveTwoRows: true,
        rowCount: cbs.length,
        a0, b0,
        a1, b1,
        a2, b2,
        aR, bR,
        aToggledOnFirstClick: a1 !== a0,
        bUnchangedAfterAToggle: b1 === b0,
        bToggledOnSecondClick: b2 !== b1,
        aUnchangedAfterBToggle: a2 === a1,
        aRestored: aR === a0,
        bRestored: bR === b0,
      };
    });
    // Scenario asserts a "second library checkbox independently" toggle
    // implying ≥2 libraries are available. On a degenerate single-library
    // FileShare we degrade to documenting via the stepErrors fail — but
    // dev.datagrok.ai consistently exposes ≥4 libraries (per bio.md L480).
    expect(result.haveTwoRows).toBe(true);
    // Independence assertion: toggling A changes A but not B; toggling B
    // changes B but not A. Together this is per-row state isolation.
    expect(result.aToggledOnFirstClick).toBe(true);
    expect(result.bUnchangedAfterAToggle).toBe(true);
    expect(result.bToggledOnSecondClick).toBe(true);
    expect(result.aUnchangedAfterBToggle).toBe(true);
    // Idempotent restore so subsequent runs/users don't inherit toggled state.
    expect(result.aRestored).toBe(true);
    expect(result.bRestored).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
