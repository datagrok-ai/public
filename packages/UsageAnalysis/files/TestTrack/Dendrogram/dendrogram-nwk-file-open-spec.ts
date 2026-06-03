/* ---
sub_features_covered: [dendrogram.fileviewer.newick, dendrogram.api.tree-helper.newick-to-df, dendrogram.viewer.dendrogram-app]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (no inference applied; per the prompt's missing-field
//     policy the constraint extraction uses non-pyramid defaults — JS-API
//     substitution permitted broadly, ≥1 DOM-driving call still REQUIRED to
//     satisfy E-LAYER-COMPLIANCE-01).
//   sub_features_covered: [dendrogram.fileviewer.newick,
//     dendrogram.api.tree-helper.newick-to-df,
//     dendrogram.viewer.dendrogram-app]
//   ui_coverage_responsibility: absent
//   related_bugs: []
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   dendrogram.yaml#critical_paths[dendrogram.cp.newick-file-open-via-files-browser]
//     derived_from: public/packages/Dendrogram/src/package.ts#L285
//   dendrogram.yaml#interactions[dendrogram.cross.newick-file-to-viewer]
//     derived_from: public/packages/Dendrogram/src/package.ts#L274-L294
//   dendrogram.yaml#sub_features[dendrogram.fileviewer.newick]
//     derived_from: public/packages/Dendrogram/src/package.ts#L297-L317
//   dendrogram.yaml#sub_features[dendrogram.api.tree-helper.newick-to-df]
//     derived_from: public/packages/Dendrogram/src/utils/tree-helper.ts#L30
//   dendrogram.yaml#sub_features[dendrogram.viewer.dendrogram-app]
//     derived_from: public/packages/Dendrogram/src/apps/dendrogram-app.ts#L9
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser reference):
//   `grok.shell.v.root.querySelector('canvas')` — PhylocanvasGL preview view's
//     canvas reached via `grok.shell.addView(view)` on the DG.View returned
//     by previewNewick. The view DOES NOT carry a `[name="viewer-PhylocanvasGL"]`
//     attribute — DG.View.fromRoot(viewerRoot) wraps the PhylocanvasGL plot
//     root without surfacing a `[name=...]` selector. Observed live
//     2026-06-03 via chrome-devtools MCP (evaluate_script on
//     dev.datagrok.ai/oahadzhanian). The dendrogram.md ref's "Forward TODO"
//     anticipates this — PhylocanvasGL preview selectors were not captured
//     in the 2026-06-03 recon round. The spec asserts canvas presence
//     under `grok.shell.v.root` (the active view's root) rather than a
//     named-selector match.
//   `grok.shell.tv.path === '/func/Dendrogram.dendrogramApp'` — the
//     TableView path set by DendrogramApp.buildView() (apps/dendrogram-app.ts#L84).
//     Used as a structural witness that DendrogramApp.init(df) ran and
//     opened the app's table view (in addition to the
//     [name="viewer-Dendrogram"] DOM probe). Observed live 2026-06-03.
//
// Selectors per .claude/skills/grok-browser/references/dendrogram.md (rev
// 2026-06-03 live-MCP-validated):
//   [name="viewer-Dendrogram"] (host container; ARIA region "Dendrogram"),
//   [name="viewer-Dendrogram"] canvas (single canvas; the rendered tree).
//
// MCP recon evidence (live 2026-06-03 on https://dev.datagrok.ai/, user
// oahadzhanian, Dendrogram plugin v1.4.11, viewport 1920x1080):
//   - The registered fileHandler function NAME is `Dendrogram:importNewick`
//     (NOT `Dendrogram:importNwk`). The `name:` metadata in the @fileHandler
//     decorator declares `importNwk` as the file-handler EXTENSION KEY, but
//     the function method itself is registered as `importNewick`. The
//     scenario .md mentions `Dendrogram:importNwk` in passing — the
//     authoritative function name on the server is `importNewick`. This is
//     consistent with `DG.Func.find({package: 'Dendrogram'})` enumeration.
//   - `grok.functions.call('Dendrogram:importNewick', {fileContent: NEWICK})`
//     returns `[]` (Promise<DG.DataFrame[]>; the handler always returns
//     empty per package.ts#L293 because DendrogramApp opens its own
//     TableView via grok.shell.addTableView and the result is not
//     surfaced back through the fileHandler return contract).
//   - Within ~3-5s of the call: `grok.shell.tv` becomes a 7-row DataFrame
//     named `Table` with columns ['node','parent','leaf','distance'];
//     tv.viewers enumerates ['Grid','Dendrogram']; the underlying df has
//     tag .newick = the original literal '((a:1,b:1):1,(c:1,d:1):1);' and
//     .newickJson populated (length ~252).
//   - The 4 leaf rows (where df.col('leaf').get(i) === true) carry
//     node names ['a','b','c','d'] (sorted set). 3 internal nodes.
//   - tv.path === '/func/Dendrogram.dendrogramApp' (set by
//     DendrogramApp.buildView at apps/dendrogram-app.ts#L84) — structural
//     witness for the app-mount.
//   - The [name="viewer-Dendrogram"] host mounts inline; canvas size
//     573x1000 with non-zero bounding rect.
//   - For Scenario 2: `grok.functions.call('Dendrogram:previewNewick',
//     {file: <FileInfo>})` returns a DG.View. The FileInfo MUST be
//     obtained via `grok.dapi.files.list('System:AppData/Dendrogram')` and
//     filtered for `fileName === 'sample.nwk'` — passing a path-string-only
//     synthetic produces a "FileSystemException: Not a directory" server
//     error. DG.FileInfo.fromString(text, name) is NOT a substitute for a
//     server-resident file because previewNewick calls file.readAsString()
//     which expects a backing path.
//   - The returned DG.View's root is `<div class="ui-panel grok-view">`
//     wrapping a single `<canvas>` (the PhylocanvasGL WebGL canvas). After
//     `grok.shell.addView(view)`, `grok.shell.v.root.querySelectorAll(
//     'canvas').length === 1` and the canvas's bounding rect is non-zero
//     (1123x1000 in viewport). NO [name="viewer-PhylocanvasGL"] attribute
//     is present anywhere — the wrapper is plain DG.View.fromRoot.
//   - `TreeHelper.newickToDf('((a:1,b:1):1,(c:1,d:1):1);', 'tmp')` returns
//     a 7-row DataFrame whose 4 leaf rows are ['a','b','c','d'] — set
//     semantics preserved.
//   - `grok.dapi.files.delete('System:AppData/Dendrogram/sample.nwk')`
//     succeeds; `grok.dapi.files.exists(...)` returns false post-delete.
//
// Reference template:
//   public/packages/UsageAnalysis/files/TestTrack/Dendrogram/dendrogram-viewer-from-newick-prop-spec.ts
//     (same section, same target_layer: playwright, same MCP-evidence
//      comment block pattern, softStep + stepErrors + cleanup pattern,
//      scoped helper functions for grok-shell setup).
//
// Authoring choice — JS-API substitution for the Files-browser open gesture:
//   pyramid_layer is absent so JS-API substitution is permitted (FORBIDDEN
//   "JS API substitution for owned UI flows" only applies when pyramid_layer
//   == ui-smoke). The scenario .md Step 2 explicitly licenses the JS-API
//   substitution via `grok.functions.call('Dendrogram:importNwk', ...)` as
//   "the exact same `importNewick` handler the file-open gesture dispatches
//   to". The canonical Files-browser UI gestures (double-click row to open,
//   single-click row to preview) have no documented [name=...] selectors in
//   `dendrogram.md` (the Forward TODO in that ref calls out this exact gap).
//   We exercise the registered handler functions directly — same code path
//   reached, with deterministic timing for `grok test`. E-LAYER-COMPLIANCE-01
//   (the >=1 DOM-driving-call invariant) is satisfied via
//   page.locator('[name="viewer-Dendrogram"]').waitFor() in Scenario 1 and
//   via the active-view canvas DOM probe in Scenario 2.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Canonical 4-leaf newick fixture, matching the scenario .md Setup.
const NEWICK = '((a:1,b:1):1,(c:1,d:1):1);';
const EXPECTED_LEAVES: ReadonlyArray<string> = ['a', 'b', 'c', 'd'];
const FIXTURE_PATH = 'System:AppData/Dendrogram/sample.nwk';

interface ImportNewickWitness {
  callReturnedArray: boolean;
  callReturnLen: number | null;
  tvPath: string | null;
  tableName: string | null;
  viewerTypes: string[];
  rowCount: number;
  colNames: string[];
  tagNewick: string | null;
  tagNewickJsonNonEmpty: boolean;
  leafRows: string[];
  canvasWidth: number;
  canvasHeight: number;
  canvasBoundingW: number;
  canvasBoundingH: number;
  dendrogramHostCount: number;
}

async function runImportNewick(page: Page, newick: string): Promise<ImportNewickWitness> {
  return await page.evaluate(async (nwk: string) => {
    const ret: any = await grok.functions.call('Dendrogram:importNewick', {fileContent: nwk});
    // Wait for DendrogramApp.buildView (addTableView + dock) to settle.
    await new Promise(r => setTimeout(r, 3500));
    const tv: any = grok.shell.tv;
    const df: any = tv?.dataFrame;
    const canvas = document.querySelector('[name="viewer-Dendrogram"] canvas') as HTMLCanvasElement | null;
    const br = canvas ? canvas.getBoundingClientRect() : {width: 0, height: 0};
    const viewerTypes = tv ? Array.from(tv.viewers || []).map((v: any) => v.type) : [];
    const colNames = df ? df.columns.names() : [];
    const leafRows: string[] = [];
    if (df && df.col('leaf') && df.col('node')) {
      const leafCol = df.col('leaf');
      const nodeCol = df.col('node');
      for (let i = 0; i < df.rowCount; i++)
        if (leafCol.get(i)) leafRows.push(String(nodeCol.get(i)));
    }
    return {
      callReturnedArray: Array.isArray(ret),
      callReturnLen: Array.isArray(ret) ? ret.length : null,
      tvPath: tv?.path ?? null,
      tableName: df?.name ?? null,
      viewerTypes,
      rowCount: df?.rowCount ?? 0,
      colNames,
      tagNewick: df?.getTag?.('.newick') ?? null,
      tagNewickJsonNonEmpty: !!(df?.getTag?.('.newickJson')?.length),
      leafRows: leafRows.sort(),
      canvasWidth: canvas?.width ?? 0,
      canvasHeight: canvas?.height ?? 0,
      canvasBoundingW: Math.round(br.width),
      canvasBoundingH: Math.round(br.height),
      dendrogramHostCount: document.querySelectorAll('[name="viewer-Dendrogram"]').length,
    };
  }, newick);
}

interface PreviewNewickWitness {
  callOk: boolean;
  callErr: string | null;
  returnedViewType: string | null;
  returnedRootClass: string | null;
  // After grok.shell.addView(view) the active view root carries the canvas:
  activeViewRootCanvases: number;
  activeViewCanvasBoundingW: number;
  activeViewCanvasBoundingH: number;
  activeViewRootClass: string | null;
}

async function runPreviewNewick(page: Page, filePath: string): Promise<PreviewNewickWitness> {
  return await page.evaluate(async (path: string) => {
    // Resolve a real FileInfo via dapi.files.list — pass-by-path strings
    // server-side fail with FileSystemException; the FileInfo from list()
    // carries the backing fileShare reference needed by file.readAsString().
    const folder = path.replace(/\/[^/]+$/, '');
    const baseName = path.split('/').pop() ?? '';
    const arr = await grok.dapi.files.list(folder);
    const fi = arr.find((f: any) => (f.fileName || f.name) === baseName);
    if (!fi) throw new Error(`Fixture not found in ${folder}: ${baseName}`);
    let callOk = false;
    let callErr: string | null = null;
    let ret: any = null;
    try {
      ret = await grok.functions.call('Dendrogram:previewNewick', {file: fi});
      callOk = true;
    } catch (e) {
      callErr = String(e).slice(0, 400);
    }
    let returnedRootClass: string | null = null;
    let returnedViewType: string | null = null;
    if (ret) {
      returnedViewType = ret.type ?? null;
      returnedRootClass = ret.root?.className ?? null;
      // Mount the returned DG.View into the shell so the canvas lays out.
      try { grok.shell.addView(ret); } catch (e) { callErr = (callErr || '') + ' | addView: ' + String(e).slice(0,200); }
    }
    await new Promise(r => setTimeout(r, 2500));
    const v: any = grok.shell.v;
    const root = v?.root as HTMLElement | undefined;
    const canvases = root ? root.querySelectorAll('canvas') : ({length: 0} as any);
    let firstBr = {width: 0, height: 0};
    if (canvases.length > 0) {
      const c = canvases[0] as HTMLCanvasElement;
      firstBr = c.getBoundingClientRect();
    }
    return {
      callOk, callErr,
      returnedViewType, returnedRootClass,
      activeViewRootCanvases: canvases.length,
      activeViewCanvasBoundingW: Math.round(firstBr.width),
      activeViewCanvasBoundingH: Math.round(firstBr.height),
      activeViewRootClass: root?.className ?? null,
    };
  }, filePath);
}

async function leafSetViaTreeHelper(page: Page, newick: string): Promise<string[]> {
  return await page.evaluate(async (nwk: string) => {
    const helper: any = await grok.functions.call('Dendrogram:getTreeHelper');
    const df: any = helper.newickToDf(nwk, 'leaves_probe');
    const leafCol = df.col('leaf');
    const nodeCol = df.col('node');
    const leaves: string[] = [];
    for (let i = 0; i < df.rowCount; i++)
      if (leafCol.get(i)) leaves.push(String(nodeCol.get(i)));
    return leaves.sort();
  }, newick);
}

test('Dendrogram: Open .nwk via importNewick (DendrogramApp) + preview via previewNewick (PhylocanvasGL)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // ── Setup ───────────────────────────────────────────────────────────────
  // Write the synthetic 4-leaf newick fixture under
  // System:AppData/Dendrogram/sample.nwk so that previewNewick has a
  // server-resident file to read via file.readAsString(). Setup runs outside
  // any softStep — Setup failures MUST surface immediately, not be swallowed.
  await page.evaluate(async (args: {path: string; newick: string}) => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 800));
    // Clean up any leftover from a prior aborted run, then write.
    try { await grok.dapi.files.delete(args.path); } catch (e) { /* tolerate */ }
    await grok.dapi.files.writeAsText(args.path, args.newick);
    const exists = await grok.dapi.files.exists(args.path);
    if (!exists) throw new Error(`Fixture not written: ${args.path}`);
  }, {path: FIXTURE_PATH, newick: NEWICK});

  // ── Scenario 1 — importNewick opens DendrogramApp on the parsed tree ────

  await softStep('1.1 Fixture sample.nwk is present in System:AppData/Dendrogram', async () => {
    const present = await page.evaluate(async (path: string) => {
      const folder = path.replace(/\/[^/]+$/, '');
      const baseName = path.split('/').pop() ?? '';
      const arr = await grok.dapi.files.list(folder);
      return arr.some((f: any) => (f.fileName || f.name) === baseName);
    }, FIXTURE_PATH);
    expect(present, 'sample.nwk visible in Dendrogram AppData listing').toBe(true);
  });

  await softStep('1.2 Dendrogram:importNewick runs without throwing; DendrogramApp mounts a TableView with viewer-Dendrogram', async () => {
    const w = await runImportNewick(page, NEWICK);
    // The fileHandler return contract is Promise<DG.DataFrame[]> and the
    // handler returns [] — verify the call did not throw AND the side-effect
    // (DendrogramApp.init -> addTableView + dock viewer-Dendrogram) ran.
    expect(w.callReturnedArray, 'importNewick returned an array').toBe(true);
    expect(w.callReturnLen, 'importNewick return length').toBe(0);
    // E-LAYER-COMPLIANCE-01: the ≥1 DOM-driving call for this spec.
    await page.locator('[name="viewer-Dendrogram"]').waitFor({timeout: 30_000});
    await page.locator('[name="viewer-Dendrogram"] canvas').waitFor({timeout: 30_000});
    expect(w.dendrogramHostCount, '[name="viewer-Dendrogram"] mounted exactly once').toBe(1);
    // DendrogramApp tv.path witness (apps/dendrogram-app.ts#L84).
    expect(w.tvPath, 'tv.path set by DendrogramApp.buildView').toBe('/func/Dendrogram.dendrogramApp');
    expect(w.viewerTypes, 'tv.viewers includes Grid + Dendrogram').toEqual(
      expect.arrayContaining(['Grid', 'Dendrogram']),
    );
  });

  await softStep('1.3 The active view\'s DataFrame is the parsed tree DataFrame (7 rows, [node,parent,leaf,distance])', async () => {
    const w = await runImportNewick(page, NEWICK);
    // 4 leaves + 3 internal nodes for a balanced binary 4-leaf newick.
    expect(w.rowCount, 'parsed tree DataFrame row count').toBe(7);
    expect(w.colNames, 'parsed tree DataFrame columns').toEqual(
      expect.arrayContaining(['node', 'parent', 'leaf', 'distance']),
    );
  });

  await softStep('1.4 The DataFrame is tagged with .newick (literal) and .newickJson; treeNewick round-trip', async () => {
    const w = await runImportNewick(page, NEWICK);
    // Per scenario Step 4: the .newick tag value equals the file content
    // exactly. The .newickJson tag is populated as the parsed JSON form.
    expect(w.tagNewick, '.newick tag equals literal file content').toBe(NEWICK);
    expect(w.tagNewickJsonNonEmpty, '.newickJson tag is populated').toBe(true);
  });

  await softStep('1.5 Leaf set = {a, b, c, d} via the leaf column on the parsed DataFrame', async () => {
    const w = await runImportNewick(page, NEWICK);
    expect(w.leafRows, 'leaf rows where leaf==true').toEqual([...EXPECTED_LEAVES]);
  });

  await softStep('1.6 Dendrogram viewer canvas has non-zero bounding rect (renderer mounted)', async () => {
    const w = await runImportNewick(page, NEWICK);
    expect(w.canvasBoundingW, 'canvas bounding rect width > 0').toBeGreaterThan(0);
    expect(w.canvasBoundingH, 'canvas bounding rect height > 0').toBeGreaterThan(0);
  });

  // ── Scenario 2 — previewNewick renders a PhylocanvasGL preview view ─────

  await softStep('2.1 Close the importNewick view; ensure fixture is still on disk before preview', async () => {
    await page.evaluate(async (path: string) => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 800));
      // Idempotency: scenario .md Step 2 of Scenario 2 calls out
      // re-writing the fixture at the start in case cleanup ran between
      // scenarios. The setup already wrote it and Scenario 1 did not delete
      // it, but we re-write defensively.
      const exists = await grok.dapi.files.exists(path);
      if (!exists)
        await grok.dapi.files.writeAsText(path, '((a:1,b:1):1,(c:1,d:1):1);');
    }, FIXTURE_PATH);
  });

  await softStep('2.2 Dendrogram:previewNewick returns a DG.View with a single PhylocanvasGL canvas', async () => {
    const w = await runPreviewNewick(page, FIXTURE_PATH);
    expect(w.callOk, `previewNewick call ok (err: ${w.callErr ?? ''})`).toBe(true);
    expect(w.returnedViewType, 'previewNewick returns a DG.View').toBe('view');
    // The returned root MAY have an empty or "ui-panel" className depending
    // on whether the DG.View has been hosted yet; class-2 recon-noted —
    // assert against the post-addView active-view root.
    expect(w.activeViewRootClass, 'active view root class').toContain('grok-view');
    // A single PhylocanvasGL canvas is mounted inside the view root.
    expect(w.activeViewRootCanvases, 'active view root canvases').toBe(1);
    expect(w.activeViewCanvasBoundingW, 'preview canvas bounding rect width > 0').toBeGreaterThan(0);
    expect(w.activeViewCanvasBoundingH, 'preview canvas bounding rect height > 0').toBeGreaterThan(0);
  });

  await softStep('2.3 Preview parse contract: TreeHelper.newickToDf produces the same leaf set as Scenario 1', async () => {
    // The PhylocanvasGL preview view does not expose a [name="viewer-*"]
    // anchor for its underlying DataFrame; per scenario Step 3 the
    // assertion is via the canonical parse contract — feed the same bytes
    // through TreeHelper.newickToDf and verify the leaf-set invariant.
    // This is the same code path the previewNewick handler walks
    // (package.ts#L308-L310).
    const leaves = await leafSetViaTreeHelper(page, NEWICK);
    expect(leaves, 'TreeHelper.newickToDf leaf set on preview bytes').toEqual([...EXPECTED_LEAVES]);
  });

  await softStep('2.4 Cleanup: delete the fixture; AppData/Dendrogram no longer contains sample.nwk', async () => {
    const after = await page.evaluate(async (path: string) => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      let deleteErr: string | null = null;
      try { await grok.dapi.files.delete(path); } catch (e) { deleteErr = String(e).slice(0, 200); }
      const folder = path.replace(/\/[^/]+$/, '');
      const baseName = path.split('/').pop() ?? '';
      const arr = await grok.dapi.files.list(folder);
      const stillThere = arr.some((f: any) => (f.fileName || f.name) === baseName);
      return {deleteErr, stillThere};
    }, FIXTURE_PATH);
    expect(after.deleteErr, `files.delete error: ${after.deleteErr ?? ''}`).toBeNull();
    expect(after.stillThere, 'sample.nwk removed from AppData/Dendrogram').toBe(false);
  });

  // ── Final cleanup ───────────────────────────────────────────────────────
  // Defensive: if scenario-2 cleanup softStep failed mid-way (leaving the
  // fixture on disk), best-effort delete here keeps the AppData folder
  // idempotent. Swallow any error — the assertion lives in 2.4.
  await page.evaluate(async (path: string) => {
    try { await grok.dapi.files.delete(path); } catch (e) { /* tolerate */ }
    try { grok.shell.closeAll(); } catch (e) {}
  }, FIXTURE_PATH);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
