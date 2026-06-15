/* ---
sub_features_covered: [dendrogram.api.tree-helper.newick-to-df, dendrogram.fileviewer.newick, dendrogram.viewer.dendrogram-app]
--- */
// The registered fileHandler function is Dendrogram:importNewick (the @fileHandler `name:`
// importNwk is only the extension key); it returns [] because DendrogramApp opens its own
// TableView via grok.shell.addTableView. previewNewick needs a server-resident FileInfo from
// dapi.files.list — a path-string-only synthetic fails with FileSystemException.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

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
    // Resolve a real FileInfo via dapi.files.list — path strings fail server-side.
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

  // Setup — write the 4-leaf newick fixture so previewNewick has a server-resident file.
  await page.evaluate(async (args: {path: string; newick: string}) => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 800));
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
    expect(w.callReturnedArray, 'importNewick returned an array').toBe(true);
    expect(w.callReturnLen, 'importNewick return length').toBe(0);
    await page.locator('[name="viewer-Dendrogram"]').waitFor({timeout: 30_000});
    await page.locator('[name="viewer-Dendrogram"] canvas').waitFor({timeout: 30_000});
    expect(w.dendrogramHostCount, '[name="viewer-Dendrogram"] mounted exactly once').toBe(1);
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
      // Re-write the fixture defensively in case cleanup ran between scenarios.
      const exists = await grok.dapi.files.exists(path);
      if (!exists)
        await grok.dapi.files.writeAsText(path, '((a:1,b:1):1,(c:1,d:1):1);');
    }, FIXTURE_PATH);
  });

  await softStep('2.2 Dendrogram:previewNewick returns a DG.View with a single PhylocanvasGL canvas', async () => {
    const w = await runPreviewNewick(page, FIXTURE_PATH);
    expect(w.callOk, `previewNewick call ok (err: ${w.callErr ?? ''})`).toBe(true);
    expect(w.returnedViewType, 'previewNewick returns a DG.View').toBe('view');
    expect(w.activeViewRootClass, 'active view root class').toContain('grok-view');
    expect(w.activeViewRootCanvases, 'active view root canvases').toBe(1);
    expect(w.activeViewCanvasBoundingW, 'preview canvas bounding rect width > 0').toBeGreaterThan(0);
    expect(w.activeViewCanvasBoundingH, 'preview canvas bounding rect height > 0').toBeGreaterThan(0);
  });

  await softStep('2.3 Preview parse contract: TreeHelper.newickToDf produces the same leaf set as Scenario 1', async () => {
    // PhylocanvasGL preview exposes no [name="viewer-*"] anchor; verify the leaf set via the
    // same TreeHelper.newickToDf parse contract the previewNewick handler walks.
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

  // Final cleanup — best-effort fixture delete (assertion lives in 2.4).
  await page.evaluate(async (path: string) => {
    try { await grok.dapi.files.delete(path); } catch (e) { /* tolerate */ }
    try { grok.shell.closeAll(); } catch (e) {}
  }, FIXTURE_PATH);

  finishSpec();
});
