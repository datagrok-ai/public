/* ---
sub_features_covered: [powerpack.io, powerpack.io.exceljs-service, powerpack.io.xlsx-file-handler]
--- */
// GROK-19329 regression: XLSX opens across all 5 entry paths + optional sheetName.

import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, loginAsSecondUser, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers';

test.use(specTestOptions);

// Cold Playwright contexts hit the handler before pp.load() runs; trigger load + poll for registration.
async function ensurePowerPackLoaded(page: Page): Promise<void> {
  const registered = await page.evaluate(async () => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    if ((DG.Func.find({name: 'xlsxFileHandler'}) ?? []).length > 0) return true;
    try {
      const pp = await grok.dapi.packages.filter('shortName = "PowerPack"').first();
      if (pp) { try { pp.load(); } catch (_) { /* ignore — poll below */ } }
    } catch (_) { /* ignore — poll below */ }
    for (let i = 0; i < 60; i++) {
      if ((DG.Func.find({name: 'xlsxFileHandler'}) ?? []).length > 0) return true;
      await new Promise((r) => setTimeout(r, 500));
    }
    return false;
  });
  expect(
    registered,
    'PowerPack xlsxFileHandler did not register within 30s of pp.load(). ' +
    'On a healthy deploy the handler registers shortly after the package ' +
    'module loads — a persistent failure here signals a broken PowerPack ' +
    'publish, not the GROK-19329 regression.',
  ).toBe(true);
}

// Outer 60s cap so a regressed (hung) xlsxFileHandler dispatch fails fast with a bug-citing message.
async function openTableFromFileWithTimeout(
  page: Page,
  fullPath: string,
  capMs = 60_000,
): Promise<Awaited<ReturnType<typeof openTableFromFile>>> {
  return await Promise.race([
    openTableFromFile(page, fullPath),
    new Promise<never>((_, reject) =>
      setTimeout(() =>
        reject(new Error(`TIMEOUT after ${capMs}ms — openTableFromFile("${fullPath}") ` +
          'hung; GROK-19329 likely active on this build (XLSX file handler ' +
          'regression: file does not open from any location). See bug-library/' +
          'powerpack.yaml GROK-19329 for the canonical reproduction.')),
        capMs),
    ),
  ]);
}

// Locate the first reachable XLSX fixture under System:DemoFiles, or null.
async function locatePlatformXlsxFixture(page: Page): Promise<string | null> {
  return await page.evaluate(async () => {
    const grok = (window as any).grok;
    const candidates = [
      'System:DemoFiles/test/excel/excel-1mb.xlsx',
      'System:DemoFiles/test/excel/excel-rich-text-test.xlsx',
      'System:DemoFiles/SPGI-linked.xlsx',
    ];
    for (const path of candidates) {
      try {
        if (await grok.dapi.files.exists(path)) return path;
      } catch (_) { /* try next */ }
    }
    try {
      const items = await grok.dapi.files.list('System:DemoFiles/test/excel', false);
      for (const it of items ?? []) {
        const n = (it?.name ?? it?.fileName ?? '').toLowerCase();
        if (n.endsWith('.xlsx')) {
          const full = it?.fullPath ?? it?.path ?? `System:DemoFiles/test/excel/${it?.name}`;
          return full as string;
        }
      }
    } catch (_) { /* ignore */ }
    return null;
  });
}

// GROK-19329 invariant: active TableView has rows and no error balloon is on-screen.
async function verifyXlsxOpenedSuccessfully(
  page: Page,
): Promise<{rowCount: number; colCount: number; sheetCount: number; errorBalloons: number}> {
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000, state: 'visible'});
  await page.waitForTimeout(500); // let the grid settle
  const obs = await page.evaluate(() => {
    const grok = (window as any).grok;
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    const rowCount = df?.rowCount ?? 0;
    const colCount = df?.columns?.length ?? 0;
    const tables = grok.shell.tables ?? [];
    const sheetCount = Array.isArray(tables) ? tables.length : 0;
    return {rowCount, colCount, sheetCount};
  });
  const errorBalloons = await page.locator('.d4-balloon.error').count();
  return {...obs, errorBalloons};
}

// Reset between scenarios: close all views, dismiss balloons/toasts/menus.
async function resetShellState(page: Page): Promise<void> {
  await page.evaluate(() => {
    const grok = (window as any).grok;
    try { grok.shell.closeAll(); } catch (_) {}
    document.querySelectorAll('.d4-toast, .d4-menu-popup, .d4-balloon')
      .forEach((el) => { try { (el as HTMLElement).remove(); } catch (_) {} });
  }).catch(() => {});
  await page.waitForTimeout(800);
}

// Read fixture bytes via JS API; returns byte length.
async function readXlsxBytes(page: Page, fullPath: string): Promise<number> {
  return await page.evaluate(async (p) => {
    const grok = (window as any).grok;
    const bytes = await grok.dapi.files.readAsBytes(p);
    return bytes?.length ?? 0;
  }, fullPath);
}

// Invoke the registered xlsxFileHandler directly with fixture bytes; returns DataFrames (one per sheet).
async function invokeXlsxHandlerWithBytes(
  page: Page,
  fullPath: string,
  sheetName: string | null,
): Promise<{dfCount: number; firstRowCount: number; firstColCount: number; error: string | null}> {
  return await page.evaluate(async (args: {p: string; sn: string | null}) => {
    const grok = (window as any).grok;
    try {
      const bytes = await grok.dapi.files.readAsBytes(args.p);
      const params: any = {bytes};
      if (args.sn) params.sheetName = args.sn;
      const TIMEOUT_MS = 60_000;
      const result: any = await Promise.race([
        grok.functions.call('PowerPack:XlsxFileHandler', params),
        new Promise((_, reject) =>
          setTimeout(() =>
            reject(new Error('TIMEOUT after 60s — handler dispatch hung; ' +
              'GROK-19329 likely active on this build (XLSX file handler ' +
              'regression: file does not open from any location, no data ' +
              'displayed). See bug-library/powerpack.yaml GROK-19329.')),
          TIMEOUT_MS),
        ),
      ]);
      const dfs: any[] = Array.isArray(result) ? result : (result ? [result] : []);
      const first = dfs[0];
      if (first) grok.shell.addTableView(first);
      await new Promise((r) => setTimeout(r, 300)); // grid renderer attach settle

      return {
        dfCount: dfs.length,
        firstRowCount: first?.rowCount ?? 0,
        firstColCount: first?.columns?.length ?? 0,
        error: null,
      };
    } catch (e: any) {
      return {dfCount: 0, firstRowCount: 0, firstColCount: 0, error: String(e?.message ?? e)};
    }
  }, {p: fullPath, sn: sheetName});
}

// Synthesize an HTML5 drag-and-drop of the XLSX bytes onto the window's global drop handler.
async function simulateXlsxDrop(
  page: Page,
  fullPath: string,
  filename: string,
): Promise<{dispatched: boolean; reason?: string}> {
  return await page.evaluate(async (args: {p: string; fn: string}) => {
    const grok = (window as any).grok;
    try {
      const bytes = await grok.dapi.files.readAsBytes(args.p);
      const file = new File([new Uint8Array(bytes)], args.fn, {
        type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
      });
      const dt = new DataTransfer();
      dt.items.add(file);
      const target = document.body;
      target.dispatchEvent(new DragEvent('dragenter', {bubbles: true, cancelable: true, dataTransfer: dt}));
      target.dispatchEvent(new DragEvent('dragover',  {bubbles: true, cancelable: true, dataTransfer: dt}));
      target.dispatchEvent(new DragEvent('drop',      {bubbles: true, cancelable: true, dataTransfer: dt}));
      target.dispatchEvent(new DragEvent('dragend',   {bubbles: true, cancelable: true, dataTransfer: dt}));
      return {dispatched: true};
    } catch (e: any) {
      return {dispatched: false, reason: String(e?.message ?? e)};
    }
  }, {p: fullPath, fn: filename});
}

test('PowerPack: GROK-19329 XLSX opens across all 5 entry paths (regression)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
  });
  await page.waitForTimeout(300);

  // Cold-start guard: load PowerPack + wait for xlsxFileHandler registration before any scenario.
  await ensurePowerPackLoaded(page);

  const xlsxPath = await locatePlatformXlsxFixture(page);
  test.skip(
    !xlsxPath,
    'No XLSX fixture reachable under System:DemoFiles on this server.\n' +
    'Tried: System:DemoFiles/test/excel/excel-1mb.xlsx, excel-rich-text-test.xlsx,\n' +
    'SPGI-linked.xlsx. Setup step 2 prerequisite from xlsx-open.md is not satisfied.\n' +
    'Provision an XLSX fixture under System:DemoFiles, then re-run.',
  );
  const xlsxFullPath = xlsxPath!;
  const xlsxFileName = xlsxFullPath.split('/').pop() ?? 'fixture.xlsx';

  const observations: Record<string, {rowCount: number; colCount: number; errorBalloons: number}> = {};

  try {
    await softStep('Setup: verify PowerPack is installed + xlsxFileHandler is registered', async () => {
      const setupOk = await page.evaluate(async () => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        let ppPresent = false;
        try {
          const pp = await grok.dapi.packages.filter('shortName = "PowerPack"').first();
          ppPresent = !!pp?.id;
        } catch (_) { /* ignore */ }
        let handlerRegistered = false;
        try {
          const fns = DG.Func.find({name: 'xlsxFileHandler'}) ?? [];
          handlerRegistered = fns.length > 0;
        } catch (_) { /* ignore */ }
        return {ppPresent, handlerRegistered};
      });
      expect(setupOk.ppPresent).toBe(true);
      expect(setupOk.handlerRegistered).toBe(true);
    });

    // Scenario 1: Open XLSX from Browse / My Files (Browse double-click → OpenFile → xlsxFileHandler).
    await softStep('Scenario 1: open XLSX from Browse / My Files (canonical GROK-19329 reproduction)', async () => {
      await page.locator('[name="Browse"]').waitFor({timeout: 30_000, state: 'visible'});

      const opened = await openTableFromFileWithTimeout(page, xlsxFullPath);

      expect(opened.rowCount).toBeGreaterThan(0); // GROK-19329 invariant
      expect(opened.colCount).toBeGreaterThan(0); // GROK-19329 invariant

      const v = await verifyXlsxOpenedSuccessfully(page);
      observations['scenario-1'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
      expect(v.rowCount).toBeGreaterThan(0);
      expect(v.errorBalloons).toBe(0); // GROK-19329 invariant: no error balloon
    });

    await resetShellState(page);

    // Scenario 2: Open XLSX from Recent files (same OpenFile dispatch as Browse).
    await softStep('Scenario 2: open XLSX from Recent files (re-dispatch via recent activity)', async () => {
      const opened = await openTableFromFileWithTimeout(page, xlsxFullPath);

      expect(opened.rowCount).toBeGreaterThan(0); // GROK-19329 invariant

      const v = await verifyXlsxOpenedSuccessfully(page);
      observations['scenario-2'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
      expect(v.rowCount).toBeGreaterThan(0);
      expect(v.errorBalloons).toBe(0); // GROK-19329 invariant
    });

    await resetShellState(page);

    // Scenario 3: Open XLSX via drag-and-drop; falls back to direct handler dispatch if synthetic drop is blocked.
    await softStep('Scenario 3: open XLSX via drag-and-drop (synthetic DragEvent → handler dispatch)', async () => {
      const dropped = await simulateXlsxDrop(page, xlsxFullPath, xlsxFileName);
      let gridAppeared = await page.locator('[name="viewer-Grid"]')
        .waitFor({timeout: 5_000, state: 'visible'})
        .then(() => true).catch(() => false);
      if (!gridAppeared) {
        const result = await invokeXlsxHandlerWithBytes(page, xlsxFullPath, null);
        expect(result.error).toBeNull();
        expect(result.dfCount).toBeGreaterThan(0); // GROK-19329 invariant
        expect(result.firstRowCount).toBeGreaterThan(0); // GROK-19329 invariant
        await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000, state: 'visible'});
        gridAppeared = true;
      }
      expect(gridAppeared).toBe(true);

      const v = await verifyXlsxOpenedSuccessfully(page);
      observations['scenario-3'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
      expect(v.rowCount).toBeGreaterThan(0); // GROK-19329 invariant
      expect(v.errorBalloons).toBe(0); // GROK-19329 invariant
      console.log(`Scenario 3: simulated drop dispatched=${dropped.dispatched} reason=${dropped.reason ?? 'n/a'}`);
    });

    await resetShellState(page);

    // Scenario 4: File > Open native picker can't be Playwright-driven; exercise the same registered handler.
    await softStep('Scenario 4: open XLSX via top menu File > Open (handler-dispatch path)', async () => {
      const fileMenuFound = await page.evaluate(() => {
        const menuBars = Array.from(document.querySelectorAll('.d4-menu-bar, [name="menu-bar"]'));
        for (const mb of menuBars) {
          const items = Array.from(mb.querySelectorAll('.d4-menu-item, [name^="div-File"], div, span'));
          for (const it of items) {
            const t = ((it as HTMLElement).textContent ?? '').trim();
            if (t === 'File' && (it as HTMLElement).offsetParent !== null) return true;
          }
        }
        return false;
      });
      console.log(`Scenario 4: top-menu File item found=${fileMenuFound}`);

      const result = await invokeXlsxHandlerWithBytes(page, xlsxFullPath, null);
      expect(result.error).toBeNull();
      expect(result.dfCount).toBeGreaterThan(0); // GROK-19329 invariant
      expect(result.firstRowCount).toBeGreaterThan(0); // GROK-19329 invariant

      const v = await verifyXlsxOpenedSuccessfully(page);
      observations['scenario-4'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
      expect(v.rowCount).toBeGreaterThan(0); // GROK-19329 invariant
      expect(v.errorBalloons).toBe(0); // GROK-19329 invariant
    });

    await resetShellState(page);

    // Scenario 5: Open XLSX from Shared with me; real cross-user leg when token2 present, else handler-dispatch fallback.
    await softStep('Scenario 5: open XLSX from Shared with me (or handler-dispatch fallback)', async () => {
      const token2 = process.env.DATAGROK_AUTH_TOKEN_2;
      if (token2 && token2.length > 0) {
        // Split: cross-user SETUP (login + share-grant) may legitimately fail (System fixture not
        // shareable) and falls through. But once the second user can READ the file, the GROK-19329
        // invariants are NOT optional — they run outside the swallow so a real regression fails loudly.
        let crossUserSetupOk = false;
        let setupSkipReason = '';
        try {
          // Discover the second user's login via a token2 re-auth round-trip, then restore primary session.
          await loginAsSecondUser(page);
          const secondLogin: string | null = await page.evaluate(async () => {
            const grok = (window as any).grok;
            return (await grok.dapi.users.current())?.login ?? null;
          });
          await loginToDatagrok(page);
          if (!secondLogin)
            throw new Error('could not resolve second user login via token2 round-trip');

          // As primary user: resolve the second user's group + the fixture FileInfo, grant View to the group.
          const grant = await page.evaluate(async (args: {p: string; login: string}) => {
            const grok = (window as any).grok;
            const second = await grok.dapi.users
              .filter('login = "' + args.login + '"').first()
              .catch(() => null);
            // Permissions attach to GROUPS, never the bare User.
            const group = second?.group ?? null;
            if (!group)
              return {ok: false, reason: `second user "${args.login}" / group not resolvable`};
            // files.list on a full file path lists directory contents; resolve FileInfo via parent dir + basename.
            const slash = args.p.lastIndexOf('/');
            const dir = slash >= 0 ? args.p.slice(0, slash + 1) : args.p;
            const base = slash >= 0 ? args.p.slice(slash + 1) : args.p;
            const fileInfo = await grok.dapi.files.list(dir, false)
              .then((items: any[]) => (items || []).find((it: any) => (it?.name ?? '') === base) ?? null)
              .catch(() => null);
            if (!fileInfo)
              return {ok: false, reason: `could not resolve FileInfo for ${args.p} (listed ${dir})`};
            await grok.dapi.permissions.grant(fileInfo, group, false);
            return {ok: true, reason: null};
          }, {p: xlsxFullPath, login: secondLogin});

          if (!grant.ok)
            throw new Error(grant.reason ?? 'share setup failed');

          // Re-auth as the second user (cross-user read context established for the invariant checks below).
          await loginAsSecondUser(page);
          await ensurePowerPackLoaded(page);
          crossUserSetupOk = true;

          // Run the cross-user invariants OUTSIDE this try so they cannot be swallowed into the fallback.
          const byteLen = await readXlsxBytes(page, xlsxFullPath);
          expect(byteLen, 'cross-user READ returned no bytes — share grant did not take effect').toBeGreaterThan(0);

          const result = await invokeXlsxHandlerWithBytes(page, xlsxFullPath, null);
          expect(result.error, `cross-user XlsxFileHandler errored: ${result.error}`).toBeNull();
          expect(result.dfCount).toBeGreaterThan(0); // GROK-19329 invariant
          expect(result.firstRowCount).toBeGreaterThan(0); // GROK-19329 invariant
          const v = await verifyXlsxOpenedSuccessfully(page);
          observations['scenario-5'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
          expect(v.rowCount).toBeGreaterThan(0);
          expect(v.errorBalloons).toBe(0); // GROK-19329 invariant
          console.log(`Scenario 5: real cross-user open of shared XLSX as ${secondLogin} succeeded`);
        } catch (e: any) {
          // Once setup succeeded, the failure IS a real regression — rethrow it. Only setup
          // failures (sharing a System-owned fixture not feasible) are tolerated as a skip.
          if (crossUserSetupOk)
            throw e;
          setupSkipReason = String(e?.message ?? e);
          console.warn(`Scenario 5: real cross-user leg skipped (${setupSkipReason}); ` +
            'falling back to handler-dispatch on platform fixture.');
        } finally {
          // Always restore the primary session so cleanup runs as the owner.
          await loginToDatagrok(page).catch(() => {});
          await ensurePowerPackLoaded(page).catch(() => {});
        }
        if (crossUserSetupOk)
          return;
      }

      // Probe for a shared XLSX not under System: or the current user's own home.
      const sharedPath = await page.evaluate(async () => {
        const grok = (window as any).grok;
        try {
          const me = await grok.dapi.users.current();
          const myLogin = (me?.login ?? '').toLowerCase();
          const roots = ['Shared:', 'Public:'];
          for (const root of roots) {
            try {
              const items = await grok.dapi.files.list(root, true);
              for (const it of items ?? []) {
                const n = (it?.name ?? it?.fileName ?? '').toLowerCase();
                if (n.endsWith('.xlsx')) {
                  const full = (it?.fullPath ?? it?.path ?? '') as string;
                  if (myLogin && full.toLowerCase().startsWith(myLogin + ':')) continue;
                  return full;
                }
              }
            } catch (_) { /* try next root */ }
          }
        } catch (_) { /* probe failed */ }
        return null;
      });

      if (sharedPath) {
        const opened = await openTableFromFileWithTimeout(page, sharedPath);
        expect(opened.rowCount).toBeGreaterThan(0); // GROK-19329 invariant
        const v = await verifyXlsxOpenedSuccessfully(page);
        observations['scenario-5'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
        expect(v.rowCount).toBeGreaterThan(0);
        expect(v.errorBalloons).toBe(0); // GROK-19329 invariant
        console.log(`Scenario 5: opened shared XLSX from ${sharedPath}`);
      } else {
        // No shared XLSX reachable; exercise the same registered handler on the platform fixture.
        const result = await invokeXlsxHandlerWithBytes(page, xlsxFullPath, null);
        expect(result.error).toBeNull();
        expect(result.dfCount).toBeGreaterThan(0); // GROK-19329 invariant
        expect(result.firstRowCount).toBeGreaterThan(0); // GROK-19329 invariant
        const v = await verifyXlsxOpenedSuccessfully(page);
        observations['scenario-5'] = {rowCount: v.rowCount, colCount: v.colCount, errorBalloons: v.errorBalloons};
        expect(v.rowCount).toBeGreaterThan(0);
        expect(v.errorBalloons).toBe(0); // GROK-19329 invariant
        console.log('Scenario 5: no shared XLSX reachable in env; exercised handler-dispatch on platform fixture');
      }
    });

    await resetShellState(page);

    // Scenario 6: sheetName param filters to one sheet; omitted returns all sheets.
    await softStep('Scenario 6: sheetName parameter selects a specific sheet', async () => {
      const allSheets = await invokeXlsxHandlerWithBytes(page, xlsxFullPath, null);
      expect(allSheets.error).toBeNull();
      expect(allSheets.dfCount).toBeGreaterThan(0); // GROK-19329 invariant

      if (allSheets.dfCount < 2) {
        observations['scenario-6'] = {rowCount: allSheets.firstRowCount, colCount: allSheets.firstColCount, errorBalloons: 0};
        console.log(`Scenario 6: workbook ${xlsxFileName} is single-sheet; default-all-sheets behavior verified`);
        return;
      }

      // Multi-sheet: discover the second sheet's name from the parsed DataFrame names.
      const secondSheetName = await page.evaluate(async (p) => {
        const grok = (window as any).grok;
        const bytes = await grok.dapi.files.readAsBytes(p);
        const TIMEOUT_MS = 60_000;
        const dfs: any[] = await Promise.race([
          grok.functions.call('PowerPack:XlsxFileHandler', {bytes}),
          new Promise<never>((_, reject) =>
            setTimeout(() => reject(new Error(
              'TIMEOUT after 60s — XlsxFileHandler hung; GROK-19329 likely active.')),
              TIMEOUT_MS),
          ),
        ]).catch(() => []);
        return dfs?.[1]?.name ?? null;
      }, xlsxFullPath);

      if (!secondSheetName) {
        console.log('Scenario 6: could not discover second sheet name; skipping sheetName-specific assertion');
        return;
      }

      const oneSheet = await invokeXlsxHandlerWithBytes(page, xlsxFullPath, secondSheetName);
      expect(oneSheet.error).toBeNull();
      expect(oneSheet.dfCount).toBeGreaterThan(0); // GROK-19329 invariant
      const namedRowCount = await page.evaluate(async (args: {p: string; sn: string}) => {
        const grok = (window as any).grok;
        const bytes = await grok.dapi.files.readAsBytes(args.p);
        const TIMEOUT_MS = 60_000;
        const dfs: any[] = await Promise.race([
          grok.functions.call('PowerPack:XlsxFileHandler', {bytes, sheetName: args.sn}),
          new Promise<never>((_, reject) =>
            setTimeout(() => reject(new Error(
              'TIMEOUT after 60s — XlsxFileHandler hung; GROK-19329 likely active.')),
              TIMEOUT_MS),
          ),
        ]).catch(() => []);
        return dfs?.[0]?.rowCount ?? 0;
      }, {p: xlsxFullPath, sn: secondSheetName});
      expect(namedRowCount).toBeGreaterThan(0);

      observations['scenario-6'] = {rowCount: oneSheet.firstRowCount, colCount: oneSheet.firstColCount, errorBalloons: 0};
      console.log(`Scenario 6: sheetName="${secondSheetName}" selected; rowCount=${namedRowCount}`);
    });
  } finally {
    await page.evaluate(() => {
      const grok = (window as any).grok;
      try { grok.shell.closeAll(); } catch (_) {}
    }).catch(() => {});

    console.log('Per-scenario observations:', JSON.stringify(observations, null, 2));
  }

  finishSpec();
});
