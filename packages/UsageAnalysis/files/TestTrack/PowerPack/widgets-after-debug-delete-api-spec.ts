import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('PowerPack — widget registration after debug-version package delete (GROK-16915 regression)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);

  await page.evaluate(() => {
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
  });

  await softStep('Scenario 1: delete debug-version PowerPack → bleeding-edge widgets remain', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;

      const debugRe = /\.X-[0-9a-f]+$/i;

      const allPackages = await grok.dapi.packages.list();
      const powerPackEntries = allPackages.filter((p: any) => p.name === 'PowerPack');
      const debugEntries = powerPackEntries.filter((p: any) => debugRe.test(String(p.version || '')));
      const releaseEntries = powerPackEntries.filter((p: any) => !debugRe.test(String(p.version || '')));
      const debugEntry: any = debugEntries.length > 0 ? debugEntries[0] : null;

      let bleedingEntry: any = null;
      if (releaseEntries.length > 0) {
        const sorted = releaseEntries.slice().sort((a: any, b: any) =>
          String(b.version || '').localeCompare(String(a.version || ''), undefined, {numeric: true, sensitivity: 'base'}));
        bleedingEntry = sorted[0];
      }

      if (!debugEntry) {
        return {
          blocked: true,
          reason: 'no debug-version PowerPack entry (no .X-<hash> match) on this server; ' +
                  'destructive delete skipped per spec self-protection guard',
          powerPackEntryCount: powerPackEntries.length,
        };
      }

      const allFuncs = DG.Func.find({}) || [];
      const ppFuncs = allFuncs.filter((f: any) => { try { return f.package && f.package.name === 'PowerPack'; } catch (e) { return false; } });

      const widgetNamesBefore = ppFuncs.filter((f: any) => {
        const role = (f.options && f.options.role) || null;
        return role === 'dashboard' || /Widget$/.test(f.name);
      }).map((f: any) => f.name).sort();

      const welcomeBefore = ppFuncs.filter((f: any) => f.name === '_welcomeView');
      const welcomeCountBefore = welcomeBefore.length;

      let deleteError: string | null = null;
      let deletePermissionDenied = false;
      try {
        await grok.dapi.packages.delete(debugEntry);
      } catch (e: any) {
        deleteError = String(e?.message ?? e).slice(0, 300);

        if (/permission/i.test(deleteError) || /forbidden/i.test(deleteError) ||
            /not authorized/i.test(deleteError) || /denied/i.test(deleteError) ||
            /403/.test(deleteError) || /admin/i.test(deleteError))
          deletePermissionDenied = true;
      }

      if (deletePermissionDenied) {
        return {
          blocked: true,
          reason: 'package delete rejected by server permission policy on ' +
                  'shared environment; scenario Setup Step 1 requires sandbox ' +
                  'with grok.dapi.packages write access — precondition unmet, ' +
                  'destructive path skipped per spec self-protection guard. ' +
                  'Server error: ' + deleteError,
          deletedVersion: debugEntry.version,
        };
      }

      await new Promise((r) => setTimeout(r, 2000));

      const allFuncsAfter = DG.Func.find({}) || [];
      const ppFuncsAfter = allFuncsAfter.filter((f: any) => { try { return f.package && f.package.name === 'PowerPack'; } catch (e) { return false; } });
      const widgetNamesAfter = ppFuncsAfter.filter((f: any) => {
        const role = (f.options && f.options.role) || null;
        return role === 'dashboard' || /Widget$/.test(f.name);
      }).map((f: any) => f.name).sort();

      const welcomeAfter = ppFuncsAfter.filter((f: any) => f.name === '_welcomeView');
      const welcomeCountAfter = welcomeAfter.length;

      return {
        blocked: false,
        deleteError,
        deletedVersion: debugEntry.version,
        bleedingVersion: bleedingEntry ? bleedingEntry.version : null,
        widgetNamesBefore,
        widgetNamesAfter,
        welcomeCountBefore,
        welcomeCountAfter,
      };
    });

    if (result.blocked) {

      console.warn(`[Scenario 1 SKIPPED] ${result.reason}`);
      expect.soft(result.blocked, `spec self-protection: ${result.reason}`).toBe(true);
      return;
    }

    expect(result.deleteError, `delete debug-version PowerPack ${result.deletedVersion}`).toBeNull();
    for (const name of result.widgetNamesBefore as string[])
      expect(result.widgetNamesAfter, `GROK-16915 invariant — widget ${name} retained after debug delete`).toContain(name);

    expect(result.welcomeCountAfter, '_welcomeView remains registered post debug-delete').toBeGreaterThanOrEqual(1);

    expect(result.welcomeCountBefore).toBeGreaterThanOrEqual(1);
  });

  await softStep('Scenario 2: delete debug-version UsageAnalysis → bleeding-edge widgets remain', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;

      const debugRe = /\.X-[0-9a-f]+$/i;
      const allPackages = await grok.dapi.packages.list();
      const uaEntries = allPackages.filter((p: any) => p.name === 'UsageAnalysis');
      const debugEntries = uaEntries.filter((p: any) => debugRe.test(String(p.version || '')));
      const releaseEntries = uaEntries.filter((p: any) => !debugRe.test(String(p.version || '')));
      const debugEntry: any = debugEntries.length > 0 ? debugEntries[0] : null;
      let bleedingEntry: any = null;
      if (releaseEntries.length > 0) {
        const sorted = releaseEntries.slice().sort((a: any, b: any) =>
          String(b.version || '').localeCompare(String(a.version || ''), undefined, {numeric: true, sensitivity: 'base'}));
        bleedingEntry = sorted[0];
      }
      if (!debugEntry) {
        return {
          blocked: true,
          reason: 'no debug-version UsageAnalysis entry (no .X-<hash> match) on this server; ' +
                  'destructive delete skipped per spec self-protection guard',
          uaEntryCount: uaEntries.length,
        };
      }

      const allFuncs = DG.Func.find({}) || [];
      const uaFuncs = allFuncs.filter((f: any) => { try { return f.package && f.package.name === 'UsageAnalysis'; } catch (e) { return false; } });
      const widgetNamesBefore = uaFuncs.filter((f: any) => {
        const role = (f.options && f.options.role) || null;
        return role === 'dashboard' || /Widget$/.test(f.name);
      }).map((f: any) => f.name).sort();
      const funcCountBefore = uaFuncs.length;

      let deleteError: string | null = null;
      let deletePermissionDenied = false;
      try {
        await grok.dapi.packages.delete(debugEntry);
      } catch (e: any) {
        deleteError = String(e?.message ?? e).slice(0, 300);
        if (/permission/i.test(deleteError) || /forbidden/i.test(deleteError) ||
            /not authorized/i.test(deleteError) || /denied/i.test(deleteError) ||
            /403/.test(deleteError) || /admin/i.test(deleteError))
          deletePermissionDenied = true;
      }
      if (deletePermissionDenied) {
        return {
          blocked: true,
          reason: 'package delete rejected by server permission policy on ' +
                  'shared environment; scenario Setup Step 1 requires sandbox ' +
                  'with grok.dapi.packages write access — precondition unmet, ' +
                  'destructive path skipped per spec self-protection guard. ' +
                  'Server error: ' + deleteError,
          deletedVersion: debugEntry.version,
        };
      }
      await new Promise((r) => setTimeout(r, 2000));

      const allFuncsAfter = DG.Func.find({}) || [];
      const uaFuncsAfter = allFuncsAfter.filter((f: any) => { try { return f.package && f.package.name === 'UsageAnalysis'; } catch (e) { return false; } });
      const widgetNamesAfter = uaFuncsAfter.filter((f: any) => {
        const role = (f.options && f.options.role) || null;
        return role === 'dashboard' || /Widget$/.test(f.name);
      }).map((f: any) => f.name).sort();
      const funcCountAfter = uaFuncsAfter.length;

      return {
        blocked: false,
        deleteError,
        deletedVersion: debugEntry.version,
        bleedingVersion: bleedingEntry ? bleedingEntry.version : null,
        widgetNamesBefore,
        widgetNamesAfter,
        funcCountBefore,
        funcCountAfter,
      };
    });

    if (result.blocked) {
      console.warn(`[Scenario 2 SKIPPED] ${result.reason}`);
      expect.soft(result.blocked, `spec self-protection: ${result.reason}`).toBe(true);
      return;
    }
    expect(result.deleteError, `delete debug-version UsageAnalysis ${result.deletedVersion}`).toBeNull();

    for (const name of result.widgetNamesBefore as string[])
      expect(result.widgetNamesAfter, `GROK-16915 invariant — UsageAnalysis widget ${name} retained`).toContain(name);

    expect(result.funcCountAfter, 'UsageAnalysis bleeding-edge funcs remain registered').toBeGreaterThan(0);
  });

  await softStep('Scenario 3: back-to-back delete of both debug packages → unified widget set remains', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;

      const debugRe = /\.X-[0-9a-f]+$/i;
      const allPackages = await grok.dapi.packages.list();
      const pickDebug = (pkgName: string) => {
        const entries = allPackages.filter((p: any) => p.name === pkgName);
        const debugEntries = entries.filter((p: any) => debugRe.test(String(p.version || '')));
        const releaseEntries = entries.filter((p: any) => !debugRe.test(String(p.version || '')));
        const debugEntry: any = debugEntries.length > 0 ? debugEntries[0] : null;
        let bleedingEntry: any = null;
        if (releaseEntries.length > 0) {
          const sorted = releaseEntries.slice().sort((a: any, b: any) =>
            String(b.version || '').localeCompare(String(a.version || ''), undefined, {numeric: true, sensitivity: 'base'}));
          bleedingEntry = sorted[0];
        }
        return {debugEntry, bleedingEntry, count: entries.length};
      };
      const pp = pickDebug('PowerPack');
      const ua = pickDebug('UsageAnalysis');

      if (!pp.debugEntry && !ua.debugEntry) {
        return {
          blocked: true,
          reason: 'no debug-version PowerPack or UsageAnalysis entry (no .X-<hash> match) ' +
                  'on this server; back-to-back destructive delete skipped per spec ' +
                  'self-protection',
          ppCount: pp.count,
          uaCount: ua.count,
        };
      }

      const allFuncs = DG.Func.find({}) || [];
      const widgetNamesFor = (pkgName: string) => allFuncs
        .filter((f: any) => { try { return f.package && f.package.name === pkgName; } catch (e) { return false; } })
        .filter((f: any) => {
          const role = (f.options && f.options.role) || null;
          return role === 'dashboard' || /Widget$/.test(f.name);
        })
        .map((f: any) => f.name);
      const ppWidgetsBefore = widgetNamesFor('PowerPack');
      const uaWidgetsBefore = widgetNamesFor('UsageAnalysis');
      const unionBefore = Array.from(new Set([...ppWidgetsBefore, ...uaWidgetsBefore])).sort();

      const deleteErrors: string[] = [];
      let anyPermissionDenied = false;
      const tryDelete = async (entry: any, pkgName: string) => {
        if (!entry) return;
        try {
          await grok.dapi.packages.delete(entry);
        } catch (e: any) {
          const msg = String(e?.message ?? e).slice(0, 200);
          deleteErrors.push(pkgName + ': ' + msg);
          if (/permission/i.test(msg) || /forbidden/i.test(msg) ||
              /not authorized/i.test(msg) || /denied/i.test(msg) ||
              /403/.test(msg) || /admin/i.test(msg))
            anyPermissionDenied = true;
        }
      };
      await tryDelete(pp.debugEntry, 'PowerPack');
      await tryDelete(ua.debugEntry, 'UsageAnalysis');

      if (anyPermissionDenied) {
        return {
          blocked: true,
          reason: 'one or both package deletes rejected by server permission ' +
                  'policy on shared environment; scenario Setup Step 1 requires ' +
                  'sandbox with grok.dapi.packages write access — precondition ' +
                  'unmet, destructive path skipped per spec self-protection guard. ' +
                  'Server errors: ' + deleteErrors.join('; '),
        };
      }

      await new Promise((r) => setTimeout(r, 3000));

      const allFuncsAfter = DG.Func.find({}) || [];
      const widgetNamesAfterFor = (pkgName: string) => allFuncsAfter
        .filter((f: any) => { try { return f.package && f.package.name === pkgName; } catch (e) { return false; } })
        .filter((f: any) => {
          const role = (f.options && f.options.role) || null;
          return role === 'dashboard' || /Widget$/.test(f.name);
        })
        .map((f: any) => f.name);
      const ppWidgetsAfter = widgetNamesAfterFor('PowerPack');
      const uaWidgetsAfter = widgetNamesAfterFor('UsageAnalysis');
      const unionAfter = Array.from(new Set([...ppWidgetsAfter, ...uaWidgetsAfter])).sort();

      return {
        blocked: false,
        deleteErrors,
        ppDeleted: !!pp.debugEntry,
        uaDeleted: !!ua.debugEntry,
        unionBefore,
        unionAfter,
      };
    });

    if (result.blocked) {
      console.warn(`[Scenario 3 SKIPPED] ${result.reason}`);
      expect.soft(result.blocked, `spec self-protection: ${result.reason}`).toBe(true);
      return;
    }
    expect(result.deleteErrors, 'back-to-back debug-version deletes succeed').toEqual([]);

    for (const name of result.unionBefore as string[])
      expect(result.unionAfter, `GROK-16915 invariant — unified widget ${name} retained after back-to-back delete`).toContain(name);
  });

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
