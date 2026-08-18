import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer — JS API extension (viewPdbById / viewPdbByData)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => pageErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);

  await page.evaluate(() => {
    const g = (window as any).grok;
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    g.shell.closeAll();
    document.body.classList.add('selenium');
    g.shell.windows.simpleMode = true;
  });

  try {

    await softStep('Pre-flight — DG.Func.find verifies viewPdbById + viewPdbByData registration shape', async () => {
      const res = await page.evaluate(() => {
        const g = (window as any).grok;
        const D = (window as any).DG;
        const mapInputs = (fn: any) => (fn && fn.inputs)
          ? fn.inputs.map((i: any) => ({
              name: i.name,
              type: i.propertyType,
              optional: i.options?.optional ?? false,
            }))
          : [];
        const byIdFns = D.Func.find({name: 'viewPdbById', package: 'BiostructureViewer'});
        const byDataFns = D.Func.find({name: 'viewPdbByData', package: 'BiostructureViewer'});
        const byIdFn = byIdFns && byIdFns[0];
        const byDataFn = byDataFns && byDataFns[0];
        return {
          byIdRegistered: !!byIdFn,
          byIdInputs: mapInputs(byIdFn),
          byDataRegistered: !!byDataFn,
          byDataInputs: mapInputs(byDataFn),
        };
      });

      expect(res.byIdRegistered).toBe(true);
      expect(res.byIdInputs).toEqual([
        {name: 'pdbId', type: 'string', optional: false},
      ]);

      expect(res.byDataRegistered).toBe(true);
      expect(res.byDataInputs).toEqual([
        {name: 'pdbData', type: 'string', optional: false},
        {name: 'name', type: 'string', optional: false},
      ]);
    });

    await softStep('Scenario 1 — viewPdbById opens PDB by ID; shell state + DOM host mount', async () => {

      consoleErrors.length = 0;
      pageErrors.length = 0;

      const res = await page.evaluate(async () => {
        const g = (window as any).grok;
        g.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));

        const preViewCount = g.shell.views ? Array.from(g.shell.views).length : 0;

        const t0 = Date.now();
        let callResolved = false;
        let callRejectMessage: string | null = null;
        try {
          await g.functions.call('BiostructureViewer:viewPdbById', {pdbId: '1QBS'});
          callResolved = true;
        } catch (e: any) {
          callRejectMessage = String(e?.message ?? e);
        }
        const callDurationMs = Date.now() - t0;

        await new Promise((r) => setTimeout(r, 3000));

        const postViewNames: Array<{name: string, type: string}> = [];
        if (g.shell.views) {
          for (const v of g.shell.views) postViewNames.push({name: v.name, type: v.type});
        }
        const curr = g.shell.v;
        const currInfo = curr ? {name: curr.name, type: curr.type} : null;

        const hasMspPlugin = !!document.querySelector('.msp-plugin');
        const hasMspViewport = !!document.querySelector('.msp-viewport');
        const hasMspCanvas = !!document.querySelector('.msp-viewport canvas');

        return {
          preViewCount,
          callResolved,
          callRejectMessage,
          callDurationMs,
          postViewNames,
          postViewCount: postViewNames.length,
          currInfo,
          hasMspPlugin,
          hasMspViewport,
          hasMspCanvas,
        };
      });

      const acceptedReject = res.callRejectMessage === 'timeout';
      expect(res.callResolved || acceptedReject).toBe(true);

      expect(res.postViewCount).toBeGreaterThan(res.preViewCount);

      expect(res.currInfo).not.toBe(null);
      expect(res.currInfo!.name).toBe('Mol*');

      expect(res.currInfo!.type).toBe('view');

      expect(res.hasMspPlugin).toBe(true);

      const pitfallRegex = /Parsed object is empty/i;
      const pitfallConsole = consoleErrors.filter((m) => pitfallRegex.test(m));
      const pitfallPage = pageErrors.filter((m) => pitfallRegex.test(m));
      expect(pitfallConsole.length).toBe(0);
      expect(pitfallPage.length).toBe(0);
    });

    await softStep('Scenario 2 — viewPdbByData opens PDB from raw string + name; shell state + DOM host mount', async () => {
      consoleErrors.length = 0;
      pageErrors.length = 0;

      const res = await page.evaluate(async (pdbPath) => {
        const g = (window as any).grok;
        g.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));

        let pdbContent: string | null = null;
        let pdbReadError: string | null = null;
        try {
          pdbContent = await g.dapi.files.readAsText(pdbPath);
        } catch (e: any) {
          pdbReadError = String(e?.message ?? e);
        }

        const preViewCount = g.shell.views ? Array.from(g.shell.views).length : 0;

        const t0 = Date.now();
        let callResolved = false;
        let callRejectMessage: string | null = null;
        try {
          await g.functions.call(
            'BiostructureViewer:viewPdbByData',
            {pdbData: pdbContent, name: '1QBS'},
          );
          callResolved = true;
        } catch (e: any) {
          callRejectMessage = String(e?.message ?? e);
        }
        const callDurationMs = Date.now() - t0;

        await new Promise((r) => setTimeout(r, 3000));

        const postViewNames: Array<{name: string, type: string}> = [];
        if (g.shell.views) {
          for (const v of g.shell.views) postViewNames.push({name: v.name, type: v.type});
        }
        const curr = g.shell.v;
        const currInfo = curr ? {name: curr.name, type: curr.type} : null;

        const hasMspPlugin = !!document.querySelector('.msp-plugin');
        const hasMspViewport = !!document.querySelector('.msp-viewport');
        const hasMspCanvas = !!document.querySelector('.msp-viewport canvas');

        return {
          pdbReadError,
          pdbContentLen: pdbContent ? pdbContent.length : 0,
          preViewCount,
          callResolved,
          callRejectMessage,
          callDurationMs,
          postViewNames,
          postViewCount: postViewNames.length,
          currInfo,
          hasMspPlugin,
          hasMspViewport,
          hasMspCanvas,
        };
      }, samplePdbPath);

      expect(res.pdbReadError).toBe(null);
      expect(res.pdbContentLen).toBeGreaterThan(1000);

      const acceptedReject = res.callRejectMessage === 'timeout';
      expect(res.callResolved || acceptedReject).toBe(true);

      expect(res.postViewCount).toBeGreaterThan(res.preViewCount);
      expect(res.currInfo).not.toBe(null);

      expect(res.currInfo!.name).toBe('1QBS');
      expect(res.currInfo!.type).toBe('view');

      expect(res.hasMspPlugin).toBe(true);

      const pitfallRegex = /Parsed object is empty/i;
      const pitfallConsole = consoleErrors.filter((m) => pitfallRegex.test(m));
      const pitfallPage = pageErrors.filter((m) => pitfallRegex.test(m));
      expect(pitfallConsole.length).toBe(0);
      expect(pitfallPage.length).toBe(0);
    });
  } finally {

    await page.evaluate(() => {
      const g = (window as any).grok;
      g.shell.closeAll();
    });
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
