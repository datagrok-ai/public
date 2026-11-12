import {expect} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Helm monomer-funcs path — buildMonomersFuncsFromLib (contract)', async ({page}) => {

  test.setTimeout(180_000);
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
    const g = (window as any).grok;
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    g.shell.closeAll();
    document.body.classList.add('selenium');
    g.shell.windows.simpleMode = true;
  });

  await page.evaluate(async () => {
    const g = (window as any).grok;
    const hh = await g.functions.call('Helm:getHelmHelper');
    const libHelper = await g.functions.call('Bio:getMonomerLibHelper');
    const lib = libHelper.getMonomerLib();
    (window as any).__hh = hh;
    (window as any).__lib = lib;
  });

  try {
    await softStep('S1: buildMonomersFuncsFromLib — plain lookup matches lib; unknown → missing placeholder; bracketed NOT stripped; getMonomerSet no-throw', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {
          funcsResolved: false, hasGetMonomer: false, hasGetMonomerSet: false,
          plainResolved: false, plainId: null, plainMatchesLib: false, libDirectId: null,
          unknownIsPlaceholder: false, unknownId: null, unknownN: null,
          bracketedIsPlaceholder: false, bracketedId: null, bracketedN: null,
          monomerSetThrew: false, monomerSetErr: '',
          errBalloonsDelta: 0, err: null,
        };
        const errBefore = document.querySelectorAll('.d4-balloon.error').length;
        try {
          const hh = (window as any).__hh;
          const lib = (window as any).__lib;
          if (!hh || !lib) throw new Error('Pre-flight prelude missing: window.__hh or window.__lib undefined');
          const funcs = hh.buildMonomersFuncsFromLib(lib);
          result.funcsResolved = funcs != null;
          result.hasGetMonomer = typeof funcs?.getMonomer === 'function';
          result.hasGetMonomerSet = typeof funcs?.getMonomerSet === 'function';

          const fromFuncs = funcs.getMonomer('PEPTIDE', 'A');
          const fromLib = lib.getWebEditorMonomer('PEPTIDE', 'A');
          result.plainResolved = fromFuncs != null;
          result.plainId = fromFuncs?.id ?? null;
          result.libDirectId = fromLib?.id ?? null;
          result.plainMatchesLib = fromFuncs != null && fromLib != null &&
            (fromFuncs === fromLib || fromFuncs.id === fromLib.id);

          const unknown = funcs.getMonomer('PEPTIDE', 'Xz_NotInLib');
          result.unknownId = unknown?.id ?? null;
          result.unknownN = unknown?.n ?? null;
          result.unknownIsPlaceholder = unknown != null && unknown.id === 'Xz_NotInLib' && unknown.n === 'missing';

          const bracketed = funcs.getMonomer('PEPTIDE', '[meI]');
          result.bracketedId = bracketed?.id ?? null;
          result.bracketedN = bracketed?.n ?? null;
          result.bracketedIsPlaceholder = bracketed != null && bracketed.id === '[meI]' && bracketed.n === 'missing';

          try {
            funcs.getMonomerSet('PEPTIDE');
          } catch (e: any) {
            result.monomerSetThrew = true;
            result.monomerSetErr = String(e?.message ?? e).slice(0, 200);
          }
        } catch (e: any) {
          result.err = String(e?.message ?? e).slice(0, 400);
        }
        const errAfter = document.querySelectorAll('.d4-balloon.error').length;
        result.errBalloonsDelta = errAfter - errBefore;
        return result;
      });
      expect(out.funcsResolved, `buildMonomersFuncsFromLib non-null; err: ${out.err ?? ''}`).toBe(true);
      expect(out.hasGetMonomer, 'funcs.getMonomer is a function').toBe(true);
      expect(out.hasGetMonomerSet, 'funcs.getMonomerSet is a function').toBe(true);
      expect(out.plainResolved, `plain getMonomer('PEPTIDE','A') resolves non-null`).toBe(true);
      expect(out.plainId, `plain getMonomer('PEPTIDE','A').id === 'A'`).toBe('A');
      expect(out.plainMatchesLib, `funcs.getMonomer id matches lib.getWebEditorMonomer id; funcs.id=${out.plainId}, lib.id=${out.libDirectId}`).toBe(true);
      expect(out.unknownIsPlaceholder, `unknown 'Xz_NotInLib' returns {id:'Xz_NotInLib', n:'missing'}; observed id=${out.unknownId}, n=${out.unknownN}`).toBe(true);

      expect(out.bracketedIsPlaceholder, `bracketed '[meI]' returns {id:'[meI]', n:'missing'} (brackets NOT stripped in new build); observed id=${out.bracketedId}, n=${out.bracketedN}`).toBe(true);
      expect(out.monomerSetThrew, `getMonomerSet('PEPTIDE') must not throw; err: ${out.monomerSetErr}`).toBe(false);
      expect(out.errBalloonsDelta, 'new error balloon delta after S1').toBe(0);
    });

  } finally {
    await page.evaluate(() => {
      try { delete (window as any).__hh; } catch (_) {  }
      try { delete (window as any).__lib; } catch (_) {  }
      try { (window as any).grok.shell.closeAll(); } catch (_) {  }
    });
  }

  if (consoleErrors.length > 0) {
    console.warn(`[diagnostic] non-benign console errors observed: ${consoleErrors.length}; first: ${consoleErrors[0]?.slice(0, 200)}`);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
