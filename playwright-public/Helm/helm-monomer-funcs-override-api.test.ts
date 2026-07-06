import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

// Pattern 2 self-documented spec — paired scenario is
// helm-monomer-funcs-override-api.md (API-contract scenario; no UI flow drives it).
// Selector sources (grok-browser/references):
//   - .claude/skills/grok-browser/references/helm.md  (§ "Monomer-override contract (CHANGED)")

test.use(specTestOptions);

// 2026-06 editor rewrite: Pistoia removed. overrideMonomersFuncs /
// revertOriginalMonomersFuncs / buildMonomersFuncsFromLib delegate to an internal
// editorAdapter; there is NO window.org.helm.webeditor.Monomers to observe.
// Verify via the helper's own observable state (hh.originalMonomersFuncs) instead.
// Behaviour changes vs the old Pistoia build:
//   - double-apply override no longer throws (old "overridden already" guard removed)
//   - pre-revert no longer throws (old "Unable to revert" guard removed)
//   - buildMonomersFuncsFromLib no longer strips outer [...] brackets
//   - rewriteLibraries + Monomers.addOneMonomer/clear are gone (old S3 removed)
// See grok-browser references/helm.md § "Monomer-override contract (CHANGED)".
test('Helm / monomer-funcs swap path — overrideMonomersFuncs / revertOriginalMonomersFuncs / buildMonomersFuncsFromLib (contract)', async ({page}) => {
  // Pure JS-API/contract spec: no Web Editor dialog open. 180s is ample.
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
    // Clean baseline: revert any override left over from a prior run.
    try {
      if (hh?.originalMonomersFuncs != null) hh.revertOriginalMonomersFuncs();
    } catch (_) { /* ignore — best-effort baseline */ }
  });

  try {
    await softStep('S1: overrideMonomersFuncs / revert state-machine — override records, revert clears, double-apply + pre-revert no longer throw', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {
          hhReady: false,
          baselineNull: false,
          overrideReturnedNonNull: false, afterOverrideNonNull: false,
          doubleApplyThrew: false, doubleApplyMsg: '',
          revertReturnedNonNull: false, afterRevertNull: false,
          preRevertThrew: false, preRevertMsg: '',
          errBalloonsDelta: 0, err: null,
        };
        const errBefore = document.querySelectorAll('.d4-balloon.error').length;
        try {
          const hh = (window as any).__hh;
          if (!hh) throw new Error('Pre-flight prelude missing: window.__hh is undefined');
          result.hhReady = true;
          // Clean baseline
          if (hh.originalMonomersFuncs != null) hh.revertOriginalMonomersFuncs();
          result.baselineNull = hh.originalMonomersFuncs === null;

          const sentinel = {
            getMonomer: (_a: any, _name: any) => ({id: 'SENTINEL-1', m: 'sentinel-mol-1', n: 'sentinel'}),
            getMonomerSet: (_a: any) => ({id: 'SENTINEL-SET-1'}),
          };
          const previous = hh.overrideMonomersFuncs(sentinel);
          result.overrideReturnedNonNull = previous != null;
          result.afterOverrideNonNull = hh.originalMonomersFuncs != null;

          // Double-apply: must NOT throw in the new build.
          const sentinel2 = {
            getMonomer: (_a: any, _name: any) => ({id: 'SENTINEL-2'}),
            getMonomerSet: (_a: any) => null,
          };
          try {
            hh.overrideMonomersFuncs(sentinel2);
          } catch (e: any) {
            result.doubleApplyThrew = true;
            result.doubleApplyMsg = String(e?.message ?? e).slice(0, 300);
          }

          const outgoing = hh.revertOriginalMonomersFuncs();
          result.revertReturnedNonNull = outgoing != null;
          result.afterRevertNull = hh.originalMonomersFuncs === null;

          // Pre-revert: must NOT throw in the new build.
          try {
            hh.revertOriginalMonomersFuncs();
          } catch (e: any) {
            result.preRevertThrew = true;
            result.preRevertMsg = String(e?.message ?? e).slice(0, 300);
          }
        } catch (e: any) {
          result.err = String(e?.message ?? e).slice(0, 400);
        }
        // Ensure clean state for the next step.
        try {
          const hh = (window as any).__hh;
          if (hh?.originalMonomersFuncs != null) hh.revertOriginalMonomersFuncs();
        } catch (_) { /* ignore */ }
        const errAfter = document.querySelectorAll('.d4-balloon.error').length;
        result.errBalloonsDelta = errAfter - errBefore;
        return result;
      });
      expect(out.hhReady, 'window.__hh present from pre-flight').toBe(true);
      expect(out.baselineNull, `clean baseline: hh.originalMonomersFuncs === null; err: ${out.err ?? ''}`).toBe(true);
      expect(out.overrideReturnedNonNull, 'overrideMonomersFuncs returns the previous funcs (non-null)').toBe(true);
      expect(out.afterOverrideNonNull, 'after override, hh.originalMonomersFuncs is non-null (override recorded)').toBe(true);
      // Behaviour change: double-apply no longer throws.
      expect(out.doubleApplyThrew, `double-apply MUST NOT throw in the new build; observed msg: ${out.doubleApplyMsg}`).toBe(false);
      expect(out.revertReturnedNonNull, 'revertOriginalMonomersFuncs returns the outgoing funcs (non-null)').toBe(true);
      expect(out.afterRevertNull, 'after revert, hh.originalMonomersFuncs is back to null').toBe(true);
      // Behaviour change: pre-revert no longer throws.
      expect(out.preRevertThrew, `pre-revert (no override in effect) MUST NOT throw in the new build; observed msg: ${out.preRevertMsg}`).toBe(false);
      expect(out.errBalloonsDelta, 'new error balloon delta after S1').toBe(0);
    });

    await softStep('S2: buildMonomersFuncsFromLib — plain lookup matches lib; unknown → missing placeholder; bracketed NOT stripped; getMonomerSet no-throw', async () => {
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

          // Behaviour change: brackets are NOT stripped at this layer anymore;
          // '[meI]' resolves to the missing placeholder rather than the 'meI' record.
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
      // Behaviour change assertion: bracketed symbol is NOT stripped → missing placeholder.
      expect(out.bracketedIsPlaceholder, `bracketed '[meI]' returns {id:'[meI]', n:'missing'} (brackets NOT stripped in new build); observed id=${out.bracketedId}, n=${out.bracketedN}`).toBe(true);
      expect(out.monomerSetThrew, `getMonomerSet('PEPTIDE') must not throw; err: ${out.monomerSetErr}`).toBe(false);
      expect(out.errBalloonsDelta, 'new error balloon delta after S2').toBe(0);
    });

    // S3 (rewriteLibraries) REMOVED in the 2026-06 build: rewriteLibraries is not a
    // method on IHelmHelper and the org.helm.webeditor.Monomers.addOneMonomer/clear
    // primitives it drove no longer exist. The lib→editor sync is covered indirectly
    // by S2 (buildMonomersFuncsFromLib) and by the cross-feature interaction
    // helm-input-bio-monomer-lib. No automatable surface remains for a direct test.
  } finally {
    await page.evaluate(() => {
      try {
        const hh = (window as any).__hh;
        if (hh?.originalMonomersFuncs != null) hh.revertOriginalMonomersFuncs();
      } catch (_) { /* ignore — teardown must not throw */ }
      try { delete (window as any).__hh; } catch (_) { /* ignore */ }
      try { delete (window as any).__lib; } catch (_) { /* ignore */ }
      try { (window as any).grok.shell.closeAll(); } catch (_) { /* ignore */ }
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
