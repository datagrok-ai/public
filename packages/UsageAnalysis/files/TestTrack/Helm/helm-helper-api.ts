import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const HELM_PEPTIDE = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$';
const HELM_PEPTIDE_EXPECTED_ATOMS = 9;
const HELM_RNA_SHORT = 'RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$';

test('Helm HelmHelper accessor + factories + hit-test API surface — Helm:getHelmHelper / createHelmInput / createHelmWebEditor / getHoveredAtom', async ({page}) => {
  test.setTimeout(360_000);
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

  try {
    await softStep('S1: Helm:getHelmHelper resolves to a singleton IHelmHelper with full method surface', async () => {
      const out = await page.evaluate(async () => {
        const g = (window as any).grok;
        const result: any = {hhResolved: false, isSingleton: false, kind: null, methods: [], errBalloons: 0, err: null};
        try {
          const hh = await g.functions.call('Helm:getHelmHelper');
          const hh2 = await g.functions.call('Helm:getHelmHelper');
          result.hhResolved = hh != null;
          result.isSingleton = hh === hh2;
          result.kind = hh ? (hh.constructor?.name ?? typeof hh) : null;
          const methodNames = [
            'parse', 'removeGaps', 'getMolfiles', 'getHoveredAtom',
            'createHelmInput', 'createHelmWebEditor', 'createWebEditorApp',
            'overrideMonomersFuncs', 'revertOriginalMonomersFuncs',
            'buildMonomersFuncsFromLib',
          ];
          result.methods = methodNames.filter((m) => typeof hh?.[m] === 'function');
          (window as any).__hh = hh;
        } catch (e) {
          result.err = String(e).slice(0, 400);
        }
        result.errBalloons = document.querySelectorAll('.d4-balloon.error').length;
        return result;
      });
      expect(out.hhResolved, `getHelmHelper err: ${out.err ?? ''}`).toBe(true);
      // Singleton invariant
      expect(out.isSingleton, 'Helm:getHelmHelper singleton invariant').toBe(true);
      expect(
        out.methods,
        `IHelmHelper method surface mismatch; methods present: ${JSON.stringify(out.methods)}`,
      ).toEqual([
        'parse', 'removeGaps', 'getMolfiles', 'getHoveredAtom',
        'createHelmInput', 'createHelmWebEditor', 'createWebEditorApp',
        'overrideMonomersFuncs', 'revertOriginalMonomersFuncs',
        'buildMonomersFuncsFromLib',
      ]);
      expect(out.errBalloons, 'error balloon count after S1').toBe(0);
    });

    await softStep('S2: HelmHelper.createHelmInput returns a HelmInputBase whose host is a real HTMLElement', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {inputResolved: false, kind: null, hostIsHtmlElement: false, hostAppendable: false, negativePath: 'NOT_TESTED', errBalloons: 0, err: null};
        try {
          const hh = (window as any).__hh;
          if (!hh) throw new Error('Scenario 1 prelude missing: window.__hh is undefined');
          const input = hh.createHelmInput('macromolecule', {});
          result.inputResolved = input != null;
          result.kind = input ? (input.constructor?.name ?? typeof input) : null;
          const host = (input as any)?.viewerHost ?? (input as any)?.root ?? null;
          result.hostIsHtmlElement = host instanceof HTMLElement;
          if (host instanceof HTMLElement) {
            const sandbox = document.createElement('div');
            try {
              sandbox.appendChild(host);
              result.hostAppendable = sandbox.firstChild === host;
            } finally {
              try { sandbox.removeChild(host); } catch (_) { /* ignore */ }
            }
          }
        } catch (e) {
          result.err = String(e).slice(0, 400);
        }
        result.negativePath = 'SKIPPED: _libHelper module-internal; rethrow path covered by helm-helper.ts#L54 unit test';
        result.errBalloons = document.querySelectorAll('.d4-balloon.error').length;
        return result;
      });
      expect(out.inputResolved, `createHelmInput err: ${out.err ?? ''}`).toBe(true);
      expect(out.hostIsHtmlElement, `host element is HTMLElement; kind: ${out.kind}`).toBe(true);
      expect(out.hostAppendable, 'host element is appendable to a sandbox div').toBe(true);
      console.warn(`[diagnostic S2 negative-path] ${out.negativePath}`);
      expect(out.errBalloons, 'error balloon count after S2').toBe(0);
    });

    await softStep('S3: HelmHelper.createHelmWebEditor returns a view-only IHelmWebEditor; host-by-reference; no-host fallback ok', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {weResolved: false, weKind: null, hostMatches: false, editorPresent: false, editorViewonly: null, weNoHostResolved: false, weNoHostHost: null, errBalloons: 0, errWithHost: null, errNoHost: null};
        try {
          const hh = (window as any).__hh;
          if (!hh) throw new Error('Scenario 1 prelude missing: window.__hh is undefined');
          const host = document.createElement('div');
          host.style.width = '400px';
          host.style.height = '300px';
          try {
            const we = hh.createHelmWebEditor(host, {});
            result.weResolved = we != null;
            result.weKind = we ? (we.constructor?.name ?? typeof we) : null;
            result.hostMatches = we?.host === host;
            result.editorPresent = we?.editor != null;
            const ed = we?.editor;
            result.editorViewonly = ed?.options?.viewonly ?? ed?.viewonly ?? null;
          } catch (e) {
            result.errWithHost = String(e).slice(0, 400);
          }
          try {
            const we2 = hh.createHelmWebEditor(undefined, {});
            result.weNoHostResolved = we2 != null;
            result.weNoHostHost = we2?.host instanceof HTMLElement ? 'HTMLElement' :
              (we2?.host == null ? 'null' : typeof we2.host);
          } catch (e) {
            result.errNoHost = String(e).slice(0, 400);
          }
        } catch (e) {
          result.errWithHost = String(e).slice(0, 400);
        }
        result.errBalloons = document.querySelectorAll('.d4-balloon.error').length;
        return result;
      });
      expect(out.weResolved, `createHelmWebEditor(host) err: ${out.errWithHost ?? ''}`).toBe(true);
      expect(out.hostMatches, `we.host !== passed host; kind: ${out.weKind}`).toBe(true);
      expect(out.editorPresent, 'we.editor is non-null').toBe(true);
      expect(
        out.editorViewonly,
        `we.editor viewonly flag (probed options.viewonly + direct); observed: ${out.editorViewonly}`,
      ).toBe(true);
      expect(out.weNoHostResolved, `createHelmWebEditor(undefined) err: ${out.errNoHost ?? ''}`).toBe(true);
      expect(
        out.weNoHostHost,
        `no-host wrapper's host accessor; expected HTMLElement; observed: ${out.weNoHostHost}`,
      ).toBe('HTMLElement');
      expect(out.errBalloons, 'error balloon count after S3').toBe(0);
    });

    await softStep('S4: HelmHelper.getHoveredAtom — exact-coord match, far-away null, height-arg sweep, immutability', async () => {
      const out = await page.evaluate(async (helmStr) => {
        const result: any = {parseOk: false, atomsCount: 0, hoveredAtomSameRef: false, farReturnedNull: 'NOT_TESTED', heightSweepOk: false, immutableAtoms: false, errBalloons: 0, err: null};
        try {
          const hh = (window as any).__hh;
          if (!hh) throw new Error('Scenario 1 prelude missing: window.__hh is undefined');
          const mol = hh.parse(helmStr);
          result.parseOk = mol != null;
          result.atomsCount = mol?.atoms?.length ?? 0;
          if (mol && Array.isArray(mol.atoms) && mol.atoms.length > 0) {
            const initialCount = mol.atoms.length;
            const a0 = mol.atoms[0];
            const p0 = a0.p;
            const ha = hh.getHoveredAtom(p0.x, p0.y, mol, 100);
            result.hoveredAtomSameRef = ha === a0;
            const far = hh.getHoveredAtom(-10000, -10000, mol, 100);
            result.farReturnedNull = (far == null) ? 'null' : 'NOT_NULL';
            let sweepOk = true;
            for (const h of [40, 80, 200]) {
              try { hh.getHoveredAtom(p0.x, p0.y, mol, h); } catch (e) { sweepOk = false; break; }
            }
            result.heightSweepOk = sweepOk;
            result.immutableAtoms = mol.atoms.length === initialCount;
          }
        } catch (e) {
          result.err = String(e).slice(0, 400);
        }
        result.errBalloons = document.querySelectorAll('.d4-balloon.error').length;
        return result;
      }, HELM_PEPTIDE);
      expect(out.parseOk, `parse() err: ${out.err ?? ''}`).toBe(true);
      expect(out.atomsCount, 'fixture atom count').toBe(HELM_PEPTIDE_EXPECTED_ATOMS);
      expect(out.hoveredAtomSameRef, 'getHoveredAtom at atom-0 coord returns atom-0 by reference').toBe(true);
      expect(out.farReturnedNull, 'getHoveredAtom(-10000,-10000) returns null').toBe('null');
      expect(out.heightSweepOk, 'getHoveredAtom does not throw across height sweep').toBe(true);
      expect(out.immutableAtoms, 'mol.atoms not mutated by getHoveredAtom').toBe(true);
      expect(out.errBalloons, 'error balloon count after S4').toBe(0);
    });

    await softStep('S4-supp: HelmHelper.parse on RNA fixture yields a non-empty HelmMol (smoke)', async () => {
      const out = await page.evaluate(async (helmStr) => {
        const result: any = {parseOk: false, atomsCount: 0, errBalloons: 0, err: null};
        try {
          const hh = (window as any).__hh;
          if (!hh) throw new Error('Scenario 1 prelude missing: window.__hh is undefined');
          const mol = hh.parse(helmStr);
          result.parseOk = mol != null;
          result.atomsCount = mol?.atoms?.length ?? 0;
        } catch (e) {
          result.err = String(e).slice(0, 400);
        }
        result.errBalloons = document.querySelectorAll('.d4-balloon.error').length;
        return result;
      }, HELM_RNA_SHORT);
      expect(out.parseOk, `parse(RNA) err: ${out.err ?? ''}`).toBe(true);
      expect(out.atomsCount, 'RNA fixture atom count > 0').toBeGreaterThan(0);
      expect(out.errBalloons, 'error balloon count after S4-supp').toBe(0);
    });
  } finally {
    await page.evaluate(() => {
      try { delete (window as any).__hh; } catch (_) { /* ignore */ }
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
