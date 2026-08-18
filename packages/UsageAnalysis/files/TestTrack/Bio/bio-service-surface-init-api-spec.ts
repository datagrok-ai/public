/* ---
  - bio.api.get-bio-lib
  - bio.api.get-seq-handler
  - bio.api.get-helm-monomers
  - bio.api.get-seq-helper
  - bio.api.get-monomer-lib-helper
  - bio.lifecycle.init
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const HELM_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';
const FASTA_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
test('Bio service-surface init — getSeqHelper / getMonomerLibHelper / getBioLib / getSeqHandler / getHelmMonomers', async ({page}) => {
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

    await softStep('S1: Bio:getSeqHelper + Bio:getMonomerLibHelper + Bio:getBioLib resolve after init', async () => {
      const out = await page.evaluate(async () => {
        const g = (window as any).grok;
        const result: any = {seqHelper: null, monomerLibHelper: null, bioLib: null, errBalloons: 0};

        async function tryCall(candidates: string[], params: Record<string, unknown>): Promise<{ok: boolean; value: any; name: string | null; err: string | null}> {
          let lastErr: string | null = null;
          for (const n of candidates) {
            try {
              const v = await g.functions.call(n, params);
              return {ok: true, value: v, name: n, err: null};
            } catch (e) {
              lastErr = String(e).slice(0, 300);
            }
          }
          return {ok: false, value: null, name: null, err: lastErr};
        }

        const r1 = await tryCall(['Bio:getSeqHelper'], {});
        if (r1.ok) {
          const sh = r1.value;
          result.seqHelper = {
            resolved: !!sh,
            kind: sh ? (sh.constructor?.name ?? typeof sh) : null,

            methods: ['getSeqHandler', 'getSeqMonomers', 'helmToAtomicLevel',
              'setUnitsToFastaColumn']
              .filter((m) => typeof sh?.[m] === 'function'),
          };
        } else {
          result.seqHelper = {err: r1.err};
        }

        const r2 = await tryCall(['Bio:getMonomerLibHelper'], {});
        if (r2.ok) {
          const mlh = r2.value;
          result.monomerLibHelper = {
            resolved: !!mlh,
            kind: mlh ? (mlh.constructor?.name ?? typeof mlh) : null,
            methods: ['getMonomerLib', 'getBioLib', 'loadMonomerLib', 'awaitLoaded',
              'getInstance']
              .filter((m) => typeof mlh?.[m] === 'function'),
          };
        } else {
          result.monomerLibHelper = {err: r2.err};
        }

        const r3 = await tryCall(['Bio:getBioLib'], {});
        if (r3.ok) {
          const bl = r3.value;
          result.bioLib = {
            resolved: !!bl,
            kind: bl ? (bl.constructor?.name ?? typeof bl) : null,
            methods: ['getMonomer', 'getMonomerSet', 'getMonomerNames',
              'getPolymerTypes', 'getMonomerSymbolsByType']
              .filter((m) => typeof bl?.[m] === 'function'),
          };
        } else {
          result.bioLib = {err: r3.err};
        }

        result.errBalloons = document.querySelectorAll('.d4-balloon-error').length;
        return result;
      });

      expect(out.seqHelper?.resolved, `getSeqHelper err: ${out.seqHelper?.err ?? ''}`).toBe(true);
      expect(out.seqHelper?.methods).toContain('getSeqHandler');
      expect(out.monomerLibHelper?.resolved, `getMonomerLibHelper err: ${out.monomerLibHelper?.err ?? ''}`).toBe(true);

      const hasMlhAccessor = (out.monomerLibHelper?.methods ?? []).some((m: string) =>
        m === 'getMonomerLib' || m === 'getBioLib');
      expect(hasMlhAccessor, `monomerLibHelper accessor methods: ${JSON.stringify(out.monomerLibHelper?.methods)}`).toBe(true);
      expect(out.bioLib?.resolved, `getBioLib err: ${out.bioLib?.err ?? ''}`).toBe(true);
      expect(out.bioLib?.methods).toContain('getMonomer');

      expect(out.errBalloons, 'error balloon count after S1').toBe(0);
    });

    await softStep('S2: Bio:getSeqHandler returns notation-correct per-column handlers (HELM + FASTA distinct)', async () => {
      const out = await page.evaluate(async ([helmPath, fastaPath]) => {
        const g = (window as any).grok;
        const result: any = {helm: null, fasta: null, distinctInstances: null, errBalloons: 0};

        const dfHelm = await g.dapi.files.readCsv(helmPath);
        const tvHelm = g.shell.addTableView(dfHelm);

        await new Promise<void>((resolve) => {
          const sub = dfHelm.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => resolve(), 4000);
        });

        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 3000));
        const helmCols = Array.from({length: dfHelm.columns.length},
          (_, i) => dfHelm.columns.byIndex(i));
        const helmCol: any = helmCols.find((c: any) => c.semType === 'Macromolecule');
        if (!helmCol) {
          result.helm = {err: 'no Macromolecule column in filter_HELM.csv'};
        } else {
          let handlerHelm: any = null;
          let invokeErr: string | null = null;
          try {
            handlerHelm = await g.functions.call('Bio:getSeqHandler', {sequence: helmCol});
          } catch (e) { invokeErr = String(e).slice(0, 300); }
          result.helm = handlerHelm == null ? {err: invokeErr} : {
            resolved: true,
            kind: handlerHelm.constructor?.name ?? typeof handlerHelm,

            notation: handlerHelm.notation ?? handlerHelm.units ?? null,
            unitsTag: helmCol.getTag?.('units') ?? null,
            hasGetSplitter: typeof handlerHelm.getSplitter === 'function',
            hasGetRegion: typeof handlerHelm.getRegion === 'function',
            methods: ['getSplitter', 'getRegion', 'getStats', 'isHelm', 'isFasta',
              'convertToNotation', 'getJoiner']
              .filter((m) => typeof handlerHelm[m] === 'function'),
          };
        }

        const dfFasta = await g.dapi.files.readCsv(fastaPath);
        g.shell.addTableView(dfFasta);
        await new Promise<void>((resolve) => {
          const sub = dfFasta.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => resolve(), 4000);
        });
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 3000));
        const fastaCols = Array.from({length: dfFasta.columns.length},
          (_, i) => dfFasta.columns.byIndex(i));
        const fastaCol: any = fastaCols.find((c: any) => c.semType === 'Macromolecule');
        if (!fastaCol) {
          result.fasta = {err: 'no Macromolecule column in filter_FASTA.csv'};
        } else {
          let handlerFasta: any = null;
          let invokeErr: string | null = null;
          try {
            handlerFasta = await g.functions.call('Bio:getSeqHandler', {sequence: fastaCol});
          } catch (e) { invokeErr = String(e).slice(0, 300); }
          result.fasta = handlerFasta == null ? {err: invokeErr} : {
            resolved: true,
            kind: handlerFasta.constructor?.name ?? typeof handlerFasta,
            notation: handlerFasta.notation ?? handlerFasta.units ?? null,
            unitsTag: fastaCol.getTag?.('units') ?? null,
            hasGetSplitter: typeof handlerFasta.getSplitter === 'function',
          };
        }

        if (helmCol && fastaCol) {
          const hh = await g.functions.call('Bio:getSeqHandler', {sequence: helmCol});
          const fh = await g.functions.call('Bio:getSeqHandler', {sequence: fastaCol});
          result.distinctInstances = hh !== fh;
        }
        result.errBalloons = document.querySelectorAll('.d4-balloon-error').length;
        return result;
      }, [HELM_PATH, FASTA_PATH]);

      expect(out.helm?.resolved, `helm handler err: ${out.helm?.err ?? ''}`).toBe(true);
      expect(out.fasta?.resolved, `fasta handler err: ${out.fasta?.err ?? ''}`).toBe(true);

      expect(out.helm?.unitsTag).toBe('helm');
      expect(out.fasta?.unitsTag).toBe('fasta');

      expect(out.helm?.hasGetSplitter).toBe(true);
      expect(out.fasta?.hasGetSplitter).toBe(true);

      expect(out.distinctInstances, 'handlers should be distinct instances per column').toBe(true);

      expect(out.errBalloons, 'error balloon count after S2').toBe(0);
    });

    await softStep('S3: Bio:getHelmMonomers returns the HELM-column monomer list (set consistency)', async () => {
      const out = await page.evaluate(async (helmPath) => {
        const g = (window as any).grok;
        const result: any = {monomers: null, columnSymbols: null, errBalloons: 0};

        let dfHelm: any = null;
        for (const tv of (g.shell.tableViews ?? [])) {
          const df = tv.dataFrame;
          if (df && df.name && /helm/i.test(df.name)) { dfHelm = df; break; }
        }
        if (!dfHelm) {
          dfHelm = await g.dapi.files.readCsv(helmPath);
          g.shell.addTableView(dfHelm);
          await new Promise<void>((resolve) => {
            const sub = dfHelm.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
            setTimeout(() => resolve(), 4000);
          });
          for (let i = 0; i < 50; i++) {
            if (document.querySelector('[name="viewer-Grid"] canvas')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
          await new Promise((r) => setTimeout(r, 3000));
        }
        const helmCols = Array.from({length: dfHelm.columns.length},
          (_, i) => dfHelm.columns.byIndex(i));
        const helmCol: any = helmCols.find((c: any) => c.semType === 'Macromolecule');
        if (!helmCol) {
          result.monomers = {err: 'no Macromolecule column in filter_HELM.csv'};
          return result;
        }

        const candidates = ['Bio:getHelmMonomers', 'Bio: getHelmMonomers'];
        let monomers: any = null;
        let invokeName: string | null = null;
        let lastErr: string | null = null;
        for (const n of candidates) {
          try {
            monomers = await g.functions.call(n, {sequence: helmCol});
            invokeName = n;
            break;
          } catch (e) { lastErr = String(e).slice(0, 300); }
        }
        if (monomers == null) {
          result.monomers = {err: lastErr, candidatesTried: candidates};
          return result;
        }

        const monomerList: string[] = Array.isArray(monomers)
          ? monomers.map((x) => String(x))
          : (typeof monomers?.[Symbol.iterator] === 'function'
            ? Array.from(monomers as Iterable<unknown>).map((x) => String(x))
            : []);

        const symbolsInColumn = new Set<string>();

        const HELM_TOKEN_RE = /\[([^\]]+)\]|(?<![\w\[])([A-Z])(?![\w\]])/g;
        for (let i = 0; i < helmCol.length; i++) {
          const s = helmCol.get(i);
          if (typeof s !== 'string' || s.length === 0) continue;

          const blockRe = /\{([^}]*)\}/g;
          let bm: RegExpExecArray | null;
          while ((bm = blockRe.exec(s)) !== null) {
            const inner = bm[1];
            HELM_TOKEN_RE.lastIndex = 0;
            let tm: RegExpExecArray | null;
            while ((tm = HELM_TOKEN_RE.exec(inner)) !== null) {
              const sym = tm[1] ?? tm[2];
              if (sym && sym.length > 0) symbolsInColumn.add(sym);
            }
          }
        }
        result.monomers = {
          resolved: true,
          invokeName,
          count: monomerList.length,
          listType: Array.isArray(monomers) ? 'array' : typeof monomers,
          sample: monomerList.slice(0, 10),
        };
        result.columnSymbols = {
          unique: Array.from(symbolsInColumn),
          count: symbolsInColumn.size,
        };

        const returnedSet = new Set(monomerList);
        const missingInReturned: string[] = [];
        for (const sym of symbolsInColumn) {
          if (!returnedSet.has(sym)) missingInReturned.push(sym);
        }
        const extraInReturned: string[] = [];
        for (const sym of returnedSet) {
          if (!symbolsInColumn.has(sym)) extraInReturned.push(sym);
        }
        result.columnSymbols.missingInReturned = missingInReturned;
        result.columnSymbols.extraInReturned = extraInReturned.slice(0, 20);
        result.errBalloons = document.querySelectorAll('.d4-balloon-error').length;
        return result;
      }, HELM_PATH);

      expect(out.monomers?.resolved, `getHelmMonomers err: ${out.monomers?.err ?? ''} tried=${JSON.stringify(out.monomers?.candidatesTried)}`).toBe(true);
      expect(out.monomers?.listType).toBe('array');
      expect(out.monomers?.count, `returned monomer list count`).toBeGreaterThan(0);

      expect(out.columnSymbols?.count, 'unique symbols extracted from HELM column').toBeGreaterThan(0);
      expect(
        out.columnSymbols?.missingInReturned ?? [],
        `column ⊆ returned violated; column-symbols absent from returned list: ${JSON.stringify(out.columnSymbols?.missingInReturned)}`,
      ).toEqual([]);

      expect(out.errBalloons, 'error balloon count after S3').toBe(0);
    });
  } finally {

    await page.evaluate(() => {
      try { (window as any).grok.shell.closeAll(); } catch (e) {  }
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
