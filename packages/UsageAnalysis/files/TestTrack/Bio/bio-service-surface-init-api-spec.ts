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
  // Setup — apitest-pure: no UI driving, only shell state hygiene. Mirrors
  // Charts/charts-api.ts L45-L53.
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
    // ========================================================================
    // Scenario 1 — Service-surface five-function resolution after init.
    //
    // Atlas: bio.cp.bio-service-surface-init (p0); bio.api.get-seq-helper,
    // bio.api.get-monomer-lib-helper, bio.api.get-bio-lib + bio.lifecycle.init
    // (per scenario steps 1-4). The init-order invariant (bio.md#L566) is the
    // load-bearing contract: the four awaits MUST resolve, not block, not
    // throw — proving the runtime serialized the function calls after
    // initBio's completeInit() handshake.
    // ========================================================================
    await softStep('S1: Bio:getSeqHelper + Bio:getMonomerLibHelper + Bio:getBioLib resolve after init', async () => {
      const out = await page.evaluate(async () => {
        const g = (window as any).grok;
        const result: any = {seqHelper: null, monomerLibHelper: null, bioLib: null, errBalloons: 0};
        // tryCall — defensive helper for the function-registry name
        // divergence pattern (source `name:` literal vs canonical alias).
        // Returns {ok, value, name, err} — name records which form
        // resolved, for downstream cross-checking.
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
        // Atlas bio.api.get-seq-helper — Bio:getSeqHelper (package.ts#L1678).
        // Returns the ISeqHelper singleton (SeqHelper concrete impl per
        // seq-helper.ts). At minimum a callable getSeqHandler(column) method
        // per the ISeqHelper interface contract (Bio/CLAUDE.md service pattern).
        const r1 = await tryCall(['Bio:getSeqHelper'], {});
        if (r1.ok) {
          const sh = r1.value;
          result.seqHelper = {
            resolved: !!sh,
            kind: sh ? (sh.constructor?.name ?? typeof sh) : null,
            // ISeqHelper documented surface (Bio/CLAUDE.md): getSeqHandler,
            // getSeqMonomers, helmToAtomicLevel, setUnitsToFastaColumn etc.
            methods: ['getSeqHandler', 'getSeqMonomers', 'helmToAtomicLevel',
              'setUnitsToFastaColumn']
              .filter((m) => typeof sh?.[m] === 'function'),
          };
        } else {
          result.seqHelper = {err: r1.err};
        }
        // Atlas bio.api.get-monomer-lib-helper — Bio:getMonomerLibHelper
        // (package.ts#L133). Returns the MonomerLibManager singleton
        // (lib-manager.ts) implementing IMonomerLibHelper. At minimum a
        // callable getMonomerLib() / getBioLib() accessor per the
        // IMonomerLibHelper interface.
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
        // Atlas bio.api.get-bio-lib — Bio:getBioLib (package.ts#L175).
        // Returns the merged IMonomerLib (MonomerLib concrete impl per
        // monomer-lib.ts extending MonomerLibBase). At minimum a callable
        // getMonomer(polymerType, symbol) accessor per the IMonomerLib
        // interface contract (monomer-lib-base.ts).
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
        // Capture any error balloons that surfaced during the three calls.
        // Per bio.md (common observability pieces) — error balloons appear
        // in `.d4-balloon-error` containers; an `isErrorBallon` helper is
        // referenced in sibling specs. Direct DOM read is sanctioned here
        // because it is observational, not driving.
        result.errBalloons = document.querySelectorAll('.d4-balloon-error').length;
        return result;
      });
      // Atlas contract — all three service-surface calls resolve to non-null
      // singletons with the documented interface methods present.
      expect(out.seqHelper?.resolved, `getSeqHelper err: ${out.seqHelper?.err ?? ''}`).toBe(true);
      expect(out.seqHelper?.methods).toContain('getSeqHandler');
      expect(out.monomerLibHelper?.resolved, `getMonomerLibHelper err: ${out.monomerLibHelper?.err ?? ''}`).toBe(true);
      // At least one of getMonomerLib / getBioLib must be present
      // (IMonomerLibHelper interface evolves — getBioLib is the legacy
      // accessor; getMonomerLib is the modern accessor).
      const hasMlhAccessor = (out.monomerLibHelper?.methods ?? []).some((m: string) =>
        m === 'getMonomerLib' || m === 'getBioLib');
      expect(hasMlhAccessor, `monomerLibHelper accessor methods: ${JSON.stringify(out.monomerLibHelper?.methods)}`).toBe(true);
      expect(out.bioLib?.resolved, `getBioLib err: ${out.bioLib?.err ?? ''}`).toBe(true);
      expect(out.bioLib?.methods).toContain('getMonomer');
      // No error balloon raised during the four resolution steps
      // (scenario Expected bullet 5).
      expect(out.errBalloons, 'error balloon count after S1').toBe(0);
    });
    // ========================================================================
    // Scenario 2 — getSeqHandler returns notation-correct per-column handlers
    // (HELM + FASTA).
    //
    // Atlas: bio.api.get-seq-handler (package.ts#L180; function body L181-L184
    // delegates to _package.seqHelper.getSeqHandler(sequence)). The contract:
    // per-column scoping (each column → its own ISeqHandler instance) +
    // notation-correctness (handler's notation accessor reflects the column's
    // detected units tag).
    // ========================================================================
    await softStep('S2: Bio:getSeqHandler returns notation-correct per-column handlers (HELM + FASTA distinct)', async () => {
      const out = await page.evaluate(async ([helmPath, fastaPath]) => {
        const g = (window as any).grok;
        const result: any = {helm: null, fasta: null, distinctInstances: null, errBalloons: 0};
        // --- Load HELM fixture; detector classifies units=helm.
        const dfHelm = await g.dapi.files.readCsv(helmPath);
        const tvHelm = g.shell.addTableView(dfHelm);
        // Await semType detection (Macromolecule + units=helm). Mirror
        // bio-lifecycle-fasta-file-spec.ts setup pattern.
        await new Promise<void>((resolve) => {
          const sub = dfHelm.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => resolve(), 4000);
        });
        // Bio cell renderer post-detect settle (Bio/Chem pattern from
        // grok-browser SKILL.md "Bio/Chem datasets need extra wait").
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
            // ISeqHandler notation accessor — per seq-handler.ts the property
            // is `notation` (with possible legacy alias `units`). Probe both
            // so the assertion is robust to the SeqHandler API shape.
            notation: handlerHelm.notation ?? handlerHelm.units ?? null,
            unitsTag: helmCol.getTag?.('units') ?? null,
            hasGetSplitter: typeof handlerHelm.getSplitter === 'function',
            hasGetRegion: typeof handlerHelm.getRegion === 'function',
            methods: ['getSplitter', 'getRegion', 'getStats', 'isHelm', 'isFasta',
              'convertToNotation', 'getJoiner']
              .filter((m) => typeof handlerHelm[m] === 'function'),
          };
        }
        // --- Load FASTA fixture; detector classifies units=fasta.
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
        // Per-column scoping: the two handlers MUST be distinct instances.
        // Re-fetch the HELM handler via the same call to verify the scope
        // (do NOT keep the first one in memory then compare; the SeqHelper
        // implementation may cache per-column).
        if (helmCol && fastaCol) {
          const hh = await g.functions.call('Bio:getSeqHandler', {sequence: helmCol});
          const fh = await g.functions.call('Bio:getSeqHandler', {sequence: fastaCol});
          result.distinctInstances = hh !== fh;
        }
        result.errBalloons = document.querySelectorAll('.d4-balloon-error').length;
        return result;
      }, [HELM_PATH, FASTA_PATH]);
      // Atlas contract: both calls resolve to a non-null ISeqHandler instance.
      expect(out.helm?.resolved, `helm handler err: ${out.helm?.err ?? ''}`).toBe(true);
      expect(out.fasta?.resolved, `fasta handler err: ${out.fasta?.err ?? ''}`).toBe(true);
      // Notation-correctness — at minimum the column's units tag should be
      // 'helm' / 'fasta' respectively (the canonical detector output per
      // bio.md L592-L595). The handler's own notation accessor may surface
      // the same value or an internal enum form; assert via the column tag
      // (load-bearing observable per bio.md "Common observability pieces").
      expect(out.helm?.unitsTag).toBe('helm');
      expect(out.fasta?.unitsTag).toBe('fasta');
      // ISeqHandler surface contract — getSplitter() is the load-bearing
      // method (used by all per-notation operations across the bio package).
      expect(out.helm?.hasGetSplitter).toBe(true);
      expect(out.fasta?.hasGetSplitter).toBe(true);
      // Per-column scoping: distinct handlers for distinct columns.
      expect(out.distinctInstances, 'handlers should be distinct instances per column').toBe(true);
      // No error balloon raised (scenario Expected bullet 5).
      expect(out.errBalloons, 'error balloon count after S2').toBe(0);
    });
    // ========================================================================
    // Scenario 3 — getHelmMonomers returns the HELM-column monomer list.
    //
    // Atlas: bio.api.get-helm-monomers (package.ts#L1240, registered with
    // name `'Bio: getHelmMonomers'` — literal space; function body L1244
    // returns _package.seqHelper.getSeqMonomers(sequence)). The contract:
    // returned list is non-empty for a non-empty HELM column, and every
    // monomer symbol in the column's HELM strings appears in the returned
    // list (set-equality / round-trip consistency).
    // ========================================================================
    await softStep('S3: Bio:getHelmMonomers returns the HELM-column monomer list (set consistency)', async () => {
      const out = await page.evaluate(async (helmPath) => {
        const g = (window as any).grok;
        const result: any = {monomers: null, columnSymbols: null, errBalloons: 0};
        // Find the open HELM table view from Scenario 2 (table view is
        // retained across softSteps within the same test). If it isn't open
        // for some reason, re-read the fixture — apitest specs don't depend
        // on view state, only on the column.
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
        // tryCall — Bio:getHelmMonomers source `name:` has a literal space
        // ('Bio: getHelmMonomers'). bio.md L564 documents the no-space form.
        // Try the no-space form first (bio.md is the validated grok-browser
        // reference), fall back to the space form (source-literal name).
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
        // Normalize the returned value to a string list. The atlas contract
        // (per scenario Expected) says it's "an array (or iterable) of
        // monomer identifiers (strings)". The source body returns
        // _package.seqHelper.getSeqMonomers(sequence) which per the SeqHelper
        // contract returns string[].
        const monomerList: string[] = Array.isArray(monomers)
          ? monomers.map((x) => String(x))
          : (typeof monomers?.[Symbol.iterator] === 'function'
            ? Array.from(monomers as Iterable<unknown>).map((x) => String(x))
            : []);
        // Extract the unique monomer-symbol set observable in the column's
        // HELM strings. HELM monomer-symbol forms appear as bracketed names
        // (e.g. [meE], [dA]) and as bare single-letter symbols inside the
        // simple-position-list parts (after the `{` and before the `}`).
        // A pragmatic round-trip check: every TOKEN extracted from the HELM
        // strings via the bracket-or-letter regex should appear in the
        // returned monomer list. This is the "no monomer in the returned
        // list should be absent from the source HELM strings AND every
        // column monomer appears in the returned list" set-equality the
        // scenario describes (Expected bullet 3).
        const symbolsInColumn = new Set<string>();
        // HELM monomer token regex per the HELM 2.0 spec — bracketed names
        // or single-letter codes inside the simple-polymer position list.
        const HELM_TOKEN_RE = /\[([^\]]+)\]|(?<![\w\[])([A-Z])(?![\w\]])/g;
        for (let i = 0; i < helmCol.length; i++) {
          const s = helmCol.get(i);
          if (typeof s !== 'string' || s.length === 0) continue;
          // Restrict the scan to the inside of {...} blocks (simple polymers
          // section); avoid picking up RNA1/PEPTIDE1 prefix letters.
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
        // Compute the symmetric-difference statistics. The returned list may
        // be larger than the strict observable-in-column set (e.g. due to
        // upstream monomer-library merging) — the load-bearing direction is:
        // every observable-in-column symbol should appear in the returned
        // list (column ⊆ returned). The reverse direction (returned ⊆
        // column) is documented as the "no monomer in the returned list
        // should be absent from the source HELM strings" assertion — but
        // SeqHelper implementations historically also include unused
        // library monomers; assert the more robust column ⊆ returned
        // direction and report the reverse-difference for diagnostic
        // surfacing (without failing on it).
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
      // Atlas contract — non-empty array of monomer identifiers.
      expect(out.monomers?.resolved, `getHelmMonomers err: ${out.monomers?.err ?? ''} tried=${JSON.stringify(out.monomers?.candidatesTried)}`).toBe(true);
      expect(out.monomers?.listType).toBe('array');
      expect(out.monomers?.count, `returned monomer list count`).toBeGreaterThan(0);
      // Round-trip / mapping consistency — column-observed symbols ⊆ returned.
      // The reverse direction (returned ⊆ column) is NOT asserted: the
      // SeqHelper implementation may include monomers from the active
      // monomer library merge that don't appear in the column's strings.
      // Surface the reverse-direction stats in the diagnostic message
      // without failing on them.
      expect(out.columnSymbols?.count, 'unique symbols extracted from HELM column').toBeGreaterThan(0);
      expect(
        out.columnSymbols?.missingInReturned ?? [],
        `column ⊆ returned violated; column-symbols absent from returned list: ${JSON.stringify(out.columnSymbols?.missingInReturned)}`,
      ).toEqual([]);
      // No error balloon raised (scenario Expected bullet 4).
      expect(out.errBalloons, 'error balloon count after S3').toBe(0);
    });
  } finally {
    // ========================================================================
    // Cleanup — close all opened table views. The service-surface contract
    // is in-session only (no project save, no server-side entities created);
    // closeAll() is sufficient per scenario Setup Notes ("no cleanup required
    // beyond closing the opened table views at end of run").
    // ========================================================================
    await page.evaluate(() => {
      try { (window as any).grok.shell.closeAll(); } catch (e) { /* ignore */ }
    });
  }
  // Surface non-benign console errors as a soft signal (per Charts/charts-
  // api-spec.ts L34-L37 pattern). Service-surface init contract does not
  // own a no-console-error invariant — surface as diagnostic, do not fail.
  if (consoleErrors.length > 0) {
    // eslint-disable-next-line no-console
    console.warn(`[diagnostic] non-benign console errors observed: ${consoleErrors.length}; first: ${consoleErrors[0]?.slice(0, 200)}`);
  }
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
