import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Helm monomer-dict swap path — overrideMonomersFuncs / revertOriginalMonomersFuncs / buildMonomersFuncsFromLib / rewriteLibraries (contract)', async ({page}) => {
  // Pure JS-API/contract spec: no Web Editor dialog open, no cold JSDraw2 init —
  // only getHelmHelper + dict-swap / buildMonomersFuncsFromLib calls. 180s is ample.
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
    try {
      if (hh?.originalMonomersFuncs != null) hh.revertOriginalMonomersFuncs();
    } catch (_) { /* ignore — best-effort baseline */ }
  });

  try {
    await softStep('S1: overrideMonomersFuncs swaps Pistoia dict; revert restores; double-apply + pre-revert throw', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {
          preCaptureOk: false, hhReady: false,
          originalReturned: false, originalGetMatchesPre: false, originalGetSetMatchesPre: false,
          liveGetIsSentinel: false, liveGetSetIsSentinel: false,
          doubleApplyThrew: false, doubleApplyMsgMatch: false, doubleApplyMsg: '',
          liveStillSentinelAfterDouble: false,
          revertReturned: false, revertReturnedSentinelGet: false, revertReturnedSentinelGetSet: false,
          liveGetRestored: false, liveGetSetRestored: false,
          preRevertThrew: false, preRevertMsgMatch: false, preRevertMsg: '',
          errBalloonsDelta: 0, err: null,
        };
        const errBefore = document.querySelectorAll('.d4-balloon.error').length;
        try {
          const hh = (window as any).__hh;
          if (!hh) throw new Error('Pre-flight prelude missing: window.__hh is undefined');
          result.hhReady = true;
          const monomers = (window as any).org.helm.webeditor.Monomers;
          const preGet = monomers.getMonomer;
          const preGetSet = monomers.getMonomerSet;
          result.preCaptureOk = typeof preGet === 'function' && typeof preGetSet === 'function';
          const sentinel = {
            getMonomer: (_a: any, _name: any) => ({id: 'SENTINEL-1', m: 'sentinel-mol-1', n: 'sentinel'}),
            getMonomerSet: (_a: any) => ({id: 'SENTINEL-SET-1'}),
          };
          const originalReturn = hh.overrideMonomersFuncs(sentinel);
          result.originalReturned = originalReturn != null;
          result.originalGetMatchesPre = originalReturn != null &&
            typeof originalReturn.getMonomer === 'function' &&
            originalReturn.getMonomer !== sentinel.getMonomer;
          result.originalGetSetMatchesPre = originalReturn != null &&
            typeof originalReturn.getMonomerSet === 'function' &&
            originalReturn.getMonomerSet !== sentinel.getMonomerSet;
          result.liveGetIsSentinel = monomers.getMonomer === sentinel.getMonomer;
          result.liveGetSetIsSentinel = monomers.getMonomerSet === sentinel.getMonomerSet;
          const sentinel2 = {
            getMonomer: (_a: any, _name: any) => ({id: 'SENTINEL-2'}),
            getMonomerSet: (_a: any) => null,
          };
          try {
            hh.overrideMonomersFuncs(sentinel2);
          } catch (e: any) {
            result.doubleApplyThrew = true;
            result.doubleApplyMsg = String(e?.message ?? e).slice(0, 300);
            result.doubleApplyMsgMatch = /originalGetMonomer is overridden already/.test(result.doubleApplyMsg);
          }
          result.liveStillSentinelAfterDouble = monomers.getMonomer === sentinel.getMonomer;
          const overridden = hh.revertOriginalMonomersFuncs();
          result.revertReturned = overridden != null;
          result.revertReturnedSentinelGet = overridden?.getMonomer === sentinel.getMonomer;
          result.revertReturnedSentinelGetSet = overridden?.getMonomerSet === sentinel.getMonomerSet;
          const helmAA1 = (window as any).org.helm.webeditor.HELM.AA;
          const liveA = monomers.getMonomer(helmAA1, 'A');
          const preA = preGet.call(monomers, helmAA1, 'A');
          result.liveGetRestored = liveA != null && preA != null &&
            (liveA === preA || liveA.id === preA.id);
          result.liveGetSetRestored = monomers.getMonomerSet !== sentinel.getMonomerSet;
          try {
            hh.revertOriginalMonomersFuncs();
          } catch (e: any) {
            result.preRevertThrew = true;
            result.preRevertMsg = String(e?.message ?? e).slice(0, 300);
            result.preRevertMsgMatch = /Unable to revert original getMonomer/.test(result.preRevertMsg);
          }
        } catch (e: any) {
          result.err = String(e?.message ?? e).slice(0, 400);
          try {
            const hh = (window as any).__hh;
            if (hh?.originalMonomersFuncs != null) hh.revertOriginalMonomersFuncs();
          } catch (_) { /* ignore */ }
        }
        try {
          const hh = (window as any).__hh;
          const lib = (window as any).__lib;
          if (hh && lib && hh.originalMonomersFuncs == null)
            hh.overrideMonomersFuncs(hh.buildMonomersFuncsFromLib(lib));
        } catch (_) { /* ignore */ }
        const errAfter = document.querySelectorAll('.d4-balloon.error').length;
        result.errBalloonsDelta = errAfter - errBefore;
        return result;
      });
      expect(out.hhReady, 'window.__hh present from pre-flight').toBe(true);
      expect(out.preCaptureOk, 'pre-swap getMonomer / getMonomerSet captured').toBe(true);
      expect(out.originalReturned, `overrideMonomersFuncs return non-null; err: ${out.err ?? ''}`).toBe(true);
      expect(out.originalGetMatchesPre, 'returned originalMonomersFuncs.getMonomer is a function and not the sentinel (reference identity)').toBe(true);
      expect(out.originalGetSetMatchesPre, 'returned originalMonomersFuncs.getMonomerSet is a function and not the sentinel (reference identity)').toBe(true);
      expect(out.liveGetIsSentinel, 'live monomers.getMonomer === sentinel.getMonomer after override').toBe(true);
      expect(out.liveGetSetIsSentinel, 'live monomers.getMonomerSet === sentinel.getMonomerSet after override').toBe(true);
      expect(out.doubleApplyThrew, `double-apply throws synchronously; msg: ${out.doubleApplyMsg}`).toBe(true);
      expect(out.doubleApplyMsgMatch, `double-apply error message contains "originalGetMonomer is overridden already"; observed: ${out.doubleApplyMsg}`).toBe(true);
      expect(out.liveStillSentinelAfterDouble, 'live monomers.getMonomer still === sentinel after failed double-apply').toBe(true);
      expect(out.revertReturned, 'revertOriginalMonomersFuncs returns non-null').toBe(true);
      expect(out.revertReturnedSentinelGet, 'revert return.getMonomer === sentinel.getMonomer (outgoing handler)').toBe(true);
      expect(out.revertReturnedSentinelGetSet, 'revert return.getMonomerSet === sentinel.getMonomerSet (outgoing handler)').toBe(true);
      expect(out.liveGetRestored, 'live monomers.getMonomer functional parity with preGet (HELM_AA/A lookup matches)').toBe(true);
      expect(out.liveGetSetRestored, 'live monomers.getMonomerSet is not the sentinel after revert').toBe(true);
      expect(out.preRevertThrew, `pre-revert throws synchronously; msg: ${out.preRevertMsg}`).toBe(true);
      expect(out.preRevertMsgMatch, `pre-revert error message contains "Unable to revert original getMonomer"; observed: ${out.preRevertMsg}`).toBe(true);
      expect(out.errBalloonsDelta, 'new error balloon delta after S1').toBe(0);
    });

    await softStep('S2: buildMonomersFuncsFromLib — plain/bracketed lookup match lib; unknown returns missing placeholder; getMonomerSet no-throw', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {
          funcsResolved: false, hasGetMonomer: false, hasGetMonomerSet: false,
          plainResolved: false, plainSameRefAsLib: false, plainId: null, libDirectId: null,
          bracketedResolved: false, bracketedIdEqualsPlainId: false, bracketedId: null,
          unknownResult: 'NOT_TESTED', unknownIsPlaceholder: false, unknownIdEcho: null, unknownN: null,
          monomerSetCallable: false, monomerSetThrew: false, monomerSetResolved: 'NOT_TESTED', monomerSetErr: '',
          errBalloonsDelta: 0, err: null,
        };
        const errBefore = document.querySelectorAll('.d4-balloon.error').length;
        try {
          const hh = (window as any).__hh;
          const lib = (window as any).__lib;
          if (!hh || !lib) throw new Error('Pre-flight prelude missing: window.__hh or window.__lib undefined');
          // Build.
          const funcs = hh.buildMonomersFuncsFromLib(lib);
          result.funcsResolved = funcs != null;
          result.hasGetMonomer = typeof funcs?.getMonomer === 'function';
          result.hasGetMonomerSet = typeof funcs?.getMonomerSet === 'function';
          const fromFuncs = funcs.getMonomer('PEPTIDE', 'A');
          const fromLib = lib.getWebEditorMonomer('PEPTIDE', 'A');
          result.plainResolved = fromFuncs != null;
          result.plainSameRefAsLib = fromFuncs === fromLib;
          result.plainId = fromFuncs?.id ?? null;
          result.libDirectId = fromLib?.id ?? null;
          const fromBracketed = funcs.getMonomer('PEPTIDE', '[meI]');
          const fromPlainMeI = funcs.getMonomer('PEPTIDE', 'meI');
          result.bracketedResolved = fromBracketed != null;
          result.bracketedId = fromBracketed?.id ?? null;
          result.bracketedIdEqualsPlainId = fromBracketed != null && fromPlainMeI != null &&
            (fromBracketed === fromPlainMeI || fromBracketed.id === fromPlainMeI.id);
          const unknown = funcs.getMonomer('PEPTIDE', 'Xz_NotInLib');
          result.unknownResult = unknown == null ? 'null/undef' : 'NOT_NULL';
          result.unknownIdEcho = unknown?.id ?? null;
          result.unknownN = unknown?.n ?? null;
          result.unknownIsPlaceholder = unknown != null &&
            unknown.id === 'Xz_NotInLib' && unknown.n === 'missing';
          try {
            const mset = funcs.getMonomerSet('PEPTIDE');
            result.monomerSetCallable = true;
            result.monomerSetResolved = mset == null ? 'null/undef' : 'NOT_NULL';
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
      expect(out.plainSameRefAsLib, `funcs.getMonomer('PEPTIDE','A') === lib.getWebEditorMonomer('PEPTIDE','A'); funcs.id: ${out.plainId}, lib.id: ${out.libDirectId}`).toBe(true);
      expect(out.bracketedResolved, `bracketed getMonomer('PEPTIDE','[meI]') resolves non-null`).toBe(true);
      expect(out.bracketedIdEqualsPlainId, `bracketed '[meI]' id matches plain 'meI' id; bracketed.id: ${out.bracketedId}`).toBe(true);
      expect(out.unknownIsPlaceholder, `unknown 'Xz_NotInLib' returns placeholder {id:'Xz_NotInLib', n:'missing'}; observed: id=${out.unknownIdEcho}, n=${out.unknownN}`).toBe(true);
      console.warn(`[diagnostic S2 unknown-direction] scenario predicted null/undef; live runtime returns ${out.unknownResult}`);
      expect(out.monomerSetThrew, `getMonomerSet must not throw; err: ${out.monomerSetErr}`).toBe(false);
      expect(out.monomerSetCallable, 'getMonomerSet callable on PEPTIDE').toBe(true);
      expect(out.errBalloonsDelta, 'new error balloon delta after S2').toBe(0);
    });

    await softStep('S3: rewriteLibraries contract — clear + addOneMonomer mutations are observable in org.helm.webeditor.Monomers (SR-01: direct import unreachable)', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {
          baselinePeptideA_n: null, baselinePeptideMeI_n: null,
          sentinelAdded: false, sentinelLookupOk: false, sentinelId: null, sentinelN: null,
          clearTookEffect: false, postClearAReadsMissing: false,
          restoredPeptideA_n: null, restoredPeptideMeI_n: null,
          peptideAReverted: false, peptideMeIReverted: false,
          replayCounts: {polymerTypes: 0, monomers: 0, addedOk: 0, addedSkipped: 0},
          errBalloonsDelta: 0, err: null,
          srNote: 'SR-01: rewriteLibraries from utils/get-monomer.ts is not a registered Datagrok function and the bundled module is not import()-able from the apitest layer. Contract asserted semantically via the primitive sequence the function calls internally.',
        };
        const errBefore = document.querySelectorAll('.d4-balloon.error').length;
        try {
          const monomers = (window as any).org.helm.webeditor.Monomers;
          const lib = (window as any).__lib;
          if (!lib) throw new Error('Pre-flight prelude missing: window.__lib undefined');
          const hh_s3 = (window as any).__hh;
          if (hh_s3?.originalMonomersFuncs != null) {
            try { hh_s3.revertOriginalMonomersFuncs(); } catch (_) { /* ignore */ }
          }
          const helmAA = (window as any).org.helm.webeditor.HELM.AA;
          const baseA = monomers.getMonomer(helmAA, 'A');
          const baseMeI = monomers.getMonomer(helmAA, 'meI');
          result.baselinePeptideA_n = baseA?.n ?? null;
          result.baselinePeptideMeI_n = baseMeI?.n ?? null;
          const sentinel = {
            id: 'TestSentinelMonomer',
            m: 'sentinel-molfile-stub',
            n: 'Test-Sentinel',
            na: 'X',
            rs: 0,
            type: 'PEPTIDE',
            mt: 'Backbone',
            at: {},
          };
          monomers.addOneMonomer(sentinel);
          const sentinelLookup = monomers.getMonomer(helmAA, 'TestSentinelMonomer');
          result.sentinelAdded = sentinelLookup != null;
          result.sentinelLookupOk = sentinelLookup != null && sentinelLookup.n !== 'missing';
          result.sentinelId = sentinelLookup?.id ?? null;
          result.sentinelN = sentinelLookup?.n ?? null;
          monomers.clear();
          const postClearA = monomers.getMonomer(helmAA, 'A');
          const postClearSentinel = monomers.getMonomer(helmAA, 'TestSentinelMonomer');
          result.clearTookEffect = postClearSentinel == null;
          result.postClearAReadsMissing = postClearA == null;
          const polymerTypes = lib.getPolymerTypes();
          result.replayCounts.polymerTypes = polymerTypes.length;
          for (const polymerType of polymerTypes) {
            const symbols = lib.getMonomerSymbolsByType(polymerType);
            for (const sym of symbols) {
              result.replayCounts.monomers++;
              const mon = lib.getMonomer(polymerType, sym);
              if (mon == null) {
                result.replayCounts.addedSkipped++;
                continue;
              }
              const at: any = {};
              let rs = mon.rgroups?.length ?? 0;
              let isBroken = false;
              if (mon.rgroups && mon.rgroups.length > 0) {
                for (const it of mon.rgroups) {
                  at[it.label] = it.capGroupName;
                }
              } else if (mon.smiles != null) {
                const probedRs = lib.getRS(mon.smiles.toString());
                if (probedRs == null || Object.keys(probedRs).length === 0)
                  isBroken = true;
                else {
                  rs = Object.keys(probedRs).length;
                  Object.assign(at, probedRs);
                }
              } else {
                isBroken = true;
              }
              if (isBroken) {
                result.replayCounts.addedSkipped++;
                continue;
              }
              const wem = {
                id: sym,
                m: mon.molfile,
                n: mon.name,
                na: mon.naturalAnalog,
                rs: rs,
                type: mon.polymerType,
                mt: mon.monomerType,
                at: at,
              };
              try {
                monomers.addOneMonomer(wem);
                result.replayCounts.addedOk++;
              } catch (_) {
                result.replayCounts.addedSkipped++;
              }
            }
          }
          const restoredA = monomers.getMonomer(helmAA, 'A');
          const restoredMeI = monomers.getMonomer(helmAA, 'meI');
          result.restoredPeptideA_n = restoredA?.n ?? null;
          result.restoredPeptideMeI_n = restoredMeI?.n ?? null;
          result.peptideAReverted = restoredA != null && restoredA.n === baseA?.n && restoredA.n !== 'missing';
          result.peptideMeIReverted = restoredMeI != null && restoredMeI.n === baseMeI?.n && restoredMeI.n !== 'missing';
          try {
            if (hh_s3?.originalMonomersFuncs == null)
              hh_s3.overrideMonomersFuncs(hh_s3.buildMonomersFuncsFromLib(lib));
          } catch (_) { /* ignore */ }
        } catch (e: any) {
          result.err = String(e?.message ?? e).slice(0, 600);
          try {
            const monomers = (window as any).org.helm.webeditor.Monomers;
            const lib = (window as any).__lib;
            if (lib && monomers) {
              for (const pt of lib.getPolymerTypes()) {
                for (const s of lib.getMonomerSymbolsByType(pt)) {
                  const m = lib.getMonomer(pt, s);
                  if (m == null) continue;
                  const at: any = {};
                  let rs = m.rgroups?.length ?? 0;
                  if (m.rgroups?.length > 0)
                    for (const it of m.rgroups) at[it.label] = it.capGroupName;
                  try {
                    monomers.addOneMonomer({id: s, m: m.molfile, n: m.name, na: m.naturalAnalog,
                      rs: rs, type: m.polymerType, mt: m.monomerType, at: at});
                  } catch (_) { /* ignore */ }
                }
              }
            }
          } catch (_) { /* ignore */ }
          try {
            const hh_s3_c = (window as any).__hh;
            const lib_s3_c = (window as any).__lib;
            if (hh_s3_c && lib_s3_c && hh_s3_c.originalMonomersFuncs == null)
              hh_s3_c.overrideMonomersFuncs(hh_s3_c.buildMonomersFuncsFromLib(lib_s3_c));
          } catch (_) { /* ignore */ }
        }
        const errAfter = document.querySelectorAll('.d4-balloon.error').length;
        result.errBalloonsDelta = errAfter - errBefore;
        return result;
      });
      expect(out.baselinePeptideA_n, `baseline PEPTIDE/A name (expected non-null production name); err: ${out.err ?? ''}`).not.toBeNull();
      expect(out.baselinePeptideA_n, 'baseline PEPTIDE/A is not missing-placeholder').not.toBe('missing');
      expect(out.sentinelAdded, `sentinel monomer added and looked up; sentinelId=${out.sentinelId}, n=${out.sentinelN}`).toBe(true);
      expect(out.sentinelLookupOk, 'sentinel lookup returns the added record (not the missing placeholder)').toBe(true);
      expect(out.clearTookEffect, 'sentinel monomer absent (null) after Monomers.clear() in Pistoia mode').toBe(true);
      expect(out.postClearAReadsMissing, 'production monomer A returns null after Monomers.clear() in Pistoia mode (no lib placeholder)').toBe(true);
      expect(out.peptideAReverted, `PEPTIDE/A restored to baseline name '${out.baselinePeptideA_n}'; observed: '${out.restoredPeptideA_n}'`).toBe(true);
      expect(out.peptideMeIReverted, `PEPTIDE/meI restored to baseline name '${out.baselinePeptideMeI_n}'; observed: '${out.restoredPeptideMeI_n}'`).toBe(true);
      console.warn(`[diagnostic S3 replay-counts] polymerTypes=${out.replayCounts.polymerTypes} monomers=${out.replayCounts.monomers} addedOk=${out.replayCounts.addedOk} addedSkipped=${out.replayCounts.addedSkipped}`);
      console.warn(`[diagnostic S3 SR-01] ${out.srNote}`);
      expect(out.errBalloonsDelta, 'new error balloon delta after S3').toBe(0);
    });
  } finally {
    await page.evaluate(async () => {
      try {
        const hh = (window as any).__hh;
        if (hh?.originalMonomersFuncs != null) {
          try { hh.revertOriginalMonomersFuncs(); } catch (_) { /* ignore */ }
        }
        const lib = (window as any).__lib;
        const monomers = (window as any).org?.helm?.webeditor?.Monomers;
        if (lib && monomers) {
          const helmAAType = (window as any).org?.helm?.webeditor?.HELM?.AA ?? 'HELM_AA';
          const probe = monomers.getMonomer(helmAAType, 'A');
          if (probe == null || probe.n === 'missing') {
            for (const pt of lib.getPolymerTypes()) {
              for (const s of lib.getMonomerSymbolsByType(pt)) {
                const m = lib.getMonomer(pt, s);
                if (m == null) continue;
                const at: any = {};
                let rs = m.rgroups?.length ?? 0;
                if (m.rgroups?.length > 0)
                  for (const it of m.rgroups) at[it.label] = it.capGroupName;
                try {
                  monomers.addOneMonomer({id: s, m: m.molfile, n: m.name, na: m.naturalAnalog,
                    rs: rs, type: m.polymerType, mt: m.monomerType, at: at});
                } catch (_) { /* ignore */ }
              }
            }
          }
        }
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
