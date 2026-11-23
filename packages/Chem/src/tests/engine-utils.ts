/**
 * How Chem's tests run the two chemical engines (src/engine/README.md). Every engine switch of a test goes through
 * these helpers, so that removing RDKit is a deletion in the tests too:
 * - `inEngine(ChemicalEngine.Crux, ...)`: a Crux leg, checked against pinned values (RDKit 2026.03's, a script's, the
 *   container's). It stays.
 * - `inEngine(ChemicalEngine.RDKit, ...)`: an RDKit leg (RDKit mode's own path, or the reference a Crux leg is
 *   compared with), and `inBothEngines(...)`: the same run in RDKit mode, then in Crux mode, compared. The RDKit leg
 *   goes; a Crux leg compared with it gets pinned values or goes too.
 * - `expectNoFallback(...)`: RDKit did none of a Crux leg's calls. It goes (there is no fallback then), its run stays.
 * - `expectFallback(...)`: RDKit did a call once, on purpose (a `...ForTests` failure); kept off the platform's log.
 *   It goes, with its run.
 * - `withCruxUnreadable(...)` and the other `...ForTests` failures of src/engine and src/crux (a failed or stalled
 *   load, a trap, a failed worker job): fallback tests, RDKit doing what Crux could not. They go.
 * `setChemicalEngine` reaches this bundle's code only: a function called through the platform (`grok.functions.call`,
 * `DG.Func.apply`) runs the published package's code, in the stand's engine; call `PackageFunctions` here instead.
 */
import {expect} from '@datagrok-libraries/test/src/test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine, cruxReady, setChemicalEngine} from '../engine/chem-engine';
import {cruxFallbackCount, cruxFallbackReason, quietFallbacksForTests} from '../engine/fallback';
import {setCruxUnreadableForTests} from '../engine/worker-jobs';
import {getCruxService} from '../crux/crux-searches';

/** `run()` in `engine`'s mode (Crux: its module loaded), the stand's engine again after. */
export async function inEngine<T>(engine: ChemicalEngine, run: () => Promise<T>): Promise<T> {
  setChemicalEngine(engine);
  try {
    if (engine === ChemicalEngine.Crux)
      expect(await cruxReady(), true, 'the Crux module loads');
    return await run();
  } finally {
    setChemicalEngine(null);
  }
}

/** What the Crux leg of {@link inBothEngines} must do alone: none of the `noFallback` calls by RDKit, and
 * (`noRdkitPool`) no RDKit worker asked for. */
export interface CruxAlone {
  noFallback?: string[];
  noRdkitPool?: boolean;
}

/** `run(engine)` in RDKit mode, then in Crux mode: both results. */
export async function inBothEngines<T>(run: (engine: ChemicalEngine) => Promise<T>, alone: CruxAlone = {}):
  Promise<{rdkit: T, crux: T}> {
  const rdkit = await inEngine(ChemicalEngine.RDKit, () => run(ChemicalEngine.RDKit));
  const crux = await inEngine(ChemicalEngine.Crux, async () => {
    if (alone.noRdkitPool)
      await getCruxService(); // the Crux pool up first
    const asks = rdkitPoolAsks();
    const result = await expectNoFallback(alone.noFallback ?? [], () => run(ChemicalEngine.Crux));
    if (alone.noRdkitPool)
      expect(rdkitPoolAsks(), asks, 'the Crux leg asked for the RDKit workers');
    return result;
  });
  return {rdkit, crux};
}

/** `run()`, RDKit doing none of `calls` for Crux (their fallback counts unchanged). */
export async function expectNoFallback<T>(calls: string[], run: () => Promise<T>): Promise<T> {
  const counts = () => calls.map((c) => `${c} ${cruxFallbackCount(c)}`).join(', ');
  const before = counts();
  const result = await run();
  expect(counts(), before, 'RDKit did it for Crux');
  return result;
}

/** `run()`, RDKit doing `call` for Crux exactly once (a fallback the test makes on purpose, counted but kept off the
 * platform's log: fallback.ts `quietFallbacksForTests`); its result and that fallback's reason. */
export async function expectFallback<T>(call: string, run: () => Promise<T>): Promise<{result: T, reason: string}> {
  const before = cruxFallbackCount(call);
  const result = await quietFallbacksForTests(run);
  expect(cruxFallbackCount(call), before + 1, `the fallbacks of ${call}`);
  return {result, reason: cruxFallbackReason(call) ?? ''};
}

/** `run()` with the Crux workers reading none of `molecules` (RDKit does them, each counted as a fallback, kept off the
 * platform's log: they are the test's, not real use's). */
export async function withCruxUnreadable<T>(molecules: string[], run: () => Promise<T>): Promise<T> {
  setCruxUnreadableForTests(molecules);
  try {
    return await quietFallbacksForTests(run);
  } finally {
    setCruxUnreadableForTests(null);
  }
}

/** How many times the RDKit worker pool was asked for (each `getRdKitService()` counts). */
export const rdkitPoolAsks = (): number => chemCommonRdKit._rdKitService?.timesInitialized ?? 0;

/** The owner's query of 2026-11-20, as the substructure filter keeps it (FILTER_SCAFFOLD_TAG: the Crux sketcher's
 * molfile, aliases translated): propane drawn aromatic / aliphatic / aromatic, a plain C grown on its third atom. The
 * filter in Crux mode reads it by run20's drawn-query meanings (a plain bond at an aromatic-marked atom also matches an
 * aromatic bond: 149 mol1K rows); read literally it matches none of them. */
export const OWNER_PROPANE_QUERY = [
  '', '     RDKit          2D', '', '  0  0  0  0  0  0  0  0  0  0999 V3000', 'M  V30 BEGIN CTAB',
  'M  V30 COUNTS 4 3 3 0 0', 'M  V30 BEGIN ATOM', 'M  V30 1 C -1.299038 -0.250000 0.000000 0',
  'M  V30 2 C 0.000000 0.500000 0.000000 0', 'M  V30 3 C 1.299038 -0.250000 0.000000 0',
  'M  V30 4 C 2.598076 0.500000 0.000000 0', 'M  V30 END ATOM',
  'M  V30 BEGIN BOND', 'M  V30 1 1 1 2', 'M  V30 2 1 2 3', 'M  V30 3 1 3 4', 'M  V30 END BOND', 'M  V30 BEGIN SGROUP',
  'M  V30 1 DAT 0 ATOMS=(1 1) QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA="[#6&a]"',
  'M  V30 2 DAT 0 ATOMS=(1 2) QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA="[#6&A]"',
  'M  V30 3 DAT 0 ATOMS=(1 3) QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA="[#6&a]"',
  'M  V30 END SGROUP', 'M  V30 END CTAB', 'M  END', ''].join('\n');
/** mol1K's row 444, one of the 149 rows {@link OWNER_PROPANE_QUERY} finds in the filter. */
export const MOL1K_ROW_444 = 'O=C(O)C(=O)CC(=O)c1ccsc1Cc1ccccc1';
