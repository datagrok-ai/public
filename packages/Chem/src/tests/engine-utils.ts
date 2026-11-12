/**
 * How Chem's tests run the two chemical engines (src/engine/README.md). Every engine switch of a test goes through
 * these helpers, so that removing RDKit is a deletion in the tests too:
 * - `inEngine(ChemicalEngine.Crux, ...)`: a Crux leg, checked against pinned values (RDKit 2026.03's, a script's, the
 *   container's). It stays.
 * - `inEngine(ChemicalEngine.RDKit, ...)`: an RDKit leg (RDKit mode's own path, or the reference a Crux leg is
 *   compared with), and `inBothEngines(...)`: the same run in RDKit mode, then in Crux mode, compared. The RDKit leg
 *   goes; a Crux leg compared with it gets pinned values or goes too.
 * - `expectNoFallback(...)`: RDKit did none of a Crux leg's calls. It goes (there is no fallback then), its run stays.
 * - `withCruxUnreadable(...)` and the other `...ForTests` failures of src/engine and src/crux (a failed or stalled
 *   load, a trap, a failed worker job): fallback tests, RDKit doing what Crux could not. They go.
 * `setChemicalEngine` reaches this bundle's code only: a function called through the platform (`grok.functions.call`,
 * `DG.Func.apply`) runs the published package's code, in the stand's engine; call `PackageFunctions` here instead.
 */
import {expect} from '@datagrok-libraries/test/src/test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine, cruxReady, setChemicalEngine} from '../engine/chem-engine';
import {cruxFallbackCount} from '../engine/fallback';
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

/** `run()` with the Crux workers reading none of `molecules` (RDKit does them, each counted as a fallback). */
export async function withCruxUnreadable<T>(molecules: string[], run: () => Promise<T>): Promise<T> {
  setCruxUnreadableForTests(molecules);
  try {
    return await run();
  } finally {
    setCruxUnreadableForTests(null);
  }
}

/** How many times the RDKit worker pool was asked for (each `getRdKitService()` counts). */
export const rdkitPoolAsks = (): number => chemCommonRdKit._rdKitService?.timesInitialized ?? 0;
