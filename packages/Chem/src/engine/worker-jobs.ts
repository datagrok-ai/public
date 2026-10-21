/**
 * Worker jobs by the "Chemical engine" setting (README.md beside this file): in Crux mode the Crux pool runs them, with
 * RDKit's own worker code on the Crux module, so the results have the RDKit pool's shape and options; the RDKit pool
 * does the molecules Crux cannot read, and the whole call when Crux fails. In RDKit mode, the RDKit pool. In Crux mode
 * the RDKit pool starts only when a molecule or a call needs it.
 */
import {Fingerprint} from '../utils/chem-common';
import type {IFpResult} from '../rdkit-service/rdkit-service-worker-similarity';
import {getRdKitService} from '../utils/chem-common-rdkit';
import {getCruxService} from '../crux/crux-searches';
import {ChemicalEngine, isCruxRequested} from './chem-engine';
import {cruxOrRdKitAsync, warnCruxFallback} from './fallback';
import type {RuleId} from '../panels/structural-alerts';
import {FingerprintOptions, fingerprintOptions, needsCrux} from '../utils/fingerprint-options';

/** Tests: molecules the Crux workers get as unreadable (see {@link setCruxUnreadableForTests}). */
let unreadableForTests: Set<string> | null = null;

/** Tests: the Crux workers get these molecules as one they cannot read (null: none); RDKit gets them as they are. */
export function setCruxUnreadableForTests(molecules: string[] | null): void {
  unreadableForTests = molecules ? new Set(molecules) : null;
}

/** The molecules as the Crux workers get them. */
function forCrux(molecules: string[]): string[] {
  const unreadable = unreadableForTests;
  return unreadable ? molecules.map((m) => unreadable.has(m) ? 'C1CC(' : m) : molecules;
}

/** The suffix of a fingerprint cache's key: none for RDKit's fingerprints, '.Crux' for Crux's. */
export function fingerprintsKey(engine: ChemicalEngine): string {
  return engine === ChemicalEngine.Crux ? '.Crux' : '';
}

/**
 * `RdKitService.getFingerprints(type, molecules, smiles, options)` by the engine in effect (by RDKit when `engine` says
 * so: the RDKit substructure search's own fingerprints), and the engine that answered. FCFP and chirality
 * (`options.useFeatures`, `options.useChirality`) are Crux's in both modes, with no RDKit behind them: Chem's RDKit.js
 * has neither (utils/fingerprint-options.ts `needsCrux`).
 */
export async function fingerprintsByEngine(type: Fingerprint, molecules: string[], smiles: boolean,
  engine: ChemicalEngine = isCruxRequested() ? ChemicalEngine.Crux : ChemicalEngine.RDKit,
  options?: FingerprintOptions | null,
): Promise<IFpResult & {engine: ChemicalEngine}> {
  const opts = fingerprintOptions(type, options);
  const cruxOnly = needsCrux(opts);
  const rdkit = async () => {
    if (cruxOnly) { // Crux refused the options or failed (logged by the fallback); RDKit.js has neither
      throw new Error(`Chem | Crux could not compute these ${type} fingerprints (${JSON.stringify(opts)}), and ` +
        'only Crux computes FCFP and chirality: see the warning before');
    }
    return {...await (await getRdKitService()).getFingerprints(type, molecules, smiles, opts),
      engine: ChemicalEngine.RDKit};
  };
  if (engine === ChemicalEngine.RDKit && !cruxOnly)
    return rdkit();
  return cruxOrRdKitAsync('fingerprints', async () => {
    const result = await (await getCruxService()).getFingerprints(type, forCrux(molecules), smiles, opts);
    const {fps, smiles: canonical, unread} = result;
    if (unread.length && !cruxOnly) { // the molecules Crux cannot read: the RDKit pool's fingerprints for them
      const patch = await (await getRdKitService()).getFingerprints(type, unread.map((i) => molecules[i]), smiles,
        opts);
      unread.forEach((row, k) => {
        fps[row] = patch.fps[k];
        if (canonical && patch.smiles)
          canonical[row] = patch.smiles[k];
      });
      const read = unread.filter((_, k) => patch.fps[k] !== null); // (one RDKit cannot read either is no fallback)
      if (read.length) {
        warnCruxFallback('fingerprints', `RDKit computed ${read.length} molecules Crux cannot read, e.g. ` +
          molecules[read[0]]);
      }
    }
    return {fps, smiles: canonical, engine: ChemicalEngine.Crux};
  }, rdkit);
}

/** `RdKitService.getStructuralAlerts(alerts, molecules)` by the engine in effect. */
export async function structuralAlertsByEngine(alerts: {[rule in RuleId]?: string[]}, molecules: string[]):
  Promise<[RuleId, boolean[]][]> {
  const rdkit = async () => (await getRdKitService()).getStructuralAlerts(alerts, molecules);
  if (!isCruxRequested())
    return rdkit();
  return cruxOrRdKitAsync('structural alerts', async () => {
    const result = await (await getCruxService()).getStructuralAlerts(alerts, forCrux(molecules));
    const {alerts: found, unread} = result;
    if (unread.length) { // the molecules Crux cannot read: the RDKit pool's matches for them
      const patch = new Map(await (await getRdKitService()).getStructuralAlerts(alerts,
        unread.map((i) => molecules[i])));
      for (const [rule, flags] of found)
        unread.forEach((row, k) => flags[row] = patch.get(rule)?.[k] ?? flags[row]);
      const matched = unread.filter((_, k) => [...patch.values()].some((flags) => flags[k]));
      if (matched.length) {
        warnCruxFallback('structural alerts', `RDKit matched ${matched.length} molecules Crux cannot read, e.g. ` +
          molecules[matched[0]]);
      }
    }
    return found;
  }, rdkit);
}

/** What structural alerts and pharmacophore features detect with: the engine in effect, RDKit behind it. */
export const alertsByEngine = {getStructuralAlerts: structuralAlertsByEngine};
