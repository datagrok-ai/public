/**
 * Worker jobs by the "Chemical engine" setting (README.md beside this file): in Crux mode the Crux pool runs them, with
 * RDKit's own worker code on the Crux module, so the results have the RDKit pool's shape and options; the RDKit pool
 * does the molecules Crux cannot read, and the whole call when Crux fails. In RDKit mode, the RDKit pool. In Crux mode
 * the RDKit pool starts only when a molecule or a call needs it.
 */
import type * as DG from 'datagrok-api/dg';
import {Fingerprint} from '../utils/chem-common';
import type {IFpResult} from '../rdkit-service/rdkit-service-worker-similarity';
import {getRdKitService} from '../utils/chem-common-rdkit';
import {getCruxService} from '../crux/crux-searches';
import {ChemicalEngine, isCruxRequested} from './chem-engine';
import {cruxOrRdKitAsync, warnCruxFallback} from './fallback';
import type {RuleId} from '../panels/structural-alerts';
import {FingerprintOptions, fingerprintOptions, needsCrux} from '../utils/fingerprint-options';
import {cruxReadsCxExtension} from '../crux/crux-service';
import {hasNewLines} from '../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';

/** Tests: molecules the Crux workers get as unreadable (see {@link setCruxUnreadableForTests}). */
let unreadableForTests: Set<string> | null = null;

/** Tests: the Crux workers get these molecules as one they cannot read (null: none); RDKit gets them as they are. */
export function setCruxUnreadableForTests(molecules: string[] | null): void {
  unreadableForTests = molecules ? new Set(molecules) : null;
}

/** The molecules as the Crux workers get them (also analysis-jobs.ts). */
export function forCrux(molecules: string[]): string[] {
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

/**
 * `RdKitService.convertMolNotation(molecules, targetNotation, kekulize)` by the engine in effect: in Crux mode RDKit.js
 * 2026.03.6's strings (the RDKit Crux reproduces), the RDKit pool's for the rows Crux leaves unconverted and for a
 * CXSMILES whose query features Crux does not read (a SMILES over Chem's length limit is no one's).
 */
export async function convertMolNotationByEngine(molecules: string[], targetNotation: DG.chem.Notation,
  kekulize: boolean): Promise<string[]> {
  const rdkit = async (rows: number[] | null) => (await getRdKitService()).convertMolNotation(
    rows ? rows.map((i) => molecules[i]) : molecules, targetNotation, kekulize);
  if (!isCruxRequested())
    return rdkit(null);
  return cruxOrRdKitAsync('convertMolNotation', async () => {
    const cx = molecules.map((m) => !!m && !cruxReadsCxExtension(m));
    const {results, none, unconverted} = await (await getCruxService()).convertMolNotation(
      forCrux(molecules).map((m, i) => cx[i] ? '' : m), targetNotation, kekulize);
    const tooLong = (m: string) => !hasNewLines(m) && m.length > MAX_SMILES_LENGTH;
    const rest = [...unconverted.filter((i) => !tooLong(molecules[i])), ...cx.flatMap((c, i) => c ? [i] : [])];
    if (rest.length) {
      const patch = await rdkit(rest);
      rest.forEach((row, k) => results[row] = patch[k]);
      const done = rest.filter((_, k) => patch[k] !== none);
      if (done.length) {
        warnCruxFallback('convertMolNotation', `RDKit converted ${done.length} molecules Crux did not, e.g. ` +
          molecules[done[0]]);
      }
    }
    return results;
  }, () => rdkit(null));
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
