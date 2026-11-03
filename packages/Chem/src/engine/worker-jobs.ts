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
import {ChemicalEngine, byEngineAsync, isCruxRequested} from './chem-engine';
import {cruxOrRdKitAsync, rdkitForRows} from './fallback';
import type {RuleId} from '../panels/structural-alerts';
import {FingerprintOptions, fingerprintOptions, needsCrux} from '../utils/fingerprint-options';
import {cruxReadsCxExtension} from '../crux/crux-service';
import {hasNewLines} from '../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';
import {inchisFromInputs} from '../crux/crux-inchi';

/** Tests: molecules the Crux workers get as unreadable (see {@link setCruxUnreadableForTests}). */
let unreadableForTests: Set<string> | null = null;

/** Tests: the Crux workers get these molecules (and alert SMARTS) as ones they cannot read (null: none); RDKit gets
 * them as they are. */
export function setCruxUnreadableForTests(molecules: string[] | null): void {
  unreadableForTests = molecules ? new Set(molecules) : null;
}

/** The molecules (or SMARTS) as the Crux workers get them: tests mark some unreadable. */
export function forCrux(molecules: string[]): string[] {
  const unreadable = unreadableForTests;
  return unreadable ? molecules.map((m) => unreadable.has(m) ? 'C1CC(' : m) : molecules;
}

/** A column job's molecules as the Crux workers get them ({@link forCrux}), '' for a CXSMILES whose query features
 * Crux's reader skips: those rows (`skipped`) are the RDKit side's. */
export function cruxInput(molecules: string[]): {input: string[], skipped: number[]} {
  const skipped = molecules.flatMap((m, i) => m && !cruxReadsCxExtension(m) ? [i] : []);
  const input = forCrux(molecules);
  if (!skipped.length)
    return {input, skipped};
  const rows = new Set(skipped);
  return {input: input.map((m, i) => rows.has(i) ? '' : m), skipped};
}

/** The suffix of a fingerprint cache's key: none for RDKit's fingerprints, '.Crux' for Crux's. */
export function fingerprintsKey(engine: ChemicalEngine): string {
  return engine === ChemicalEngine.Crux ? '.Crux' : '';
}

/**
 * `RdKitService.getFingerprints(type, molecules, smiles, options)` by the engine in effect (by RDKit when `engine` says
 * so: the RDKit substructure search's own fingerprints), and the engine that answered. FCFP and chirality
 * (`options.useFeatures`, `options.useChirality`) are Crux's in both modes, with no RDKit behind them: Chem's RDKit.js
 * has neither (utils/fingerprint-options.ts `needsCrux`). In Crux mode the RDKit pool computes the molecules Crux
 * cannot read and a CXSMILES whose query features Crux's reader skips (but FCFP and chirality).
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
    // FCFP and chirality: every molecule Crux's, as no RDKit is behind them
    const {input, skipped} = cruxOnly ? {input: forCrux(molecules), skipped: []} : cruxInput(molecules);
    const result = await (await getCruxService()).getFingerprints(type, input, smiles, opts);
    const {fps, smiles: canonical} = result;
    if (!cruxOnly) {
      await rdkitForRows('fingerprints', [...result.unread, ...skipped],
        async (rows) => (await getRdKitService()).getFingerprints(type, rows.map((i) => molecules[i]), smiles, opts),
        (patch, k, row) => {
          fps[row] = patch.fps[k];
          if (canonical && patch.smiles)
            canonical[row] = patch.smiles[k];
          return patch.fps[k] !== null;
        }, (row) => molecules[row]);
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
  return byEngineAsync('convertMolNotation', async () => {
    const {input, skipped} = cruxInput(molecules);
    const {results, none, unconverted} = await (await getCruxService()).convertMolNotation(input, targetNotation,
      kekulize);
    const tooLong = (m: string) => !hasNewLines(m) && m.length > MAX_SMILES_LENGTH;
    await rdkitForRows('convertMolNotation', [...unconverted.filter((i) => !tooLong(molecules[i])), ...skipped], rdkit,
      (patch, k, row) => (results[row] = patch[k]) !== none, (row) => molecules[row]);
    return results;
  }, () => rdkit(null));
}

/**
 * `RdKitService.getInchis(molecules)` (`getInchiKeys` when `keys`) by the engine in effect: in Crux mode InChI's input
 * by the Crux workers and the library by crux-sketch's InChI module (crux/crux-inchi.ts), the RDKit pool for a molecule
 * Crux cannot read and a CXSMILES whose query features Crux's reader skips.
 */
export async function inchisByEngine(molecules: string[], keys: boolean): Promise<string[]> {
  const rdkit = async (rows: string[]) => {
    const service = await getRdKitService();
    return keys ? service.getInchiKeys(rows) : service.getInchis(rows);
  };
  return byEngineAsync('InChI', async () => {
    const {input, skipped} = cruxInput(molecules);
    const inputs = await (await getCruxService()).inchiInputs(input);
    for (const row of skipped)
      inputs[row] = null; // (Crux's input of '' is '')
    const results = await inchisFromInputs(inputs, keys);
    await rdkitForRows('InChI', results.flatMap((r, i) => r === null ? [i] : []),
      (rows) => rdkit(rows.map((i) => molecules[i])), (patch, k, row) => !!(results[row] = patch[k]),
      (row) => molecules[row]);
    return results as string[];
  }, () => rdkit(molecules));
}

/**
 * `RdKitService.getStructuralAlerts(alerts, molecules)` by the engine in effect. In Crux mode the RDKit pool matches
 * the molecules Crux cannot read, a CXSMILES whose query features Crux's reader skips, and a SMARTS Crux reads nothing
 * from (on the molecules no other SMARTS of its rule matched).
 */
export async function structuralAlertsByEngine(alerts: {[rule in RuleId]?: string[]}, molecules: string[]):
  Promise<[RuleId, boolean[]][]> {
  const rdkit = async () => (await getRdKitService()).getStructuralAlerts(alerts, molecules);
  return byEngineAsync('structural alerts', async () => {
    const {input, skipped} = cruxInput(molecules);
    const forCruxAlerts = Object.fromEntries(Object.entries(alerts).map(([rule, smarts]) => [rule, forCrux(smarts!)]));
    const result = await (await getCruxService()).getStructuralAlerts(forCruxAlerts, input);
    const {alerts: found} = result;
    const unread = [...result.unread, ...skipped];
    const rdkitAlerts = async (some: {[rule in RuleId]?: string[]}, rows: number[]) =>
      new Map(await (await getRdKitService()).getStructuralAlerts(some, rows.map((i) => molecules[i])));
    const unreadRows = new Set(unread);
    for (const [rule, flags] of found) { // a SMARTS Crux reads nothing from: RDKit's matches of it
      const smarts = result.unreadSmarts[rule]?.map((k) => alerts[rule]![k]);
      if (!smarts)
        continue;
      await rdkitForRows('structural alerts', flags.flatMap((f, i) => f || unreadRows.has(i) ? [] : [i]),
        (rows) => rdkitAlerts({[rule]: smarts}, rows),
        (patch, k, row) => !!patch.get(rule)?.[k] && (flags[row] = true), () => `the SMARTS ${smarts[0]}`);
    }
    // the molecules Crux cannot read: the RDKit pool's matches for them
    await rdkitForRows('structural alerts', unread, (rows) => rdkitAlerts(alerts, rows), (patch, k, row) => {
      for (const [rule, flags] of found)
        flags[row] = patch.get(rule)?.[k] ?? flags[row];
      return [...patch.values()].some((flags) => flags[k]);
    }, (row) => molecules[row]);
    return found;
  }, rdkit);
}

/** What structural alerts and pharmacophore features detect with: the engine in effect, RDKit behind it. */
export const alertsByEngine = {getStructuralAlerts: structuralAlertsByEngine};
