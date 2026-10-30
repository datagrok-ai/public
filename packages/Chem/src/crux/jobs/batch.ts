/**
 * The Crux worker's whole-column chemistry (crux.worker.ts runs it; engine/worker-jobs.ts sends it): structural alerts
 * and pharmacophore features by crux-core's `ChemSmartsSet`, fingerprints by `chemFingerprints` (chem-wasm's
 * chem_batch.rs), a chunk of molecules per call instead of one RDKit.js-shaped call per molecule (and per SMARTS). The
 * results are those the RDKit worker's own code gives on the Crux module (getStructuralAlerts, getFingerprints in
 * rdkit-service/): each molecule read as getMolSafe reads it, a group's flag set by its first matching SMARTS, the same
 * fingerprint calls; identical on ChEMBL 100k. The SMARTS sets skip a SMARTS on a molecule that lacks the atoms or
 * bonds it needs, which skips no match (crux-core's chem_batch.rs says why).
 */
import type * as ChemWasm from '../../../vendor/chem-wasm/chem_wasm.js';
import type {IFpResult} from '../../rdkit-service/rdkit-service-worker-similarity';
import {Fingerprint, hasNewLines} from '../../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../../utils/chem-constants';
import {FingerprintOptions, assertFingerprintOptions, fingerprintDetails, fingerprintOptions}
  from '../../utils/fingerprint-options';
import type {WorkerJob} from '../crux-api';

type Crux = typeof ChemWasm;

/** The molecules as getMolSafe reads them: a SMILES over Chem's length limit is none (no match, no fingerprint). */
function readable(molecules: string[]): string[] {
  return molecules.map((m) => m && !hasNewLines(m) && m.length > MAX_SMILES_LENGTH ? '' : m);
}

/** A compiled set and, by rule, the SMARTS of it Crux reads nothing from (their indices in the rule's list). */
interface CompiledRules {
  set: ChemWasm.ChemSmartsSet;
  rules: string[];
  unreadSmarts: {[rule: string]: number[]};
}

/** The last few rule sets compiled (a column's chunks, the panels' repeated calls), the most recently used last. */
const compiled = new Map<string, CompiledRules>();
const MAX_COMPILED = 4;

function compile(crux: Crux, alerts: {[rule: string]: string[]}): CompiledRules {
  const key = JSON.stringify(alerts);
  let c = compiled.get(key);
  if (c)
    compiled.delete(key);
  else {
    const rules = Object.keys(alerts);
    const smarts = rules.flatMap((r) => alerts[r]);
    const set = new crux.ChemSmartsSet(smarts, new Uint32Array(rules.map((r) => alerts[r].length)));
    const unreadSmarts: {[rule: string]: number[]} = {};
    const at = rules.flatMap((r) => alerts[r].map((_, k) => [r, k] as const));
    for (const k of set.unreadSmarts())
      (unreadSmarts[at[k][0]] ??= []).push(at[k][1]);
    c = {set, rules, unreadSmarts};
    while (compiled.size >= MAX_COMPILED) {
      const [oldest, entry] = compiled.entries().next().value!;
      entry.set.free();
      compiled.delete(oldest);
    }
  }
  compiled.set(key, c);
  return c;
}

/**
 * Which molecules match each rule's SMARTS, as the RDKit worker's `getStructuralAlerts` matches them; the rows with a
 * molecule Crux cannot read, and by rule the SMARTS Crux reads nothing from, as indices into its list (RDKit's to
 * match: engine/worker-jobs.ts).
 */
export function structuralAlerts(crux: Crux, alerts: {[rule: string]: string[]}, molecules: string[]):
  {alerts: {[rule: string]: boolean[]}, unread: number[], unreadSmarts: {[rule: string]: number[]}} {
  const {set, rules, unreadSmarts} = compile(crux, alerts);
  const {flags, unread} = set.match(readable(molecules));
  const found: {[rule: string]: boolean[]} = {};
  rules.forEach((rule, g) => {
    const f = new Array<boolean>(molecules.length);
    for (let i = 0; i < molecules.length; i++)
      f[i] = flags[i * rules.length + g] === 1;
    found[rule] = f;
  });
  return {alerts: found, unread: Array.from(unread), unreadSmarts};
}

/** The module's fingerprint kinds (`get_<kind>_fp_as_uint8array`). */
const KIND: {[type in Fingerprint]: string} = {
  [Fingerprint.Morgan]: 'morgan',
  [Fingerprint.Pattern]: 'pattern',
  [Fingerprint.RDKit]: 'rdkit',
  [Fingerprint.MACCS]: 'maccs',
  [Fingerprint.AtomPair]: 'atom_pair',
  [Fingerprint.TopologicalTorsion]: 'topological_torsion',
};

/**
 * Fingerprints as the RDKit worker's `getFingerprints` computes them (Morgan after `remove_hs_in_place`, none for a
 * query molecule), the canonical SMILES when asked, and the rows with a molecule Crux cannot read. Each fingerprint
 * its own buffer, as before.
 */
export function fingerprints(crux: Crux, type: Fingerprint, molecules: string[], smiles: boolean,
  options?: FingerprintOptions | null): IFpResult & {unread: number[]} {
  const details = fingerprintDetails(type, fingerprintOptions(type, options));
  const r = crux.chemFingerprints(KIND[type], details, readable(molecules), smiles);
  const fps = new Array<Uint8Array | null>(molecules.length);
  for (let i = 0; i < molecules.length; i++)
    fps[i] = r.present[i] ? r.bytes.slice(r.offsets[i], r.offsets[i + 1]) : null;
  return {fps, smiles: r.smiles, unread: Array.from(r.unread)};
}

/** The jobs by name (the Crux worker's `op`): the RDKit pool's calls of the same names, a chunk in one Crux call. */
export const BATCH_JOBS: {[op: string]: WorkerJob} = {
  getFingerprints: ({crux, module}, type: Fingerprint, molecules: string[], smiles: boolean,
    options?: FingerprintOptions | null) => {
    assertFingerprintOptions(module, type, options); // a refusal: the RDKit pool's call
    return fingerprints(crux, type, molecules, smiles, options);
  },
  getStructuralAlerts: ({crux}, alerts: {[rule: string]: string[]}, molecules: string[]) =>
    structuralAlerts(crux, alerts, molecules),
};
