/**
 * Chem's Synthon Space Search in Crux mode, the worker's side (crux-service.ts sends the jobs): chem-wasm's
 * ChemSynthonSpace, crux-core's exact port of RDKit 2025.09.3's SynthonSpaceSearch as scripts/synthon_search.py runs
 * it. A worker keeps the libraries it read (by the caller's key) for the next searches.
 */
import type {ChemSynthonSpace} from '../../../vendor/chem-wasm/chem_wasm.js';
import type {WorkerJob} from '../crux-api';

type SynthonModule = {ChemSynthonSpace: new (text: string) => ChemSynthonSpace};

/** The libraries this worker holds, least recently used first. */
const spaces = new Map<string, ChemSynthonSpace>();
/** A library holds its synthons' molecules and fingerprints (a few MB for Chem's 5,567 synthons). */
const MAX_SPACES = 2;

/** The library of `key`, read from `text` when this worker does not hold it. */
function spaceOf(module: SynthonModule, key: string, text: string): ChemSynthonSpace {
  let space = spaces.get(key);
  if (space)
    spaces.delete(key);
  else
    space = new module.ChemSynthonSpace(text);
  spaces.set(key, space);
  for (const [k, s] of spaces) {
    if (spaces.size <= MAX_SPACES)
      break;
    s.free();
    spaces.delete(k);
  }
  return space;
}

/** The jobs by name (the Crux worker's `op`). */
export const SYNTHON_JOBS: {[op: string]: WorkerJob} = {
  /** Whether the library's similarity fingerprints are built on this worker (read when it holds none). */
  synthonHasFingerprints: ({crux}, key: string, text: string): boolean =>
    spaceOf(crux, key, text).hasFingerprints(),
  /** Part `part` of `parts` of the library's fingerprint build, as bytes (transferred); a worker that does not hold
   * the library reads it for this and lets it go. */
  synthonFingerprintPart: ({crux, transfer}, key: string, text: string, part: number, parts: number): Uint8Array => {
    const held = spaces.get(key);
    const space = held ?? new crux.ChemSynthonSpace(text);
    try {
      const bytes = space.fingerprintPart(part, parts);
      transfer.push(bytes.buffer);
      return bytes;
    } finally {
      if (!held)
        space.free();
    }
  },
  /** The fingerprints from every part (each after its length, u32 little-endian). */
  synthonInstallFingerprints: ({crux}, key: string, text: string, parts: Uint8Array): void =>
    spaceOf(crux, key, text).installFingerprints(parts),
  /** scripts/synthon_search.py on the library: its dataframe's columns as JSON (a "not supported by Crux" refusal for
   * a query Crux does not take: the caller runs the script). */
  synthonSearch: ({crux}, key: string, text: string, molecule: string, searchType: string, maxHits: number,
    similarityCutoff: number, includeSynthons: boolean): string =>
    spaceOf(crux, key, text).search(molecule, searchType, maxHits, similarityCutoff, includeSynthons),
};
