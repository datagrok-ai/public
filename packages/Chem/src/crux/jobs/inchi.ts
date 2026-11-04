/**
 * A column's InChIs in the Crux workers: InChI's input from the worker's Crux module, as before, and the InChI library
 * of crux-sketch's InChI module (vendor/crux-sketch, the module the main thread's sketchers and panels use) on it in
 * the same worker, each worker loading the module on its first InChI job. The strings are crux/crux-inchi.ts'
 * (RDKit.js's `get_inchi`, `get_inchikey_for_inchi`; the container's for Map Identifiers), the library run across the
 * workers instead of molecule by molecule on the main thread. A failure of the library fails the job, so the caller's
 * RDKit path (the RDKit pool, the container) does the column, as a failure on the main thread did. Map Identifiers'
 * InChI to SMILES and to InChIKey run here too.
 */
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {ChemWasm} from '../../engine/crux-module';
import {WorkerJob, throwIfTrap} from '../crux-api';
import {inchiInputs} from './calculations';

/** crux-sketch's InChI module, as its `start()` gives it (the calls Chem makes of it). */
interface InchiLibrary {
  inchiFromInput(input: Uint8Array): {inchi: string};
  inchiToInchiKey(inchi: string): string;
  inchiToMol(inchi: string): {smiles: string};
}

let library: Promise<InchiLibrary> | null = null;

/** crux-sketch's InChI module, loaded once per worker (a failed load is tried again by the next job). */
function inchiLibrary(): Promise<InchiLibrary> {
  // the module file of vendor/crux-sketch that crux-sketch.js imports for its InChI calls (its wasm beside it in dist/)
  // @ts-ignore: a vendored ES module without declarations
  library ??= import('../../../vendor/crux-sketch/inchi-CzWcatSB.js').then((m) => m.start() as Promise<InchiLibrary>)
    .catch((e: unknown) => {
      library = null;
      throw e;
    });
  return library;
}

/** An InChI (its key when `key`) of `input`, as crux/crux-inchi.ts' `inchiOf` writes it: '' for ''. */
function inchiOf(lib: InchiLibrary, input: Uint8Array | '', key: boolean): string {
  const inchi = input ? lib.inchiFromInput(input).inchi : '';
  return key ? (inchi ? lib.inchiToInchiKey(inchi) : '') : inchi;
}

/** The jobs by name (the Crux worker's `op`). */
export const INCHI_JOBS: {[op: string]: WorkerJob} = {
  /** `inchisByEngine` (engine/worker-jobs.ts): an InChI or key per molecule, null where Crux cannot read it. */
  'inchis': async ({module}: {module: RDModule}, molecules: string[], keys: boolean) => {
    const lib = await inchiLibrary();
    return inchiInputs(module, molecules, []).map((input) => input === null ? null : inchiOf(lib, input, keys));
  },
  /** Map Identifiers' InChIs (python-free.ts): the container's, '' where it writes none (`chemInchiInputs`' null). */
  'chemInchis': async ({crux}: {crux: ChemWasm}, molecules: (string | null)[], keys: boolean) => {
    const lib = await inchiLibrary();
    return crux.chemInchiInputs(molecules).map((input) => inchiOf(lib, input ?? '', keys));
  },
  /** Map Identifiers' InChI to SMILES (python-free.ts): crux-sketch's `inchiToMol`'s SMILES, '' for an empty InChI and
   * one the library refuses (a trap fails the job). */
  'inchiToSmiles': async (_: unknown, inchis: (string | null)[]) => {
    const lib = await inchiLibrary();
    return inchis.map((inchi) => {
      if (!inchi)
        return '';
      try {
        return lib.inchiToMol(inchi).smiles;
      } catch (e) {
        throwIfTrap(e);
        return '';
      }
    });
  },
  /** Map Identifiers' InChI to InChIKey (python-free.ts): null for an empty InChI and where the library gives none. */
  'inchiToInchiKey': async (_: unknown, inchis: (string | null)[]) => {
    const lib = await inchiLibrary();
    return inchis.map((inchi) => inchi ? lib.inchiToInchiKey(inchi) || null : null);
  },
};
