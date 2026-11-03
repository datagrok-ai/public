/**
 * InChI in Crux mode, by the split (crux-core decision chem-wasm-inchi): Crux fills InChI's input from the molecule as
 * RDKit's MolToInchi fills it (chem-wasm's `get_inchi_input`), and crux-sketch's InChI module, the page's (loaded once,
 * shared with the sketchers), runs only the InChI library on it. The InChI is RDKit.js's `get_inchi()`, the key its
 * `get_inchikey_for_inchi`; Chem fetches no file of its own for InChI.
 */
import type * as CruxSketch from '../../vendor/crux-sketch/crux-sketch';
import type {JSMol} from '../../vendor/chem-wasm/chem_wasm.js';
import {ChemicalEngine, engineInEffect} from '../engine/chem-engine';
import {cruxModuleSync, noteCruxError} from '../engine/crux-module';
import {cruxOrRdKitAsync} from '../engine/fallback';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import {cruxReadsCxExtension} from './crux-service';
import {loadCruxSketch} from './crux-sketcher';

/** How many molecules run between two yields to the page. */
const YIELD_EVERY = 1000;

/** An InChI and its key of `input` ('' for ''), by crux-sketch's InChI module. Throws when the module cannot load or
 * fails (the caller's RDKit path runs). */
async function inchiOf(sketch: typeof CruxSketch, input: Uint8Array | '', key: boolean): Promise<string> {
  const inchi = input ? await sketch.inchiFromInput(input) : '';
  return key ? sketch.inchiKeyFromInchi(inchi) : inchi;
}

/** `map` of each item, one after another, yielding to the page every {@link YIELD_EVERY} items: a column's InChIs by
 * crux-sketch's InChI module, which runs on the main thread. */
export async function mapYielding<T, R>(items: readonly T[], map: (item: T) => R | Promise<R>): Promise<R[]> {
  const out: R[] = new Array(items.length);
  for (let i = 0; i < items.length; i++) {
    if (i % YIELD_EVERY === YIELD_EVERY - 1)
      await new Promise((resolve) => setTimeout(resolve, 0));
    out[i] = await map(items[i]);
  }
  return out;
}

/** The InChIs (InChI keys when `keys`) of InChI inputs from the Crux workers; null for a null input (a molecule Crux
 * cannot read: the caller's RDKit pool does it). */
export async function inchisFromInputs(inputs: (Uint8Array | '' | null)[], keys: boolean): Promise<(string | null)[]> {
  const sketch = await loadCruxSketch();
  return mapYielding(inputs, (input) => input === null ? null : inchiOf(sketch, input, keys));
}

/**
 * In Crux mode, a molecule's InChI and InChIKey: InChI's input from the main thread's Crux module, the library by
 * crux-sketch's InChI module. Null otherwise, and for a molecule Crux cannot read or a CXSMILES whose query features
 * Crux's reader skips (RDKit writes them then); a throw (RDKit's own, e.g. a molecule it cannot kekulize) leaves it
 * to RDKit too.
 */
export async function cruxInchi(molecule: string): Promise<{inchi: string, inchiKey: string} | null> {
  if (engineInEffect() !== ChemicalEngine.Crux || !cruxReadsCxExtension(molecule))
    return null;
  return cruxOrRdKitAsync('InChI', async () => {
    const mol = getMolSafe(molecule, {}, cruxModuleSync()!, true, noteCruxError).mol;
    if (!mol)
      return null;
    let input: Uint8Array | undefined;
    try {
      input = (mol as unknown as JSMol).get_inchi_input('');
    } finally {
      mol.delete();
    }
    const sketch = await loadCruxSketch();
    const inchi = await inchiOf(sketch, input ?? '', false);
    return {inchi, inchiKey: await sketch.inchiKeyFromInchi(inchi)};
  }, async () => null);
}
