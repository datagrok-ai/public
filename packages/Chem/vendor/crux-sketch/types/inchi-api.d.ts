/**
 * The InChI of InChI's input, as Datagrok Chem's chem-wasm fills it from a molecule (`mol.get_inchi_input(options)`, RDKit's
 * `MolToInchi` up to the library call): the page's InChI module runs the InChI library on it, so the InChI is RDKit.js's
 * `mol.get_inchi(options)` of that molecule, the options having travelled in the input. The InChI module loads on the first
 * call, once per page.
 *
 * @param input - the bytes `get_inchi_input` gave; `null` or `undefined` (what `get_inchi_input` gives where RDKit writes
 *   `''` without asking the library) is `''`, and loads nothing.
 * @returns `''` where RDKit.js's `get_inchi` gives `''`: an input the library writes no InChI of (a dummy atom, a library
 *   error), and `null` or `undefined`.
 * @throws a `SketchError` (the promise rejects), `format: 'inchi'`: for bytes that are not such an input (empty, or not in
 *   its format); `code: 'unavailable'` when the InChI module could not be loaded (a later call tries again); the module's
 *   failure (a trap), after which the module is called no more.
 */
export declare function inchiFromInput(input: Uint8Array | null | undefined): Promise<string>;
/**
 * The InChIKey of an InChI: RDKit's `InchiToInchiKey`, RDKit.js's `get_inchikey_for_inchi`, by the page's InChI module,
 * loaded on the first call that needs it.
 *
 * @returns `''` for `''` (which loads nothing) and for an InChI the library refuses.
 * @throws a `SketchError` (the promise rejects) as `inchiFromInput` does when the module cannot be loaded or has failed.
 */
export declare function inchiKeyFromInchi(inchi: string): Promise<string>;
