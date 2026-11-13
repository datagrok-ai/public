/* tslint:disable */
/* eslint-disable */

/**
 * A CSV held in this wasm instance's memory (copied in once; a leading BOM dropped).
 */
export class CsvInput {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Pass 2 of verified chunk `k` (on a worker) → an Arrow IPC stream of its batches.
     */
    buildChunk(options: any, verified: any, k: number): Uint8Array;
    /**
     * Pass 1 of chunk `k` (on a worker).
     */
    inferChunk(options: any, plan: any, k: number): any;
    constructor(bytes: Uint8Array);
    /**
     * Header + chunk guesses (`null` for an input without rows). `engine`: "Fast" (default) or
     * "Compat" (after `verify` returned `null`).
     */
    plan(options: any, n_chunks: number, engine?: string | null): any;
    /**
     * Checks the chunks (re-running misplaced ones here) and merges the types: the verified plan
     * for `buildChunk`, or `null` (non-standard quoting: plan again with "Compat").
     */
    verify(options: any, plan: any, infos: any): any;
    /**
     * Bytes held (after the BOM).
     */
    readonly length: number;
}

/**
 * The result of [`descriptors_table`]: the requested columns (in request order)
 * and the per-row `error`, as a `RecordBatch`.
 */
export class DescriptorTable {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Column `i` as numbers (null → NaN; see [`DescriptorTable::nulls`]).
     */
    column(i: number): Float64Array;
    /**
     * Per row: why it has no values (`""` when it has).
     */
    errors(): string[];
    /**
     * The descriptor column names, in request order (without `error`).
     */
    names(): string[];
    /**
     * Column `i`'s null mask (1 = null: RDKit has no value, or the row failed).
     */
    nulls(i: number): Uint8Array;
    /**
     * The table as Arrow IPC stream bytes.
     */
    toIpc(): Uint8Array;
    /**
     * `"double"` / `"bigint"` per column.
     */
    types(): string[];
    readonly numRows: number;
}

/**
 * RDKit-js `JSMol`.
 */
export class JSMol {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * PORT(minilib.cpp `add_hs`): the MOL block of a copy with its Hs as atoms
     * (`MolOps::addHs(mol, false, addCoords=true)`).
     */
    add_hs(): string;
    /**
     * PORT(minilib.cpp `add_hs_in_place`): `addHs` in place, the stereo perceived again
     * ([`perceived_again`]); `true`.
     */
    add_hs_in_place(): boolean;
    /**
     * PORT(common.h:894-904): `canonical_tautomer(details)` —
     * `TautomerEnumerator(params).canonicalize`.
     */
    canonical_tautomer(details?: string | null): boolean;
    /**
     * PORT(common.h:918-934): `charge_parent(details)` — `chargeParent(mol, params,
     * skipStandardize)`.
     */
    charge_parent(details?: string | null): boolean;
    /**
     * PORT(common.h:873-880): `cleanup(details)` — `MolStandardize::cleanup`.
     */
    cleanup(details?: string | null): boolean;
    /**
     * `clear_prop(key)`: `true` when it was set.
     */
    clear_prop(key: string): boolean;
    /**
     * PORT(minilib.cpp `convert_to_aromatic_form`): `updatePropertyCache`, `setAromaticity`.
     */
    convert_to_aromatic_form(): void;
    /**
     * PORT(minilib.cpp `convert_to_kekule_form`): `MolOps::Kekulize` in place (its canonical
     * traversal, `rk::mol_ops::kekulize_canonical`).
     */
    convert_to_kekule_form(): void;
    /**
     * `copy()`: a copy of the molecule.
     */
    copy(): JSMol;
    /**
     * Embind's `delete()`: frees the molecule (the same as wasm-bindgen's `free()`).
     */
    delete(): void;
    /**
     * Crux extension: `MetalDisconnector(MetalDisconnectorOptions).disconnect`.
     */
    disconnect_metals(details?: string | null): boolean;
    /**
     * Crux extension: `disconnectOrganometallics(mol, MetalDisconnectorOptions)`
     * (the Syngenta defaults, overridden by the details).
     */
    disconnect_organometallics(details?: string | null): boolean;
    /**
     * PORT(jswrapper.cpp:83-86): `draw_to_canvas(canvas, width, height)`
     * (-1: the canvas's).
     */
    draw_to_canvas(canvas: any, width: number, height: number): string;
    /**
     * PORT(jswrapper.cpp:88-96): `draw_to_canvas_with_highlights(canvas,
     * details)`.
     */
    draw_to_canvas_with_highlights(canvas: any, details: string): string;
    /**
     * PORT(jswrapper.cpp:67-81): `draw_to_canvas_with_offset`.
     */
    draw_to_canvas_with_offset(canvas: any, offsetx: number, offsety: number, width: number, height: number): string;
    /**
     * PORT(common.h:936-952): `fragment_parent(details)` — `fragmentParent(mol,
     * params, skipStandardize)`.
     */
    fragment_parent(details?: string | null): boolean;
    /**
     * PORT(common.h `generate_aligned_coords`): the molecule laid out with the part matching
     * `template` at the template's coordinates (`generateDepictionMatching2DStructure`):
     * `allowRGroups`, `acceptFailure`, `alignOnly`, `referenceSmarts`, `useCoordGen`
     * (`RDDepict::preferCoordGen` for the call: CoordGen, decision `coordgen`). The match as
     * `{"atoms":[...],"bonds":[...]}`; `""` when the template has no coordinates or nothing
     * matches (the coordinates kept); `"{}"` when `acceptFailure` laid it out freely. A layout
     * CoordGen cannot make yet (a part not ported) throws "not supported by Crux: CoordGen ...",
     * where RDKit-js would have laid it out.
     */
    generate_aligned_coords(template: JSMol, details?: string | null): string;
    /**
     * PORT(minilib.cpp `get_aromatic_form`): the MOL block of the aromatic form.
     */
    get_aromatic_form(): string;
    /**
     * MinimalLib `get_atom_pair_fp(details)`.
     */
    get_atom_pair_fp(details?: string | null): string;
    get_atom_pair_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * `get_autocorr2d()`: `rdMolDescriptors.CalcAUTOCORR2D(mol)` (192 values,
     * `AUTOCORR2D_1..192`).
     */
    get_autocorr2d(): Float64Array;
    /**
     * MinimalLib `get_avalon_fp(details)` (`nBits` 512, `avalonSSSBits`, as there); throws in
     * chem-wasm's one file ([`AVALON_AND_MAPS`]).
     */
    get_avalon_fp(details?: string | null): string;
    get_avalon_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * PORT(SmartsWrite.cpp `MolToCXSmarts`): `get_cxsmarts(details)`: the SMARTS, then
     * `getCXExtensions(mol, CX_ALL)` over the SMARTS writer's atom and bond output order.
     */
    get_cxsmarts(details?: string | null): string;
    /**
     * PORT(minilib.cpp:104-110): `get_cxsmiles(details)`: `MolToCXSmiles(mol, params,
     * flags, restoreBondDirs)`; the coordinates written when the molecule has them, as RDKit.
     */
    get_cxsmiles(details?: string | null): string;
    /**
     * PORT(common.h `get_descriptors`): RDKit's `Descriptors::Properties` as JSON.
     */
    get_descriptors(): string;
    /**
     * `get_descriptors_3d()`: `Descriptors3D.CalcMolDescriptors3D(mol)` on the
     * molecule's conformer, in [`descriptors_3d_names`] order. Throws without one.
     */
    get_descriptors_3d(): Float64Array;
    /**
     * `get_descriptors_3d_vectors()`: `CalcAUTOCORR3D` (80), `CalcRDF` (210),
     * `CalcMORSE` (224) and `CalcWHIM` (114, `thresh = 0.001`) on the molecule's
     * conformer, concatenated (628 values). Throws without one.
     */
    get_descriptors_3d_vectors(): Float64Array;
    /**
     * `get_efgs_json()`: the molecule's functional groups (module docs).
     */
    get_efgs_json(): string;
    /**
     * `get_efgs_svg(width, height)`: the molecule drawn with its functional
     * groups coloured as in Ertl's paper (EFGs `col_mol`; 300 x 250 by default).
     */
    get_efgs_svg(width?: number | null, height?: number | null): string;
    /**
     * Crux: any kind (`maccs|rdkit|tt|pattern|morgan|ap|avalon`) in any form, as JSON.
     */
    get_fp_json(kind: string, details?: string | null): string;
    /**
     * PORT(jswrapper.cpp `get_frags_helper`): `{molList, mappings}`: the fragments
     * (`getMolFrags(asMols, sanitizeFrags)`) and `{"frags":[...],"fragsMolAtomMapping":[[...]]}`.
     * `sanitizeFrags` (default true); `copyConformers` only on.
     */
    get_frags(details?: string | null): any;
    /**
     * `get_gasteiger_charges()`: `ComputeGasteigerCharges(mol)` (12 iterations),
     * each atom's `_GasteigerCharge`.
     */
    get_gasteiger_charges(): Float64Array;
    /**
     * `get_gasteiger_h_charges()`: each atom's `_GasteigerHCharge` (the summed
     * charge of its implicit hydrogens).
     */
    get_gasteiger_h_charges(): Float64Array;
    /**
     * `get_inchi(options)`: InChI lives in inchi-wasm (`molToInchi` of `get_molblock()`), the lazy
     * module Chem already ships (decision `chem-wasm`): throws `not supported by Crux`.
     */
    get_inchi(_options?: string | null): string;
    /**
     * PORT(minilib.cpp:158-160): `get_json()`: `MolInterchange::MolToJSONData(mol)`, RDKit's
     * commonchem writer with its extensions (`molfile::mol_interchange`, decision `chem-wasm`).
     */
    get_json(): string;
    /**
     * PORT(minilib.cpp `get_kekule_form`): the MOL block of the Kekulé form.
     */
    get_kekule_form(): string;
    /**
     * MinimalLib `get_maccs_fp()`: the 167 MACCS keys as a bit string.
     */
    get_maccs_fp(): string;
    get_maccs_fp_as_uint8array(): Uint8Array;
    /**
     * `get_molblock(details)`: `MolToMolBlock`; with MinimalLib's details
     * (`includeStereo`, `kekulize`, `useMolBlockWedging`, `addChiralHs`; decision
     * `chem-wasm`: [`crate::minilib::molblock_with_details`]).
     */
    get_molblock(details?: string | null): string;
    /**
     * `rdMolHash.MolHash(mol, func, useCXSmiles, cxFlagsToSkip)`.
     */
    get_molhash(func: string, use_cx_smiles: boolean, cx_flags_to_skip: number): string;
    /**
     * MinimalLib `get_morgan_fp(details)` (radius 2 by default, as there).
     */
    get_morgan_fp(details?: string | null): string;
    get_morgan_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * PORT(minilib.cpp:396-409): `get_new_coords(useCoordGen)`: the mol block of a copy
     * with new 2D coordinates, `RDDepict::preferCoordGen` = `useCoordGen` for the call.
     */
    get_new_coords(use_coord_gen?: boolean | null): string;
    /**
     * `get_num_atoms(heavyOnly)`.
     */
    get_num_atoms(heavy_only?: boolean | null): number;
    /**
     * `get_num_bonds()`.
     */
    get_num_bonds(): number;
    /**
     * MinimalLib `get_pattern_fp(details)`.
     */
    get_pattern_fp(details?: string | null): string;
    get_pattern_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * `get_prop(key)`: the molecule property as text (`""` when unset).
     */
    get_prop(key: string): string;
    /**
     * `get_prop_list(includePrivate, includeComputed)`: the molecule property names (an
     * array; private names, `_`-prefixed, only with `includePrivate`).
     */
    get_prop_list(include_private?: boolean | null, _include_computed?: boolean | null): string[];
    /**
     * `get_rdkit_descriptors()`: `Descriptors.CalcMolDescriptors(mol)` values in
     * [`rdkit_descriptor_names`] order (NaN where RDKit has none).
     */
    get_rdkit_descriptors(): Float64Array;
    /**
     * MinimalLib `get_rdkit_fp(details)`.
     */
    get_rdkit_fp(details?: string | null): string;
    get_rdkit_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * PORT(minilib.cpp:112-118): `get_smarts(details)`: `MolToSmarts(mol, params)`
     * (`doIsomericSmiles`, `rootedAtAtom`).
     */
    get_smarts(details?: string | null): string;
    /**
     * `get_smiles(details)`: RDKit's canonical SMILES; with MinimalLib's details
     * (`SmilesWriteParams` JSON, decision `chem-wasm`: [`crate::minilib::smiles_write_params`]).
     */
    get_smiles(details?: string | null): string;
    /**
     * PORT(minilib.cpp `get_stereo_tags`): `{"CIP_atoms":[[atom,"(R)"],...],"CIP_bonds":[[begin,end,"(E)"],...]}`:
     * `assignStereochemistry(cleanIt, force, flagPossible)`, then the new CIP labeler
     * (`CIPLabeler::assignCIPLabels`, crate cip-labeler) on that molecule; a possible centre it
     * leaves unlabelled is `(?)`.
     */
    get_stereo_tags(): string;
    /**
     * PORT(minilib.cpp:189-209): `get_substruct_match(qmol, details)`: the first match as
     * `{"atoms":[...],"bonds":[...]}`, `"{}"` for none.
     */
    get_substruct_match(q: JSMol, details?: string | null): string;
    /**
     * PORT(minilib.cpp:211-233): `get_substruct_matches(qmol, details)`: every match
     * (uniquified, at most 1000) as a JSON array, `"{}"` for none.
     */
    get_substruct_matches(q: JSMol, details?: string | null): string;
    /**
     * RDKit-js `get_substruct_matches(qmol)` with the query given as SMARTS
     * (there is no query `JSMol` yet): JSON
     * `[{"atoms":[...],"bonds":[...]}, ...]`, atoms in query-atom order, bonds
     * in query-bond order (PORT(minilib.cpp:189-209) `get_substruct_matches`).
     */
    get_substruct_matches_smarts(smarts: string): string;
    /**
     * PORT(minilib.cpp:130-132): `get_svg(width, height)` (250 x 200 by
     * default).
     */
    get_svg(width?: number | null, height?: number | null): string;
    /**
     * PORT(minilib.cpp:133-138): `get_svg_with_highlights(details)`.
     */
    get_svg_with_highlights(details: string): string;
    /**
     * MinimalLib `get_topological_torsion_fp(details)`.
     */
    get_topological_torsion_fp(details?: string | null): string;
    get_topological_torsion_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * `get_v2Kmolblock(details)`: `MolToV2KMolBlock` (an error when V2000 cannot hold it).
     */
    get_v2Kmolblock(details?: string | null): string;
    /**
     * PORT(minilib.cpp:150-153): `get_v3Kmolblock(details)`: `MolToV3KMolBlock`.
     */
    get_v3Kmolblock(details?: string | null): string;
    /**
     * PORT(minilib.cpp:518-523): `has_coords()`: 0, 2 or 3.
     */
    has_coords(): number;
    /**
     * `has_prop(key)`: the molecule property is set.
     */
    has_prop(key: string): boolean;
    /**
     * `is_valid()`.
     */
    is_valid(): boolean;
    /**
     * Crux extension: `isotopeParent(mol, params, skipStandardize)`.
     */
    isotope_parent(details?: string | null): boolean;
    /**
     * Crux extension: `LargestFragmentChooser(preferOrganic, useAtomCount,
     * countHeavyAtomsOnly).choose` (defaults from the params:
     * `preferOrganic`, `largestFragmentChooserUseAtomCount`,
     * `largestFragmentChooserCountHeavyAtomsOnly`).
     */
    largest_fragment(details?: string | null): boolean;
    /**
     * PORT(common.h:906-916): `neutralize(details)` —
     * `Uncharger(params.doCanonical).uncharge`.
     */
    neutralize(details?: string | null): boolean;
    /**
     * PORT(common.h:882-886): `normalize(details)` — `MolStandardize::normalize`.
     */
    normalize(details?: string | null): boolean;
    /**
     * PORT(minilib.cpp `normalize_depiction`): RDKit's `normalizeDepiction(mol, -1,
     * canonicalize, scaleFactor)` with MinimalLib's canonical transform (its build without
     * Eigen: `depictor::CanonicalTransformBuild::PowerSolver`); the scale factor, or -1
     * (nothing done). Throws where MinimalLib's does (no principal axis).
     */
    normalize_depiction(canonicalize?: number | null, scale_factor?: number | null): number;
    /**
     * Crux extension: `normalizeInPlace(mol, params)` — RDKit's in-place Normalizer
     * (`runReactant` in place), which can give another result than `normalize`.
     */
    normalize_in_place(details?: string | null): boolean;
    /**
     * PORT(common.h:888-892): `reionize(details)` — `MolStandardize::reionize`.
     */
    reionize(details?: string | null): boolean;
    /**
     * Crux extension: `removeFragments(mol, params)`; with `leave_last` /
     * `skip_if_all_match` in the details, `FragmentRemover(params' catalog,
     * leave_last, skip_if_all_match).remove` (defaults true / false).
     */
    remove_fragments(details?: string | null): boolean;
    /**
     * PORT(minilib.cpp `remove_hs`): the MOL block of a copy without its Hs: RDKit's
     * `removeAllHs` (isotopic Hs too, as RDKit-js), or with details `removeHs` with RDKit's
     * defaults (details only `{}`). A copy without coordinates is laid out with the atom ranks
     * the removal kept (RDKit's depictor reads them); the writer then reads the molecule without
     * its computed stereo ([`without_computed_stereo`]), so a potential stereo double bond is
     * written "either" (3), as RDKit-js's `remove_hs()` writes it.
     */
    remove_hs(details?: string | null): string;
    /**
     * PORT(minilib.cpp `remove_hs_in_place`): [`Self::remove_hs`]'s removal in place, the stereo
     * perceived again ([`perceived_again`]); `true`.
     */
    remove_hs_in_place(details?: string | null): boolean;
    /**
     * PORT(minilib.cpp:383-394): `set_new_coords(useCoordGen)`: `compute2DCoords` (C++
     * defaults: no canonical orientation) with `RDDepict::preferCoordGen` = `useCoordGen` for
     * the call ([`with_prefer_coord_gen`]: CoordGen, crate `coordgen`, decision `coordgen`),
     * under the stereo perception RDKit-js lays out with ([`compute_new_coords`]).
     */
    set_new_coords(use_coord_gen?: boolean | null): boolean;
    /**
     * `set_prop(key, value, computed)`: sets a molecule property (text; `_MolFileChiralFlag`
     * as the number the MOL writer reads); `true`.
     */
    set_prop(key: string, value: string, _computed?: boolean | null): boolean;
    /**
     * Crux extension: `stereoParent(mol, params, skipStandardize)`.
     */
    stereo_parent(details?: string | null): boolean;
    /**
     * PORT(minilib.cpp `straighten_depiction`): RDKit's `straightenDepiction(mol, -1,
     * minimizeRotation)`.
     */
    straighten_depiction(minimize_rotation?: boolean | null): void;
    /**
     * Crux extension: `superParent(mol, params, skipStandardize)`.
     */
    super_parent(details?: string | null): boolean;
    /**
     * Crux extension: `tautomerParent(mol, params, skipStandardize)`.
     */
    tautomer_parent(details?: string | null): boolean;
    /**
     * Crux extension: `Uncharger(canonicalOrder, force, protonationOnly).uncharge`
     * (`canonicalOrder` defaults to the params' `doCanonical`).
     */
    uncharge(details?: string | null): boolean;
    /**
     * Crux extension: `ValidationMethod::validate(mol, reportAllFailures)` — a JSON
     * array of the messages. Details: `validations` (RDKit class names, run as a
     * `CompositeValidation`; default `["RDKitValidation"]`), `reportAllFailures`
     * (default false), and the method options (`allowEmptyMolecules`, `strict`,
     * `allowedAtoms` / `disallowedAtoms` as atomic numbers, the
     * `FeaturesValidation` `allow*` flags). RDKitValidation reports valence errors
     * only on an unsanitized molecule: `get_mol(smiles, '{"sanitize":false}')`.
     */
    validate(details?: string | null): string;
}

/**
 * `MAP4Calculator` (module docs).
 */
export class Map4Calculator {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * `calculate(mol)`.
     */
    calculateMol(mol: JSMol): Uint32Array;
    /**
     * `calculate(Chem.MolFromSmiles(smiles))`.
     */
    calculate(smiles: string): Uint32Array;
    constructor(options?: string | null);
}

/**
 * `MHFPEncoder` (module docs).
 */
export class MhfpEncoder {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * `encode_mol(in_mol, …)`.
     */
    encodeMol(mol: JSMol, options?: string | null): Uint32Array;
    /**
     * `encode(in_smiles, radius, rings, kekulize, min_radius, sanitize=True)`.
     */
    encode(smiles: string, options?: string | null): Uint32Array;
    /**
     * `{"nPermutations": 2048, "seed": 42}`.
     */
    constructor(options?: string | null);
}

/**
 * RDKit-js `MolList` (MinimalLib `JSMolList`). An entry can be empty (`null`, as RDKit-js):
 * `get_mmpa_frags`'s core of a single cut.
 */
export class MolList {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * PORT(minilib.cpp:697-700): `append(mol)` — a copy of `mol`; the new size.
     */
    append(mol: JSMol): number;
    /**
     * PORT(minilib.cpp:674-683): `at(idx)`; `null` past the end or for an empty entry.
     */
    at(idx: number): JSMol | null;
    /**
     * `at_end()`.
     */
    at_end(): boolean;
    /**
     * Embind's `delete()` (wasm-bindgen's `free()`).
     */
    delete(): void;
    /**
     * PORT(minilib.cpp:702-706): `insert(idx, mol)` — at `min(idx, size)`; the new
     * size.
     */
    insert(idx: number, mol: JSMol): number;
    /**
     * `new MolList()`.
     */
    constructor();
    /**
     * PORT(minilib.cpp:666-672): `next()` — the molecule at the cursor, then the
     * cursor moves on; `null` at the end (or for an empty entry), as RDKit-js.
     */
    next(): JSMol | null;
    /**
     * PORT(minilib.cpp:685-695): `pop(idx)` — removes and returns the molecule (`null` past the
     * end); the cursor stays on the same next molecule.
     */
    pop(idx: number): JSMol | null;
    /**
     * `reset()` — the cursor back to the start.
     */
    reset(): void;
    /**
     * `size()`.
     */
    size(): number;
}

/**
 * [`info`] as JSON text: `tools/parity/wasm-modules/modules.py build` writes the loader's
 * manifest.json (function → module) from the modules it built.
 */
export function cruxModuleInfo(): string;

/**
 * The modules this wasm was built with, comma-separated (`"ml"`).
 */
export function cruxModules(): string;

/**
 * The whole file on this thread → an Arrow IPC stream.
 */
export function csvToArrow(bytes: Uint8Array, options: any): Uint8Array;

/**
 * `descriptorNames()`: every name `descriptorsTable` accepts, in order:
 * `[{name, family, type: "double" | "bigint", needs3d, doc}]` — for pickers.
 */
export function descriptorNames(): any;

/**
 * `descriptors3dNames()`: `Descriptors3D.descList` names, in order.
 */
export function descriptors3dNames(): string[];

/**
 * `descriptorsArrow(smiles, names, seed?)`: [`descriptors_table`] as Arrow IPC
 * stream bytes (`application/vnd.apache.arrow.stream`): the requested columns
 * (Float64; Int64 for counts; null where RDKit has none), then `error`.
 */
export function descriptorsArrow(smiles: string[], names: string[], seed?: number | null): Uint8Array;

/**
 * `descriptorsTable(smiles, names, seed?)`: RDKit descriptors by name
 * ([`descriptor_names_info`]; `all`, `@<family>`, `@all`) over a SMILES list —
 * only the requested ones are computed. Throws on an unknown name (with the near
 * matches). 3D names embed `AddHs(mol)` with ETKDGv3 (`seed`, default 42).
 */
export function descriptorsTable(smiles: string[], names: string[], seed?: number | null): DescriptorTable;

/**
 * Crux extension: `enumerateTautomerSmiles(smiles, params)` — a JSON array of the
 * tautomer SMILES.
 */
export function enumerate_tautomer_smiles(smiles: string, details?: string | null): string;

/**
 * The entries of a catalog selection, in match order (`entry` is the index
 * `structuralAlerts` reports).
 */
export function filterCatalogEntries(catalogs?: string | null): any;

/**
 * Register resource `name` as a dataset read in place: `manifestJson` (its
 * `crux-dataset.json`) and `base`: a base URL (`<base>/<path>` per file, HTTP range), or an
 * object mapping each file's relative path to a source (an OPFS sync access handle, a
 * `File` / `Blob`, bytes, a URL, or a synchronous `{size, read(offset, length)}` reader, e.g.
 * node's `fs.readSync`). Reads block: call functions that read it in a Web Worker.
 * Returns the dataset id (`sha256:` of the canonical manifest), its cache identity.
 */
export function fnDataset(name: string, manifest_json: string, base: any): string;

/**
 * One function's spec (JSON).
 */
export function fnDescribe(name: string): string;

/**
 * Every registered function's spec: `{"crux_fn_version", "functions": [...]}` (JSON).
 */
export function fnList(): string;

/**
 * Register resource `name` (e.g. `model`) as in-memory bytes for later
 * `callFn`s; `id` identifies them in the loaded-object cache (a URL or hash).
 */
export function fnResource(name: string, id: string, bytes: Uint8Array): void;

/**
 * Register resource `name` as one object read in place (`src`: a URL; in a worker also an
 * OPFS sync access handle, a `File` / `Blob` or bytes); `id` identifies it in caches.
 */
export function fnSource(name: string, id: string, src: any): void;

/**
 * Generic Murcko scaffold of one molecule (all atoms C, all bonds single) —
 * `MolToSmiles(MakeScaffoldGeneric(GetScaffoldForMol(MolFromSmiles(s))))`.
 */
export function genericScaffold(smiles: string): string;

/**
 * The functional groups of a column of SMILES in one call: an array of the
 * same length holding each molecule's JSON (as [`JSMol::get_efgs_json`]), or
 * `null` where it does not parse.
 */
export function getEfgsJson(smiles: string[]): any;

/**
 * RDKit-js `get_inchikey_for_inchi(inchi)`: in inchi-wasm (`inchiToInchiKey`): throws.
 */
export function get_inchikey_for_inchi(_inchi: string): string;

/**
 * RDKit-js `get_mol(input, details)` (decision `chem-wasm`): `null` when the input does
 * not parse, as RDKit-js 2026.03.6; throws for details that are not JSON or that ask
 * what Crux cannot do (`not supported by Crux: <key>`). Details keys: [`MolInputDetails`].
 */
export function get_mol(input: string, details?: string | null): JSMol | null;

/**
 * PORT(minilib.cpp:827-829): `get_mol_copy(other)`.
 */
export function get_mol_copy(other: JSMol): JSMol;

/**
 * PORT(common.h `qmol_from_input`): RDKit-js `get_qmol(input)`: a SMARTS (`SmartsToMol`) or
 * a MOL block read as a query (`MolBlockToMol(sanitize=false, removeHs=false)`, its query
 * atoms and bonds); `null` when it does not parse.
 */
export function get_qmol(input: string): JSMol | null;

/**
 * The identity key of `smiles` under the named policy (`molstandardize::identity`).
 */
export function identityKey(smiles: string, policy: string): string;

/**
 * [`identity_key`] of a column: `null` where a molecule fails.
 */
export function identityKeys(smiles: string[], policy: string): any;

/**
 * The identity policies: `[{name, description}]`.
 */
export function identityPolicies(): any;

/**
 * Load the NP_Score model (`npscore_publicnp.cxft` bytes). Once per page;
 * later calls are no-ops.
 */
export function initNpScore(bytes: Uint8Array): void;

/**
 * Load the SA_Score fragment table (`sascore_fpscores.cxft` bytes). Once per
 * page; later calls are no-ops.
 */
export function initSaScore(bytes: Uint8Array): void;

export function leafAbi(): number;

export function leafCall(name: string, table: any, params: any, opts_json: string | null | undefined, keep: any, call_id?: number | null, on_progress?: Function | null, on_log?: Function | null): any;

export function leafCallAsync(name: string, table: any, params: any, opts_json: string | null | undefined, keep: any, call_id?: number | null, on_progress?: Function | null, on_log?: Function | null): Promise<any>;

export function leafCancel(call_id: number): boolean;

/**
 * `{instance, cap}` of the handle table, set from `cfg` first where given (`instance`: the
 * host's number for this leaf instance, before any handle; `cap`: kept bytes).
 */
export function leafConfig(cfg: any): any;

export function leafDescribe(handle: number): any;

export function leafFree(handle: number): boolean;

export function leafGather(spec: any): any;

/**
 * `{abi, leaf, functions, probes}`: the module's functions (specs: its fns.json) and, per
 * variant with a device probe, whether it runs here.
 */
export function leafInfo(): any;

/**
 * Pins (`on`) or unpins a kept result: a pinned one is never evicted (a run's live data).
 */
export function leafPin(handle: number, on: boolean): boolean;

/**
 * A call checked without running it: `{fn, want, shape}` (crux_leaf::abi::plan).
 */
export function leafPlan(name: string, table: any, params: any, opts_json: string | null | undefined, resources: any, relax_variant?: boolean | null): any;

export function leafPut(table: any, keep: any): any;

export function leafRead(handle: number, view: boolean, which?: string | null): any;

/**
 * Ends lease `lease` (its stored resources evictable again); a call naming it ends it too.
 */
export function leafRelease(lease: number): boolean;

/**
 * `fnResource` under the ABI's name.
 */
export function leafResource(name: string, id: string, bytes: Uint8Array): void;

export function leafStats(): any;

/**
 * Stores resource bytes by identity: calls bind them by `resources: {name: {id}}` only.
 */
export function leafStore(id: string, bytes: Uint8Array, lease?: number | null): void;

/**
 * Which of the resource identities `ids` are not stored here (evicted, unstored, never stored);
 * the others marked used.
 */
export function leafStored(ids: any, lease?: number | null): any;

/**
 * Drops stored resource bytes `id`.
 */
export function leafUnstore(id: string): boolean;

/**
 * Load RDKit Contrib's `SubstructureFilter_HitTriaging_wPubChemExamples.csv`
 * (its text, fetched by the caller; nothing is built in). Returns the number
 * of filters; replaces a previously loaded table.
 */
export function loadNibrCatalog(csv_text: string): number;

/**
 * tmap `Minhash.get_distance(a, b)` (MAP4).
 */
export function map4Distance(a: Uint32Array, b: Uint32Array): number;

/**
 * `MHFPEncoder.distance(a, b)`.
 */
export function mhfpDistance(a: Uint32Array, b: Uint32Array): number;

/**
 * `rdMolHash.MolHash(Chem.MolFromSmiles(smiles), func, useCXSmiles, cxFlagsToSkip)`;
 * throws on a SMILES that does not read.
 */
export function molHash(smiles: string, func: string, use_cx_smiles: boolean, cx_flags_to_skip: number): string;

/**
 * The MolHash function names, RDKit's order.
 */
export function molHashNames(): string[];

/**
 * [`mol_hash`] of a column: an array of the same length, `null` where a molecule fails.
 */
export function molHashes(smiles: string[], func: string, use_cx_smiles: boolean, cx_flags_to_skip: number): any;

/**
 * Bemis–Murcko scaffold of one molecule — RDKit
 * `MurckoScaffold.MurckoScaffoldSmiles(smiles, includeChirality)`; `""` for an
 * acyclic molecule. `isomeric = false` is RDKit's default (no stereo, no
 * isotopes). Throws on a SMILES that does not parse.
 */
export function murckoScaffold(smiles: string, isomeric: boolean): string;

/**
 * Murcko (or, with `generic`, generic) scaffolds of a column of SMILES in one
 * call: an array of the same length, `null` where a molecule fails.
 */
export function murckoScaffolds(smiles: string[], generic: boolean, isomeric: boolean): any;

/**
 * RDKit Contrib `assignFilters` for one molecule, plus each matching NIBR
 * filter with its annotations and match atoms. Needs `loadNibrCatalog` first.
 */
export function nibrFilters(smiles: string): any;

/**
 * NP-likeness with confidence — `npscorer.scoreMolWConfidence`:
 * `{score, confidence, rawScore, fragmentsFound, fragmentsTotal}`, or `null`
 * for a molecule without atoms. Throws on a SMILES that does not parse, or
 * before [`init_np_score`].
 */
export function npScore(smiles: string): any;

/**
 * NP scores of a column of SMILES: an array of `{score, confidence, …}` or
 * `null`. Throws before [`init_np_score`].
 */
export function npScores(smiles: string[]): any;

/**
 * PORT(minilib.cpp:864-868): RDKit-js `prefer_coordgen(prefer)`: `RDDepict::preferCoordGen`,
 * process-wide (CoordGen, crate `coordgen`, decision `coordgen`); `true` throws in a build without
 * CoordGen.
 */
export function prefer_coordgen(prefer: boolean): void;

/**
 * QED of one molecule: the three weightings, the eight properties and each
 * property's desirability.
 */
export function qed(smiles: string): any;

/**
 * QED (mean weights) per SMILES; `NaN` where a SMILES does not parse.
 */
export function qedMany(smiles: string[]): Float64Array;

/**
 * `rdkitDescriptorNames()`: `Descriptors._descList` names, in order.
 */
export function rdkitDescriptorNames(): string[];

/**
 * `RegistrationHash.GetMolLayers(mol, enable_tautomer_hash_v2)` and `GetMolHash(layers,
 * scheme)` (`ALL_LAYERS`, `STEREO_INSENSITIVE_LAYERS`, `TAUTOMER_INSENSITIVE_LAYERS`) as
 * `{"layers": {"CANONICAL_SMILES": …, …}, "hash": "<sha1>"}`.
 */
export function registrationHash(smiles: string, scheme: string, tautomer_v2: boolean): any;

/**
 * SA score (1 easy … 10 hard) — `sascorer.calculateScore(MolFromSmiles(s))`;
 * `undefined` for a molecule without atoms. Throws on a SMILES that does not
 * parse, or before [`init_sa_score`].
 */
export function saScore(smiles: string): number | undefined;

/**
 * SA score with its terms: `{score, fragmentScore, complexityScore,
 * densityCorrection, numAtoms, numChiralCenters, numSpiroAtoms,
 * numBridgeheadAtoms, numMacrocycles}`, or `null`.
 */
export function saScoreDetails(smiles: string): any;

/**
 * SA scores of a column of SMILES: an array of the same length, `null` where a
 * molecule fails or has no atoms. Throws before [`init_sa_score`].
 */
export function saScores(smiles: string[]): any;

/**
 * `MHFPEncoder.secfp_from_smiles(in_smiles, length, radius, rings, kekulize, sanitize=False)`:
 * the on bits.
 */
export function secfp(smiles: string, options?: string | null): Uint32Array;

/**
 * `MHFPEncoder.secfp_from_mol(in_mol, …)`: the on bits.
 */
export function secfpMol(mol: JSMol, options?: string | null): Uint32Array;

/**
 * RDKit-js `set_log_capture(logName)` / `set_log_tee(logName)`: Crux keeps no RDKit log (its
 * messages are not RDKit's), so these throw `not supported by Crux` rather than being undefined.
 */
export function set_log_capture(log_name: string): any;

/**
 * See [`set_log_capture`].
 */
export function set_log_tee(log_name: string): any;

/**
 * RDKit `SimilarityMaps.GetSimilarityMapForFingerprint(ref, probe, fp, draw2d)`
 * on a `MolDraw2DSVG(width, height)`: the SVG text.
 */
export function similarityMapSvg(reference: JSMol, probe: JSMol, options?: string | null): string;

/**
 * RDKit `SimilarityMaps.GetAtomicWeightsForFingerprint(ref, probe, fp,
 * metric)` as JSON `{"weights":[...],"max_weight":x}` (raw weights; the
 * map draws them divided by `max_weight`).
 */
export function similarityMapWeights(reference: JSMol, probe: JSMol, options?: string | null): string;

/**
 * Crux extension: `standardizeSmiles(smiles)` — the SMILES read unsanitized,
 * `cleanup` with the default parameters, `MolToSmiles`.
 */
export function standardize_smiles(smiles: string): string;

/**
 * Install a panic hook that forwards Rust panics to `console.error`. Called
 * once at module init on wasm; a no-op on native targets.
 */
export function start(): void;

/**
 * Every matching FilterCatalog entry for one molecule, with its matches.
 */
export function structuralAlerts(smiles: string, catalogs?: string | null): any;

/**
 * RDKit-js `use_legacy_stereo_perception(value)`: RDKit's legacy stereo perception on or
 * off (process-wide, legacy by default as in RDKit); the value it had.
 */
export function use_legacy_stereo_perception(value: boolean): boolean;

/**
 * Crux extension: `validateSmiles(smiles)` — a JSON array of the messages.
 */
export function validate_smiles(smiles: string): string;

/**
 * RDKit-js `version()`: the RDKit release whose behaviour Crux reproduces (`RDKIT_VERSION`).
 */
export function version(): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_csvinput_free: (a: number, b: number) => void;
    readonly __wbg_descriptortable_free: (a: number, b: number) => void;
    readonly __wbg_jsmol_free: (a: number, b: number) => void;
    readonly __wbg_map4calculator_free: (a: number, b: number) => void;
    readonly __wbg_mhfpencoder_free: (a: number, b: number) => void;
    readonly __wbg_mollist_free: (a: number, b: number) => void;
    readonly cruxModuleInfo: () => [number, number];
    readonly cruxModules: () => [number, number];
    readonly csvToArrow: (a: number, b: number, c: any) => [number, number, number, number];
    readonly csvinput_buildChunk: (a: number, b: any, c: any, d: number) => [number, number, number, number];
    readonly csvinput_inferChunk: (a: number, b: any, c: any, d: number) => [number, number, number];
    readonly csvinput_length: (a: number) => number;
    readonly csvinput_new: (a: number, b: number) => number;
    readonly csvinput_plan: (a: number, b: any, c: number, d: number, e: number) => [number, number, number];
    readonly csvinput_verify: (a: number, b: any, c: any, d: any) => [number, number, number];
    readonly descriptorNames: () => [number, number, number];
    readonly descriptors3dNames: () => [number, number];
    readonly descriptorsArrow: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly descriptorsTable: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly descriptortable_column: (a: number, b: number) => [number, number];
    readonly descriptortable_errors: (a: number) => [number, number];
    readonly descriptortable_names: (a: number) => [number, number];
    readonly descriptortable_nulls: (a: number, b: number) => [number, number];
    readonly descriptortable_numRows: (a: number) => number;
    readonly descriptortable_toIpc: (a: number) => [number, number, number, number];
    readonly descriptortable_types: (a: number) => [number, number];
    readonly enumerate_tautomer_smiles: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly filterCatalogEntries: (a: number, b: number) => [number, number, number];
    readonly fnDataset: (a: number, b: number, c: number, d: number, e: any) => [number, number, number, number];
    readonly fnDescribe: (a: number, b: number) => [number, number, number, number];
    readonly fnList: () => [number, number];
    readonly fnResource: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
    readonly fnSource: (a: number, b: number, c: number, d: number, e: any) => [number, number];
    readonly genericScaffold: (a: number, b: number) => [number, number, number, number];
    readonly getEfgsJson: (a: number, b: number) => [number, number, number];
    readonly get_inchikey_for_inchi: (a: number, b: number) => [number, number, number, number];
    readonly get_mol: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly get_mol_copy: (a: number) => number;
    readonly get_qmol: (a: number, b: number) => any;
    readonly identityKey: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly identityKeys: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly identityPolicies: () => [number, number, number];
    readonly initNpScore: (a: number, b: number) => [number, number];
    readonly initSaScore: (a: number, b: number) => [number, number];
    readonly jsmol_add_hs: (a: number) => [number, number, number, number];
    readonly jsmol_add_hs_in_place: (a: number) => [number, number, number];
    readonly jsmol_canonical_tautomer: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_charge_parent: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_cleanup: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_clear_prop: (a: number, b: number, c: number) => number;
    readonly jsmol_convert_to_aromatic_form: (a: number) => [number, number];
    readonly jsmol_convert_to_kekule_form: (a: number) => [number, number];
    readonly jsmol_copy: (a: number) => number;
    readonly jsmol_delete: (a: number) => void;
    readonly jsmol_disconnect_metals: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_disconnect_organometallics: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_draw_to_canvas: (a: number, b: any, c: number, d: number) => [number, number, number, number];
    readonly jsmol_draw_to_canvas_with_highlights: (a: number, b: any, c: number, d: number) => [number, number, number, number];
    readonly jsmol_draw_to_canvas_with_offset: (a: number, b: any, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly jsmol_fragment_parent: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_generate_aligned_coords: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly jsmol_get_aromatic_form: (a: number) => [number, number, number, number];
    readonly jsmol_get_atom_pair_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_atom_pair_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_autocorr2d: (a: number) => [number, number];
    readonly jsmol_get_avalon_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_avalon_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_cxsmarts: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_cxsmiles: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_descriptors: (a: number) => [number, number, number, number];
    readonly jsmol_get_descriptors_3d: (a: number) => [number, number, number, number];
    readonly jsmol_get_descriptors_3d_vectors: (a: number) => [number, number, number, number];
    readonly jsmol_get_efgs_json: (a: number) => [number, number, number, number];
    readonly jsmol_get_efgs_svg: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_fp_json: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly jsmol_get_frags: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_get_gasteiger_charges: (a: number) => [number, number];
    readonly jsmol_get_gasteiger_h_charges: (a: number) => [number, number];
    readonly jsmol_get_inchi: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_json: (a: number) => [number, number, number, number];
    readonly jsmol_get_kekule_form: (a: number) => [number, number, number, number];
    readonly jsmol_get_maccs_fp: (a: number) => [number, number, number, number];
    readonly jsmol_get_maccs_fp_as_uint8array: (a: number) => [number, number, number, number];
    readonly jsmol_get_molblock: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_molhash: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly jsmol_get_morgan_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_morgan_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_new_coords: (a: number, b: number) => [number, number, number, number];
    readonly jsmol_get_num_atoms: (a: number, b: number) => number;
    readonly jsmol_get_num_bonds: (a: number) => number;
    readonly jsmol_get_pattern_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_pattern_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_prop: (a: number, b: number, c: number) => [number, number];
    readonly jsmol_get_prop_list: (a: number, b: number, c: number) => [number, number];
    readonly jsmol_get_rdkit_descriptors: (a: number) => [number, number];
    readonly jsmol_get_rdkit_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_rdkit_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_smarts: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_smiles: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_stereo_tags: (a: number) => [number, number, number, number];
    readonly jsmol_get_substruct_match: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly jsmol_get_substruct_matches: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly jsmol_get_substruct_matches_smarts: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_svg: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_svg_with_highlights: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_topological_torsion_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_topological_torsion_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_v2Kmolblock: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_v3Kmolblock: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_has_coords: (a: number) => number;
    readonly jsmol_has_prop: (a: number, b: number, c: number) => number;
    readonly jsmol_is_valid: (a: number) => number;
    readonly jsmol_isotope_parent: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_largest_fragment: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_neutralize: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_normalize: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_normalize_depiction: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly jsmol_normalize_in_place: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_reionize: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_remove_fragments: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_remove_hs: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_remove_hs_in_place: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_set_new_coords: (a: number, b: number) => [number, number, number];
    readonly jsmol_set_prop: (a: number, b: number, c: number, d: number, e: number, f: number) => number;
    readonly jsmol_stereo_parent: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_straighten_depiction: (a: number, b: number) => void;
    readonly jsmol_super_parent: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_tautomer_parent: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_uncharge: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_validate: (a: number, b: number, c: number) => [number, number, number, number];
    readonly leafAbi: () => number;
    readonly leafCall: (a: number, b: number, c: any, d: any, e: number, f: number, g: any, h: number, i: number, j: number) => [number, number, number];
    readonly leafCallAsync: (a: number, b: number, c: any, d: any, e: number, f: number, g: any, h: number, i: number, j: number) => any;
    readonly leafCancel: (a: number) => number;
    readonly leafConfig: (a: any) => any;
    readonly leafDescribe: (a: number) => [number, number, number];
    readonly leafFree: (a: number) => number;
    readonly leafGather: (a: any) => [number, number, number];
    readonly leafInfo: () => any;
    readonly leafPin: (a: number, b: number) => number;
    readonly leafPlan: (a: number, b: number, c: any, d: any, e: number, f: number, g: any, h: number) => [number, number, number];
    readonly leafPut: (a: any, b: any) => [number, number, number];
    readonly leafRead: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly leafRelease: (a: number) => number;
    readonly leafResource: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
    readonly leafStats: () => any;
    readonly leafStore: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly leafStored: (a: any, b: number, c: number) => any;
    readonly leafUnstore: (a: number, b: number) => number;
    readonly loadNibrCatalog: (a: number, b: number) => [number, number, number];
    readonly map4Distance: (a: number, b: number, c: number, d: number) => number;
    readonly map4calculator_calculate: (a: number, b: number, c: number) => [number, number, number, number];
    readonly map4calculator_calculateMol: (a: number, b: number) => [number, number];
    readonly map4calculator_new: (a: number, b: number) => [number, number, number];
    readonly mhfpDistance: (a: number, b: number, c: number, d: number) => number;
    readonly mhfpencoder_encode: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly mhfpencoder_encodeMol: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly mhfpencoder_new: (a: number, b: number) => [number, number, number];
    readonly molHash: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly molHashNames: () => [number, number];
    readonly molHashes: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number];
    readonly mollist_append: (a: number, b: number) => number;
    readonly mollist_at: (a: number, b: number) => any;
    readonly mollist_at_end: (a: number) => number;
    readonly mollist_delete: (a: number) => void;
    readonly mollist_insert: (a: number, b: number, c: number) => number;
    readonly mollist_new: () => number;
    readonly mollist_next: (a: number) => any;
    readonly mollist_pop: (a: number, b: number) => any;
    readonly mollist_reset: (a: number) => void;
    readonly mollist_size: (a: number) => number;
    readonly murckoScaffold: (a: number, b: number, c: number) => [number, number, number, number];
    readonly murckoScaffolds: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly nibrFilters: (a: number, b: number) => [number, number, number];
    readonly npScore: (a: number, b: number) => [number, number, number];
    readonly npScores: (a: number, b: number) => [number, number, number];
    readonly prefer_coordgen: (a: number) => [number, number];
    readonly qed: (a: number, b: number) => [number, number, number];
    readonly qedMany: (a: number, b: number) => [number, number];
    readonly rdkitDescriptorNames: () => [number, number];
    readonly registrationHash: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly saScore: (a: number, b: number) => [number, number, number, number];
    readonly saScoreDetails: (a: number, b: number) => [number, number, number];
    readonly saScores: (a: number, b: number) => [number, number, number];
    readonly secfp: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly secfpMol: (a: number, b: number, c: number) => [number, number, number, number];
    readonly set_log_capture: (a: number, b: number) => [number, number, number];
    readonly set_log_tee: (a: number, b: number) => [number, number, number];
    readonly similarityMapSvg: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly similarityMapWeights: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly standardize_smiles: (a: number, b: number) => [number, number, number, number];
    readonly start: () => void;
    readonly structuralAlerts: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly use_legacy_stereo_perception: (a: number) => number;
    readonly validate_smiles: (a: number, b: number) => [number, number, number, number];
    readonly version: () => [number, number];
    readonly wasm_bindgen__convert__closures_____invoke__h274c5cb09364b90a: (a: number, b: number, c: any, d: any) => void;
    readonly wasm_bindgen__convert__closures_____invoke__h86605f80b5250a21: (a: number, b: number, c: any) => [number, number];
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_destroy_closure: (a: number, b: number) => void;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __externref_drop_slice: (a: number, b: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
