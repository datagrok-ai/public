/* tslint:disable */
/* eslint-disable */

/**
 * Groups of SMARTS compiled once (Chem's structural alert rule sets, its pharmacophore families), matched against
 * many molecules (`match`; chem-wasm's chem_batch.rs docs).
 */
export class ChemSmartsSet {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * `lastWork()`: the last match's molecule × SMARTS pairs, those skipped by what the SMARTS needs, those run
     * (measurements).
     */
    lastWork(): Float64Array;
    /**
     * `match(molecules)`: `{flags, unread}`: `flags[i * groups + g]` 1 where some SMARTS of group `g` matches molecule
     * `i` (each read as Chem's getMolSafe reads it), `unread` the molecules none of its readings reads.
     */
    match(molecules: (string | null)[]): {flags: Uint8Array, unread: Uint32Array};
    /**
     * `new ChemSmartsSet(smarts, groupSizes, options?)`: group `g` is the next `groupSizes[g]` SMARTS, each read as
     * `get_qmol` reads it; `options`: JSON `{"prefilter": true | false | "and", "searchOrder": bool}` (on by default:
     * both shortcuts change no result; chem_batch.rs docs). Throws when the sizes do not add up to the SMARTS.
     */
    constructor(smarts: string[], group_sizes: Uint32Array, options?: string | null);
    /**
     * `unreadSmarts()`: the SMARTS (indices into the list) `get_qmol` reads nothing from.
     */
    unreadSmarts(): Uint32Array;
    /**
     * The number of groups.
     */
    readonly groups: number;
}

/**
 * A synthon library read for searching (`SynthonSpace.ReadTextFile`); the similarity
 * fingerprints are built at its first similarity search.
 */
export class ChemSynthonSpace {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Part `part` of `parts` of the fingerprint build (`BuildSynthonFingerprints`), as bytes, for
     * the copy of this library on another worker to install: the parts add up to the whole build.
     */
    fingerprintPart(part: number, parts: number): Uint8Array;
    /**
     * Whether the similarity fingerprints are built (at the first similarity search, or by
     * [`Self::install_fingerprints`]).
     */
    hasFingerprints(): boolean;
    /**
     * Installs the fingerprints from every part of [`Self::fingerprint_part`], given as one buffer:
     * each part's bytes after its length (u32, little-endian).
     */
    installFingerprints(parts: Uint8Array): void;
    /**
     * The library from its text (RDKit's synthon text format: the Enamine / ChemSpace CSV);
     * throws on a file RDKit cannot read.
     */
    constructor(text: string);
    /**
     * The number of products the library describes.
     */
    numProducts(): number;
    numReactions(): number;
    numSynthons(): number;
    /**
     * The script on `molecule` (module docs): its dataframe as JSON.
     */
    search(molecule: string, search_type: string, max_hits: number, similarity_cutoff: number, include_synthons: boolean): string;
}

/**
 * A finalised, queryable collection of molecules.
 *
 * `indexed` collections carry fingerprint sidecars and answer queries via the
 * screened matcher / sidecar Tanimoto kernels. Non-indexed ("direct")
 * collections store molecules only and answer queries by brute-force graph
 * matching / on-the-fly fingerprinting over every molecule.
 */
export class Collection {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * The engine revision this module was built from (see [`ENGINE_REV`]).
     */
    static engineRev(): string;
    /**
     * Restore a collection written by [`Collection::to_bytes`]. Rejects (with
     * an error, never a panic) blobs from another engine revision or format
     * version, and any truncated / corrupted blob (checksum + structural
     * checks of every section).
     */
    static fromBytes(bytes: Uint8Array): Collection;
    /**
     * Build an indexed collection from a single batch of SMILES (convenience
     * for callers that don't need progress reporting).
     */
    static fromSmiles(smiles: string[]): Collection;
    /**
     * Whether this collection has a fingerprint index (vs. direct matching).
     */
    isIndexed(): boolean;
    /**
     * The persisted-collection container version (see [`PERSIST_FORMAT_VERSION`]).
     */
    static persistFormatVersion(): number;
    /**
     * Number of scan positions for [`Collection::substructure_search_range`]
     * (== `size`): callers walk `[0, scanLength)` in consecutive ranges.
     */
    scanLength(): number;
    /**
     * Similarity (ECFP4 Tanimoto) search. If `top_k > 0`, returns the top
     * `top_k` hits with Tanimoto >= `threshold` (a floor); otherwise returns
     * every hit with Tanimoto >= `threshold`. Results are sorted by score
     * descending — a {@link SimilarityResult} of parallel `indices` / `scores`.
     */
    similaritySearch(query_smiles: string, threshold: number, top_k: number): SimilarityResult;
    /**
     * Number of molecules in the collection.
     */
    size(): number;
    /**
     * Substructure search with a query molecule (RDKit-js's `get_substruct_match(query)`, as
     * Datagrok's Chem searches with the query its `getQueryMolSafe` builds: a `get_qmol`, or a
     * `get_mol` with `mergeQueryHs`): every molecule `get_substruct_match` would match, as
     * [`Collection::substructure_search`] reports them. The query's atoms and bonds match as
     * RDKit matches them (query features, plain atoms by element, charge, isotope and
     * radicals, H atoms merged into counts as the query holds them); stereo is ignored, as
     * `get_substruct_match` ignores it (spike fast-search-rdkit).
     */
    substructureSearchMol(query: JSMol, limit: number): Uint32Array;
    /**
     * Resumable substructure search over scan positions `[start, end)`
     * (screen-sidecar sorted positions when indexed, input order in direct
     * mode), capped at `limit` hits (0 = unlimited). Walking consecutive
     * ranges from 0 to [`Collection::scan_length`] and concatenating the
     * results (passing the remaining `limit` each time, stopping when it hits
     * 0) yields exactly [`Collection::substructure_search`]'s output — same
     * hits, same order — while letting the host yield between ranges.
     */
    substructureSearchRange(smarts: string, start: number, end: number, limit: number): Uint32Array;
    /**
     * Substructure (SMARTS) search. Returns a `Uint32Array` of the matching
     * molecules' original input positions, capped at `limit` (0 = unlimited).
     */
    substructureSearch(smarts: string, limit: number): Uint32Array;
    /**
     * Serialise an **indexed** collection (CXMOL store, both sidecars, rarity
     * table, original input ids) into one self-validating blob that
     * [`Collection::from_bytes`] restores without re-parsing any SMILES.
     * Layout: `magic[8] | format u32 | flags u32 | rev_len u32 | rev bytes |
     * n u64 | N_SECTIONS × (len u64, bytes) | checksum u64` (all LE; the
     * checksum covers everything before it).
     */
    toBytes(): Uint8Array;
}

/**
 * Incrementally builds a [`Collection`] from SMILES strings.
 *
 * When `build_index` is true the builder also computes the screening + ECFP4
 * fingerprints and assembles the sidecars (fast queries, larger memory). When
 * false it stores molecules only — searches then run by direct graph matching
 * / on-the-fly fingerprinting over every molecule (no build cost, slower
 * queries).
 */
export class CollectionBuilder {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Parse + fingerprint a chunk of SMILES, appending each to the index.
     * Returns the number successfully added; unparseable / unsupported SMILES
     * are skipped and counted in [`CollectionBuilder::failed`]. In a build with
     * RDKit-js's `JSMol` API (chem-wasm, crux-wasm's chem module) an item may be a
     * MOL block: read as Datagrok Chem's `getMolSafe` reads it, its H atoms kept.
     */
    addMany(smiles: string[]): number;
    /**
     * Stream SD-file bytes in (any chunking): each record read so far is added like a
     * SMILES (RDKit's MOL-file reading: sanitize, removeHs). Returns the number added.
     */
    addSdf(chunk: Uint8Array): number;
    /**
     * Number of molecules accepted so far.
     */
    added(): number;
    /**
     * Number of SMILES skipped because they failed to parse / encode.
     */
    failed(): number;
    /**
     * End of the SD input: reads the last record(s). Returns the number added.
     */
    finishSdf(): number;
    /**
     * Finalise into a queryable [`Collection`]. Builds the sidecars + rarity
     * table when indexing; otherwise produces a molecules-only collection.
     */
    finish(): Collection;
    constructor(build_index: boolean);
    /**
     * Read a SMILES RDKit's default read rejects the way RDKit reads it with Kekulize left out
     * of its sanitization (MinimalLib's `get_mol(smiles, {"kekulize": false})`, the retry
     * Datagrok's Chem makes): aromatic rings written without their `[nH]`, aromatic
     * phosphazenes. Off by default.
     */
    setLenient(lenient: boolean): void;
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
    get_atom_pair_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * MinimalLib `get_atom_pair_fp(details)`.
     */
    get_atom_pair_fp(details?: string | null): string;
    get_avalon_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * MinimalLib `get_avalon_fp(details)` (`nBits` 512, `avalonSSSBits`, as there); throws in
     * chem-wasm's one file ([`AVALON_AND_MAPS`]).
     */
    get_avalon_fp(details?: string | null): string;
    /**
     * PORT(SmartsWrite.cpp `MolToCXSmarts`): `get_cxsmarts(details)`: the SMARTS, then
     * `getCXExtensions(mol, CX_ALL)` over the SMARTS writer's atom and bond output order: with a
     * conformer its coordinates and their wedges, and the link nodes, as `get_cxsmiles` writes them.
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
     * `get_inchi(options)`: InChI lives in inchi-wasm (`molToInchi` of `get_molblock()`), the lazy
     * module Chem already ships (decision `chem-wasm`): throws `not supported by Crux`.
     */
    get_inchi(_options?: string | null): string;
    /**
     * [`mol_inchi_input`]: InChI's input for crux-sketch's `inchiFromInput`; `undefined` where RDKit-js's `get_inchi`
     * gives `""` without the library, a throw where it throws.
     */
    get_inchi_input(options?: string | null): Uint8Array | undefined;
    /**
     * PORT(minilib.cpp:158-160): `get_json()`: `MolInterchange::MolToJSONData(mol)`, RDKit's
     * commonchem writer with its extensions (`molfile::mol_interchange`, decision `chem-wasm`).
     */
    get_json(): string;
    /**
     * PORT(minilib.cpp `get_kekule_form`): the MOL block of the Kekulé form.
     */
    get_kekule_form(): string;
    get_maccs_fp_as_uint8array(): Uint8Array;
    /**
     * MinimalLib `get_maccs_fp()`: the 167 MACCS keys as a bit string.
     */
    get_maccs_fp(): string;
    /**
     * PORT(minilib.cpp:592-603, jswrapper.cpp:449-458 `get_mmpa_frags`): `{cores, sidechains}`, two `MolList`s
     * of RDKit's `MMPA::fragmentMol(mol, minCuts, maxCuts, maxCutBonds)` (its default pattern) results, empty when
     * no cut has `minCuts` bonds; both `null` where `fragmentMol` returns false (no cuttable bond, or more than
     * `maxCutBonds`). A single cut's core is empty (`null` from `at` / `next`).
     */
    get_mmpa_frags(min_cuts: number, max_cuts: number, max_cut_bonds: number): any;
    /**
     * `get_molblock(details)`: `MolToMolBlock`; with MinimalLib's details
     * (`includeStereo`, `kekulize`, `useMolBlockWedging`, `addChiralHs`; decision
     * `chem-wasm`: [`crate::minilib::molblock_with_details`]).
     */
    get_molblock(details?: string | null): string;
    get_morgan_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * MinimalLib `get_morgan_fp(details)` (radius 2 by default, as there).
     */
    get_morgan_fp(details?: string | null): string;
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
    get_pattern_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * MinimalLib `get_pattern_fp(details)`.
     */
    get_pattern_fp(details?: string | null): string;
    /**
     * `get_prop(key)`: the molecule property as text (`""` when unset).
     */
    get_prop(key: string): string;
    /**
     * `get_prop_list(includePrivate, includeComputed)`: the molecule property names (an
     * array; private names, `_`-prefixed, only with `includePrivate`).
     */
    get_prop_list(include_private?: boolean | null, _include_computed?: boolean | null): string[];
    get_rdkit_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * MinimalLib `get_rdkit_fp(details)`.
     */
    get_rdkit_fp(details?: string | null): string;
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
    get_topological_torsion_fp_as_uint8array(details?: string | null): Uint8Array;
    /**
     * MinimalLib `get_topological_torsion_fp(details)`.
     */
    get_topological_torsion_fp(details?: string | null): string;
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
     * PORT(minilib.cpp `normalize_depiction`): RDKit's `normalizeDepiction(mol, -1,
     * canonicalize, scaleFactor)` with MinimalLib's canonical transform (its build without
     * Eigen: `depictor::CanonicalTransformBuild::PowerSolver`); the scale factor, or -1
     * (nothing done). Throws where MinimalLib's does (no principal axis).
     */
    normalize_depiction(canonicalize?: number | null, scale_factor?: number | null): number;
    /**
     * PORT(common.h:882-886): `normalize(details)` — `MolStandardize::normalize`.
     */
    normalize(details?: string | null): boolean;
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
 * Embind's `std::vector<JSMolList *>` as `run_reactants` returns it.
 */
export class JSMolListList {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Embind's `delete()`.
     */
    delete(): void;
    /**
     * `get(i)`; `undefined` past the end.
     */
    get(i: number): MolList | undefined;
    /**
     * `size()`.
     */
    size(): number;
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
 * MinimalLib's `RGroupDecomposition` object.
 */
export class RGroupDecomposition {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * PORT(minilib.cpp:997-999): RDKit-js `add(mol)` — the molecule's row, -1 when no core
     * matches, -2 when it matches with no side chains.
     */
    add(mol: JSMol): number;
    /**
     * `add` of a molecule given as SMILES or a MOL block.
     */
    add_text(mol: string): number;
    /**
     * Embind's `delete()` (wasm-bindgen's `free()`).
     */
    delete(): void;
    /**
     * PORT(jswrapper.cpp:489-500): RDKit-js `get_rgroups_as_columns()` — `{Core: MolList, R1:
     * MolList, ...}`, each piece a `JSMol` with RDKit's props (`_rgroupTargetAtoms`, ...).
     */
    get_rgroups_as_columns(): any;
    /**
     * `get_rgroups_as_columns()` as canonical SMILES: `{Core: [smiles, ...], R1: [...], ...}`
     * (`""` where a molecule has no group at a label).
     */
    get_rgroups_as_columns_smiles(): any;
    /**
     * PORT(jswrapper.cpp:502-518): RDKit-js `get_rgroups_as_rows()` — an array of
     * `{Core: JSMol, R1: JSMol, ...}`.
     */
    get_rgroups_as_rows(): any;
    /**
     * `get_rgroups_as_rows()` as canonical SMILES: an array of `{Core: smiles, R1: smiles, ...}`.
     */
    get_rgroups_as_rows_smiles(): any;
    /**
     * PORT(minilib.cpp:1001): `process()`.
     */
    process(): boolean;
}

/**
 * RDKit-js `Reaction` (MinimalLib `JSReaction`).
 */
export class Reaction {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Embind's `delete()`.
     */
    delete(): void;
    /**
     * See [`Self::get_svg`].
     */
    draw_to_canvas(_canvas: any, _width: number, _height: number): string;
    /**
     * See [`Self::get_svg`].
     */
    draw_to_canvas_with_highlights(_canvas: any, _details: string): string;
    /**
     * RDKit `EnumerateLibrary(rxn, reagents, strategy, params)` (not in
     * MinimalLib): an enumerator over `reagentLists` — one array of SMILES per
     * reactant template (unparsable SMILES are skipped). `details` (JSON, all
     * optional): `strategy` (`"cartesian"` default, `"random"`, `"randomAllBBs"`,
     * `"evenPairs"`), `seed` (the random strategies' `minstd_rand`; RDKit's default
     * is 1), `reagentMaxMatchCount`, `sanePartialProducts`, `state` (a
     * `get_state()` text to resume from).
     */
    enumerate(reagent_lists: any, details?: string | null): ReactionEnumerator;
    /**
     * `ChemicalReactionToRxnSmarts` (not in MinimalLib).
     */
    get_smarts(): string;
    /**
     * RDKit-js's reaction drawing (`get_svg`, `get_svg_with_highlights`, `draw_to_canvas*`) is not
     * ported (moldraw2d draws molecules only): these throw `not supported by Crux` rather than
     * being undefined (decision `chem-wasm`).
     */
    get_svg(_width?: number | null, _height?: number | null): string;
    /**
     * See [`Self::get_svg`].
     */
    get_svg_with_highlights(_details: string): string;
    /**
     * `is_valid()`.
     */
    is_valid(): boolean;
    /**
     * PORT(minilib.cpp:644-663): `run_reactants(reactants, maxProducts)` —
     * `initReactantMatchers`, then `runReactants`: one `MolList` of (unsanitized)
     * products per combination of reactant matches. `maxProducts` defaults to
     * 1000 (`JSReaction::maxProducts`).
     */
    run_reactants(reactants: MolList, max_products?: number | null): JSMolListList;
}

/**
 * RDKit `EnumerateLibrary` ([`JSReaction::enumerate`]; `rxn::enumerate`).
 */
export class ReactionEnumerator {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Embind's `delete()`.
     */
    delete(): void;
    /**
     * `getEnumerator().getNumPermutations()` (`2^64 - 1` when it overflows).
     */
    get_num_permutations(): number;
    /**
     * `getEnumerator().getPermutationIdx()`: combinations handed out so far.
     */
    get_permutation_idx(): number;
    /**
     * `getPosition()`: the indices of the last combination into each reagent
     * list the reaction accepts (`get_reagent_counts`).
     */
    get_position(): Float64Array;
    /**
     * `getReagents()` sizes: the reagents kept per reactant template.
     */
    get_reagent_counts(): Float64Array;
    /**
     * `getState()`: RDKit's state text (resume with `set_state` or the
     * `state` detail).
     */
    get_state(): string;
    /**
     * `bool(library)`: are there combinations left (the random strategies never
     * halt).
     */
    has_next(): boolean;
    /**
     * `next()`: the product sets of the next reagent combination (one `MolList`
     * per combination of reactant matches, one product per product template);
     * `undefined` when the enumeration is exhausted.
     */
    next(): JSMolListList | undefined;
    /**
     * `nextSmiles()`: [`ReactionEnumerator::next`] as a JSON array of arrays of
     * SMILES; `undefined` when exhausted.
     */
    next_smiles(): string | undefined;
    /**
     * `resetState()`: back to the start.
     */
    reset_state(): void;
    /**
     * `setState(state)`.
     */
    set_state(state: string): void;
}

/**
 * A similarity-search result as two parallel typed arrays: `indices[i]` is the
 * caller's original input position of the i-th hit, `scores[i]` its Tanimoto
 * similarity in `[0, 1]`. Hits are ordered by descending score.
 */
export class SimilarityResult {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Original input positions of the hits (a `Uint32Array` on the JS side).
     */
    readonly indices: Uint32Array;
    /**
     * Tanimoto scores, parallel to `indices` (a `Float32Array` on the JS side).
     */
    readonly scores: Float32Array;
}

/**
 * `calcDescriptors(molecules, names)`: the JSON of the container's `POST /chem/descriptors` (module docs). A
 * `null` (or any non-string) molecule is a row of nulls, as chem.py's `"M  END" in None` raises.
 */
export function calcDescriptors(molecules: (string | null)[], names: string[]): string;

/**
 * `chemAlerts(molecules, catalogs, atoms)`: the JSON of each catalog's hits per molecule (module docs); throws on a
 * catalog other than Brenk and NIH.
 */
export function chemAlerts(molecules: (string | null)[], catalogs: string[], atoms: boolean): string;

/**
 * `chemCanonicalSmiles(molecules)`: the JSON list of `molecules_to_canonical`'s SMILES (module docs).
 */
export function chemCanonicalSmiles(molecules: (string | null)[]): string;

/**
 * `chemConformers(molecule, options?)`: the JSON of the script's table (module docs).
 */
export function chemConformers(molecule: string, options?: string | null): string;

/**
 * `chemCurate(molecules, options?)`: the JSON list of the script's `curated_molecule` per molecule (module docs).
 */
export function chemCurate(molecules: (string | null)[], options?: string | null): string;

/**
 * `chemEmbed3d(molecule, options)`: the MOL block of `molecule` (a SMILES, or a MOL block when it holds `M  END`)
 * with its RDKit ETKDG coordinates, or `undefined`; throws on bad options and on a timeout (module docs).
 */
export function chemEmbed3d(molecule: string, options?: string | null): string | undefined;

/**
 * `chemEmbed3dRecommendedOptions()`: [`RECOMMENDED_OPTIONS`], the options Chem passes.
 */
export function chemEmbed3dRecommendedOptions(): string;

/**
 * `chemFilterCatalog(molecules, catalog)`: per molecule 1 where the catalog matches, else 0 (module docs).
 */
export function chemFilterCatalog(molecules: (string | null)[], catalog: string): Uint8Array;

/**
 * `chemFingerprints(kind, details, molecules, smiles)`: the fingerprints Chem's RDKit worker code computes with
 * `get_<kind>_fp_as_uint8array(details)` (`kind`: morgan, pattern, rdkit, maccs, atom_pair, topological_torsion),
 * molecule `i`'s at `bytes[offsets[i] .. offsets[i + 1]]` where `present[i]`; `unread` the molecules Chem's getMolSafe
 * reading reads none from; `smiles` the canonical SMILES when asked (chem_batch.rs docs).
 */
export function chemFingerprints(kind: string, details: string | null | undefined, molecules: (string | null)[], smiles: boolean): {bytes: Uint8Array, offsets: Uint32Array, present: Uint8Array, unread: Uint32Array, smiles: (string | null)[] | null};

/**
 * `chemGasteigerMap(molecule, contours, width?, height?)`: the JSON of the panel's picture and the charges, or
 * `null` (module docs).
 */
export function chemGasteigerMap(molecule: string, contours: number, width?: number | null, height?: number | null): string;

/**
 * `chemInchiInputs(molecules)`: per molecule `molecules_to_inchi`'s input for the InChI library, or `null` (module
 * docs).
 */
export function chemInchiInputs(molecules: (string | null)[]): (Uint8Array | null)[];

/**
 * `chemMurckoScaffolds(molecules)`: the JSON list of the script's `scaffolds` per molecule (module docs).
 */
export function chemMurckoScaffolds(molecules: (string | null)[]): string;

/**
 * `chemMutate(molecule, options?)`: the JSON list of the script's `mutations` (module docs).
 */
export function chemMutate(molecule: string, options?: string | null): string;

/**
 * `chemProperties(molecules, names)`: the JSON of Chem's property columns (module docs); throws on a name not in
 * its table.
 */
export function chemProperties(molecules: (string | null)[], names: string[]): string;

/**
 * `chemQed(molecules)`: the JSON of QED and its eight properties per molecule (module docs).
 */
export function chemQed(molecules: (string | null)[]): string;

/**
 * Which parts this build holds: `"core"`, `"heavy"`, `"core,heavy"`, ...; the one file Chem loads (feature `chem`)
 * is `"core,heavy,search,descriptors,scaffold-tree,inchi-input,chem-scripts,chem-conformers,chem-identifiers,chem-murcko,
 * chem-filter-catalogs,synthon-search"`.
 */
export function chemWasmModules(): string;

/**
 * `descriptorTree()`: the JSON of the container's `GET /chem/descriptors/tree`.
 */
export function descriptorTree(): string;

/**
 * Crux extension: `enumerateTautomerSmiles(smiles, params)` — a JSON array of the
 * tautomer SMILES.
 */
export function enumerate_tautomer_smiles(smiles: string, details?: string | null): string;

/**
 * One-shot `RGroupDecompose(cores, mols, asSmiles=True)` as JSON: `{"rows": [...],
 * "unmatched": [mol index, ...], "row_of_mol": [row or null, ...]}`. Cores are
 * SMILES / MOL blocks (or SMARTS with `cores_are_smarts`), molecules SMILES / MOL
 * blocks; `details` is RDKit's RGD parameter JSON.
 */
export function getRGroupDecomposition(cores: string[], mols: string[], details?: string | null, cores_are_smarts?: boolean | null): string;

/**
 * RDKit-js `get_inchikey_for_inchi(inchi)`: in inchi-wasm (`inchiToInchiKey`): throws.
 */
export function get_inchikey_for_inchi(_inchi: string): string;

/**
 * `get_mces(mol1, mol2, details)`: RDKit's `rascalMCES(mol1, mol2, opts)`: the Maximum
 * Common Edge Substructure of two molecules (the largest set of bonds they share, possibly
 * in several fragments), when their similarity estimate reaches `similarityThreshold`.
 *
 * `details`: a JSON object of `RascalOptions` fields by RDKit's names
 * (`similarityThreshold`, `completeAromaticRings`, `ringMatchesRingOnly`,
 * `completeSmallestRings`, `exactConnectionsMatch`, `singleLargestFrag`, `minFragSize`,
 * `maxFragSeparation`, `allBestMCESs`, `maxBestMCESs`, `timeout`, `doEquivBondPruning`,
 * `returnEmptyMCES`, `maxBondMatchPairs`, `equivalentAtoms`, `ignoreBondOrders`,
 * `ignoreAtomAromaticity`, `minCliqueSize`); each defaults to RDKit's except
 * `similarityThreshold`, 0.5 here (RDKit: 0.7). Other keys are ignored. `timeout` (seconds,
 * 60 by default, -1 for none) runs on the page's clock.
 *
 * Returns a JSON array of results, as RDKit's `FindMCES` returns a list: empty when the
 * pair is below the threshold (unless `returnEmptyMCES`), one result by default, every best
 * MCES with `allBestMCESs`. Each result, with the Python wrapper's names:
 * `{"smartsString", "atomMatches": [[atom of mol1, atom of mol2], ...], "bondMatches": [[bond
 * of mol1, bond of mol2], ...], "similarity", "tier1Sim", "tier2Sim", "timedOut",
 * "numFragments", "largestFragmentSize"}` (the indices each molecule's own: what a
 * highlight takes). Throws on details that are not JSON or a bad `equivalentAtoms`.
 */
export function get_mces(mol1: JSMol, mol2: JSMol, details?: string | null): string;

/**
 * RDKit-js `get_mcs_as_json(mols, details)`: `{"smarts","canceled","numAtoms","numBonds"}`,
 * in RDKit-js's key order.
 */
export function get_mcs_as_json(mols: MolList, details?: string | null): string;

/**
 * RDKit-js `get_mcs_as_mol(mols, details)`: the MCS as a query molecule (`get_qmol` of its
 * SMARTS); `null` when there is none.
 */
export function get_mcs_as_mol(mols: MolList, details?: string | null): JSMol | null;

/**
 * RDKit-js `get_mcs_as_smarts(mols, details)`: the MCS's SMARTS.
 */
export function get_mcs_as_smarts(mols: MolList, details?: string | null): string;

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
 * PORT(jswrapper.cpp:464-487): RDKit-js `get_rgd(core, details)`: the decomposition around
 * the core molecule (a `get_mol` / `get_qmol` `JSMol`).
 */
export function get_rgd(core: JSMol, details?: string | null): RGroupDecomposition;

/**
 * `get_rgd` over SMARTS cores (MinimalLib: `get_rgd(get_qmol(smarts), details)`).
 */
export function get_rgd_smarts(cores: string[], details?: string | null): RGroupDecomposition;

/**
 * `get_rgd` of one or more cores given as text (SMILES or MOL block each), `details` =
 * RDKit's RGD parameter JSON (`updateRGroupDecompositionParametersFromJSON`).
 */
export function get_rgd_text(cores: string[], details?: string | null): RGroupDecomposition;

/**
 * PORT(minilib.cpp:856-859): `get_rxn(input, details)` — `undefined` when the
 * input does not read.
 */
export function get_rxn(input: string, details?: string | null): Reaction | undefined;

/**
 * PORT(jswrapper.cpp:156-212, 921-928): RDKit-js `molzip`, its three forms and Crux's one-molecule form (module
 * docs): `molzip(a, b, details?)`, `molzip(row, details?)`, `molzip(mol, details?)`. A `JSMol`, or `null`.
 */
export function molzip(first: JSMol | Record<string, JSMol>, second?: JSMol | string, details?: string): JSMol | null;

/**
 * Python's `Chem.molzipFragments(mols, params)` (RDKit's `molzip(decomposition, params)`, which RDKit-js does not
 * bind): the list's molecules, the core first. A `JSMol`, or `null` for an empty list.
 */
export function molzipFragments(mols: MolList, details?: string | null): JSMol | null;

/**
 * PORT(minilib.cpp:864-868): RDKit-js `prefer_coordgen(prefer)`: `RDDepict::preferCoordGen`,
 * process-wide (CoordGen, crate `coordgen`, decision `coordgen`); `true` throws in a build without
 * CoordGen.
 */
export function prefer_coordgen(prefer: boolean): void;

/**
 * `scaffoldTree(smiles, ringCutoff, discharge)`: the JSON string Datagrok Chem's `scaffold_tree.py` returns for
 * these SMILES (`GenerateScaffoldTree(data, smilesColumn, ringCutoff, dischargeAndDeradicalize)`); throws where
 * the script fails.
 */
export function scaffoldTree(smiles: string[], ring_cutoff: number, discharge: boolean): string;

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
 * Panics as `Error`s with their message in the console (not a bare `unreachable`); `std::sort`'s ties as
 * RDKit-js's (Emscripten's libc++: `smi_compiler::cxx_sort::StdSort::Libcxx`, spike chem-rdkit-gaps-2), the
 * calls that mirror a Python RDKit service excepted ([`PythonDefaults`]).
 */
export function start(): void;

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
    readonly __wbg_chemsmartsset_free: (a: number, b: number) => void;
    readonly __wbg_chemsynthonspace_free: (a: number, b: number) => void;
    readonly __wbg_collection_free: (a: number, b: number) => void;
    readonly __wbg_collectionbuilder_free: (a: number, b: number) => void;
    readonly __wbg_jsmol_free: (a: number, b: number) => void;
    readonly __wbg_jsmollistlist_free: (a: number, b: number) => void;
    readonly __wbg_mollist_free: (a: number, b: number) => void;
    readonly __wbg_reaction_free: (a: number, b: number) => void;
    readonly __wbg_reactionenumerator_free: (a: number, b: number) => void;
    readonly __wbg_rgroupdecomposition_free: (a: number, b: number) => void;
    readonly __wbg_similarityresult_free: (a: number, b: number) => void;
    readonly calcDescriptors: (a: any, b: number, c: number) => [number, number];
    readonly chemAlerts: (a: any, b: number, c: number, d: number) => [number, number, number, number];
    readonly chemCanonicalSmiles: (a: any) => [number, number, number, number];
    readonly chemConformers: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly chemCurate: (a: any, b: number, c: number) => [number, number, number, number];
    readonly chemEmbed3d: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly chemEmbed3dRecommendedOptions: () => [number, number];
    readonly chemFilterCatalog: (a: any, b: number, c: number) => [number, number, number, number];
    readonly chemFingerprints: (a: number, b: number, c: number, d: number, e: any, f: number) => [number, number, number];
    readonly chemGasteigerMap: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly chemInchiInputs: (a: any) => [number, number, number];
    readonly chemMurckoScaffolds: (a: any) => [number, number, number, number];
    readonly chemMutate: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly chemProperties: (a: any, b: number, c: number) => [number, number, number, number];
    readonly chemQed: (a: any) => [number, number];
    readonly chemWasmModules: () => [number, number];
    readonly chemsmartsset_groups: (a: number) => number;
    readonly chemsmartsset_lastWork: (a: number) => [number, number];
    readonly chemsmartsset_match: (a: number, b: any) => [number, number, number];
    readonly chemsmartsset_new: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number];
    readonly chemsmartsset_unreadSmarts: (a: number) => [number, number];
    readonly chemsynthonspace_fingerprintPart: (a: number, b: number, c: number) => [number, number, number, number];
    readonly chemsynthonspace_hasFingerprints: (a: number) => number;
    readonly chemsynthonspace_installFingerprints: (a: number, b: number, c: number) => [number, number];
    readonly chemsynthonspace_new: (a: number, b: number) => [number, number, number];
    readonly chemsynthonspace_numProducts: (a: number) => number;
    readonly chemsynthonspace_numReactions: (a: number) => number;
    readonly chemsynthonspace_numSynthons: (a: number) => number;
    readonly chemsynthonspace_search: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number, number, number];
    readonly collection_engineRev: () => [number, number];
    readonly collection_fromBytes: (a: number, b: number) => [number, number, number];
    readonly collection_fromSmiles: (a: number, b: number) => number;
    readonly collection_isIndexed: (a: number) => number;
    readonly collection_persistFormatVersion: () => number;
    readonly collection_scanLength: (a: number) => number;
    readonly collection_similaritySearch: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly collection_size: (a: number) => number;
    readonly collection_substructureSearch: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly collection_substructureSearchMol: (a: number, b: number, c: number) => [number, number, number, number];
    readonly collection_substructureSearchRange: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly collection_toBytes: (a: number) => [number, number, number, number];
    readonly collectionbuilder_addMany: (a: number, b: number, c: number) => number;
    readonly collectionbuilder_addSdf: (a: number, b: number, c: number) => number;
    readonly collectionbuilder_added: (a: number) => number;
    readonly collectionbuilder_failed: (a: number) => number;
    readonly collectionbuilder_finish: (a: number) => number;
    readonly collectionbuilder_finishSdf: (a: number) => number;
    readonly collectionbuilder_new: (a: number) => number;
    readonly collectionbuilder_setLenient: (a: number, b: number) => void;
    readonly descriptorTree: () => [number, number];
    readonly enumerate_tautomer_smiles: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly getRGroupDecomposition: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number, number, number];
    readonly get_inchikey_for_inchi: (a: number, b: number) => [number, number, number, number];
    readonly get_mces: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly get_mcs_as_json: (a: number, b: number, c: number) => [number, number, number, number];
    readonly get_mcs_as_mol: (a: number, b: number, c: number) => [number, number, number];
    readonly get_mcs_as_smarts: (a: number, b: number, c: number) => [number, number, number, number];
    readonly get_mol: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly get_mol_copy: (a: number) => number;
    readonly get_qmol: (a: number, b: number) => any;
    readonly get_rgd: (a: number, b: number, c: number) => [number, number, number];
    readonly get_rgd_smarts: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly get_rgd_text: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly get_rxn: (a: number, b: number, c: number, d: number) => number;
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
    readonly jsmol_get_avalon_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_avalon_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_cxsmarts: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_cxsmiles: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_descriptors: (a: number) => [number, number, number, number];
    readonly jsmol_get_fp_json: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly jsmol_get_frags: (a: number, b: number, c: number) => [number, number, number];
    readonly jsmol_get_inchi: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_inchi_input: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_json: (a: number) => [number, number, number, number];
    readonly jsmol_get_kekule_form: (a: number) => [number, number, number, number];
    readonly jsmol_get_maccs_fp: (a: number) => [number, number, number, number];
    readonly jsmol_get_maccs_fp_as_uint8array: (a: number) => [number, number, number, number];
    readonly jsmol_get_mmpa_frags: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly jsmol_get_molblock: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_morgan_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_morgan_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_new_coords: (a: number, b: number) => [number, number, number, number];
    readonly jsmol_get_num_atoms: (a: number, b: number) => number;
    readonly jsmol_get_num_bonds: (a: number) => number;
    readonly jsmol_get_pattern_fp: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_pattern_fp_as_uint8array: (a: number, b: number, c: number) => [number, number, number, number];
    readonly jsmol_get_prop: (a: number, b: number, c: number) => [number, number];
    readonly jsmol_get_prop_list: (a: number, b: number, c: number) => [number, number];
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
    readonly jsmollistlist_delete: (a: number) => void;
    readonly jsmollistlist_get: (a: number, b: number) => number;
    readonly jsmollistlist_size: (a: number) => number;
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
    readonly molzip: (a: any, b: number, c: number, d: number) => [number, number, number];
    readonly molzipFragments: (a: number, b: number, c: number) => [number, number, number];
    readonly prefer_coordgen: (a: number) => [number, number];
    readonly reaction_delete: (a: number) => void;
    readonly reaction_draw_to_canvas: (a: number, b: any, c: number, d: number) => [number, number, number, number];
    readonly reaction_draw_to_canvas_with_highlights: (a: number, b: any, c: number, d: number) => [number, number, number, number];
    readonly reaction_enumerate: (a: number, b: any, c: number, d: number) => [number, number, number];
    readonly reaction_get_smarts: (a: number) => [number, number];
    readonly reaction_get_svg: (a: number, b: number, c: number) => [number, number, number, number];
    readonly reaction_get_svg_with_highlights: (a: number, b: number, c: number) => [number, number, number, number];
    readonly reaction_is_valid: (a: number) => number;
    readonly reaction_run_reactants: (a: number, b: number, c: number) => [number, number, number];
    readonly reactionenumerator_delete: (a: number) => void;
    readonly reactionenumerator_get_num_permutations: (a: number) => number;
    readonly reactionenumerator_get_permutation_idx: (a: number) => number;
    readonly reactionenumerator_get_position: (a: number) => [number, number];
    readonly reactionenumerator_get_reagent_counts: (a: number) => [number, number];
    readonly reactionenumerator_get_state: (a: number) => [number, number];
    readonly reactionenumerator_has_next: (a: number) => number;
    readonly reactionenumerator_next: (a: number) => [number, number, number];
    readonly reactionenumerator_next_smiles: (a: number) => [number, number, number, number];
    readonly reactionenumerator_reset_state: (a: number) => void;
    readonly reactionenumerator_set_state: (a: number, b: number, c: number) => [number, number];
    readonly rgroupdecomposition_add: (a: number, b: number) => [number, number, number];
    readonly rgroupdecomposition_add_text: (a: number, b: number, c: number) => [number, number, number];
    readonly rgroupdecomposition_delete: (a: number) => void;
    readonly rgroupdecomposition_get_rgroups_as_columns: (a: number) => [number, number, number];
    readonly rgroupdecomposition_get_rgroups_as_columns_smiles: (a: number) => [number, number, number];
    readonly rgroupdecomposition_get_rgroups_as_rows: (a: number) => [number, number, number];
    readonly rgroupdecomposition_get_rgroups_as_rows_smiles: (a: number) => [number, number, number];
    readonly rgroupdecomposition_process: (a: number) => number;
    readonly scaffoldTree: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly set_log_capture: (a: number, b: number) => [number, number, number];
    readonly set_log_tee: (a: number, b: number) => [number, number, number];
    readonly similarityMapSvg: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly similarityMapWeights: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly similarityresult_indices: (a: number) => [number, number];
    readonly similarityresult_scores: (a: number) => [number, number];
    readonly standardize_smiles: (a: number, b: number) => [number, number, number, number];
    readonly start: () => void;
    readonly use_legacy_stereo_perception: (a: number) => number;
    readonly validate_smiles: (a: number, b: number) => [number, number, number, number];
    readonly version: () => [number, number];
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
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
