/* tslint:disable */
/* eslint-disable */

/**
 * A row of what InChI cannot hold (crux-sketch N6), sketch-wasm's `LossCode` where the meaning is the same. The first six
 * stop the InChI ('' is written): `r-label`, `dummy` (dummy atoms and pseudoatoms), `query` (query atoms and bonds, link
 * nodes, position-variation and multicentre bonds), `dative`, `hydrogen-bond`, `zero-order`. The rest are written with
 * the warning: `stereo-groups` (OR and AND groups, written as absolute), `stereo` (atropisomer axes, not written),
 * `attachment-points`, `alias` (aliases on real atoms), `sgroup`.
 */
export type InchiLossCode = 'r-label' | 'dummy' | 'query' | 'dative' | 'hydrogen-bond' | 'zero-order'
| 'stereo-groups' | 'stereo' | 'attachment-points' | 'alias' | 'sgroup';
/**
 * A warning of the InChI library: `undefined-stereo` ("Omitted undefined stereo"), `charges-rearranged` ("Charges were
 * rearranged"), `unusual-valence` ("Accepted unusual valence(s): C(2)"), `ambiguous-stereo` ("Ambiguous stereo:
 * center(s)"), and `library` for any other, its message the library's text. The library's normalisations (metal and salt
 * disconnection, protons added or removed, charges neutralized) are no warning.
 */
export type InchiLibraryCode = 'undefined-stereo' | 'charges-rearranged' | 'unusual-valence' | 'ambiguous-stereo' | 'library';
/** What InChI cannot hold of the molecule written, in sketch-wasm's `writeWarnings` shape: one row a code. */
export interface WriteWarning {
    /** 'inchi'. */
    format: string;
    code: InchiLossCode | InchiLibraryCode;
    message: string;
    /** The atoms concerned (indices of the MOL block); empty for a library warning. */
    atoms: number[];
    /** The bonds concerned (indices of the MOL block); empty for a library warning. */
    bonds: number[];
}
/**
 * Why a conversion failed: `empty` no text; `syntax` a MOL block that does not read, an InChI the library makes no
 * structure of; `unsupported` an InChIKey, a MOL block feature the reader does not read; `invalid` a structure RDKit
 * cannot sanitize, or one `MolToInchi` cannot kekulize.
 */
export type InchiErrorKind = 'empty' | 'syntax' | 'unsupported' | 'invalid';
/** What the module throws (a plain object), in sketch-wasm's `ParseError` shape, with its kind. */
export interface ParseError {
    /** 'inchi', 'molV2000' or 'molV3000'. */
    format: string;
    kind: InchiErrorKind;
    message: string;
    /** The 1-based line of a MOL block error. */
    line?: number;
}
/** `molToInchi`'s result: RDKit's `MolToInchi` and its `ExtraInchiReturnValues`, with the write warnings. */
export interface InchiWritten {
    /** The InChI; '' for an empty molecule, one InChI cannot hold, and one the library fails on (code 2 or 3). */
    inchi: string;
    /** The AuxInfo ('' with no InChI). */
    auxInfo: string;
    /** The library's return code: 0 okay, 1 warning, 2 error, 3 fatal; null when it was not asked (empty, or refused). */
    code: number | null;
    /** The library's messages, `; ` between them. */
    message: string;
    /** The library's log. */
    log: string;
    /** What InChI cannot hold, then the library's warnings. */
    writeWarnings: WriteWarning[];
}
/** `inchiToMol`'s result. */
export interface MolRead {
    /** RDKit's canonical SMILES of the molecule `MolFromInchi` makes. */
    smiles: string;
    /** The library's return code: 0 okay, 1 warning. */
    code: number;
    message: string;
    log: string;
}



/**
 * An InChI's InChIKey, RDKit.js's `get_inchikey_for_inchi` ([`convert::inchi_to_inchi_key`]).
 */
export function inchiToInchiKey(inchi: string): string;

/**
 * An InChI's molecule, RDKit's `MolFromInchi`, as its canonical SMILES ([`convert::inchi_to_mol`]).
 */
export function inchiToMol(inchi: string): MolRead;

/**
 * A MOL block's InChI, as RDKit.js's `get_inchi(options)` writes it, with RDKit's
 * `ExtraInchiReturnValues` and the write warnings ([`convert::mol_to_inchi`]).
 */
export function molToInchi(molblock: string, options?: string | null): InchiWritten;

/**
 * A MOL block's InChIKey, RDKit's `MolToInchiKey` ([`convert::mol_to_inchi_key`]).
 */
export function molToInchiKey(molblock: string, options?: string | null): string;

/**
 * The crate version.
 */
export function version(): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly inchiToInchiKey: (a: number, b: number) => [number, number];
    readonly inchiToMol: (a: number, b: number) => [number, number, number];
    readonly molToInchi: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly molToInchiKey: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly version: () => [number, number];
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
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
