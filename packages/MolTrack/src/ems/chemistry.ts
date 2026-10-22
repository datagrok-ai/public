// Crux (vendor/crux) reproduces RDKit 2026.03.6: the legacy backend's standardization and RegistrationHash.
import initCrux, * as crux from '../../vendor/crux/crux_wasm.js';
import initInchi, * as inchi from '../../vendor/crux/inchi_wasm.js';
import * as DG from 'datagrok-api/dg';

export const MATCH_RULES = ['ALL_LAYERS', 'STEREO_INSENSITIVE_LAYERS', 'TAUTOMER_INSENSITIVE_LAYERS'] as const;

export type MatchRule = typeof MATCH_RULES[number];

export type HashColumn = 'hash_all_layers' | 'hash_stereo_insensitive' | 'hash_tautomer_insensitive';

export const HASH_COLUMNS: Record<MatchRule, HashColumn> = {
  ALL_LAYERS: 'hash_all_layers',
  STEREO_INSENSITIVE_LAYERS: 'hash_stereo_insensitive',
  TAUTOMER_INSENSITIVE_LAYERS: 'hash_tautomer_insensitive',
};

export interface CompoundIdentity {
  structure: string;
  hash_mol: string;
  match_rule: MatchRule;
  hash_all_layers: string;
  hash_stereo_insensitive: string;
  hash_tautomer_insensitive: string;
  inchi: string;
  inchikey: string;
  formula: string;
}

export const IDENTITY_COLUMNS: (keyof CompoundIdentity)[] = ['structure', 'hash_mol', 'match_rule', 'hash_all_layers',
  'hash_stereo_insensitive', 'hash_tautomer_insensitive', 'inchi', 'inchikey', 'formula'];

type Mol = crux.JSMol;

export class Chemistry {
  static readonly FINGERPRINT = '{"radius":2,"nBits":2048}';
  private static _loading: Promise<void> | null = null;
  private static _moltrack: DG.Package | undefined;

  // the package-test bundle's own `_package` is never initialized: this is the published package either way
  static package(): DG.Package {
    return Chemistry._moltrack ??= DG.Func.find({package: 'Moltrack', name: 'molTrackEmsApp'})[0].package;
  }

  static load(): Promise<void> {
    Chemistry._loading ??= (async () => {
      const dist = `${Chemistry.package().webRoot.replace(/\/?$/, '/')}dist/`;
      await Promise.all([initCrux({module_or_path: `${dist}crux_wasm_bg.wasm`}),
        initInchi({module_or_path: `${dist}inchi_wasm_bg.wasm`})]);
    })().catch((e) => {
      Chemistry._loading = null;
      throw e;
    });
    return Chemistry._loading;
  }

  static async identify(input: string, rule: MatchRule): Promise<CompoundIdentity> {
    await Chemistry.load();
    return Chemistry._identify(input, rule);
  }

  static async substructure(structures: string[], query: string): Promise<number[]> {
    await Chemistry.load();
    const q = crux.get_qmol(query);
    if (q === null)
      throw new Error(`Cannot read the query structure: ${query}`);
    try {
      return Chemistry._indices(structures, (mol) => mol.get_substruct_match(q) !== '{}');
    } finally {
      q.delete();
    }
  }

  /** "is substructure of": the indices of `structures` that `query` contains */
  static async superstructure(structures: string[], query: string): Promise<number[]> {
    await Chemistry.load();
    const target = Chemistry._read(query);
    try {
      return Chemistry._indices(structures, (mol) => target.get_substruct_match(mol) !== '{}');
    } finally {
      target.delete();
    }
  }

  static async similarity(structures: string[], query: string): Promise<Float64Array> {
    await Chemistry.load();
    const target = Chemistry._read(query);
    let fp: Uint8Array;
    try {
      fp = target.get_morgan_fp_as_uint8array(Chemistry.FINGERPRINT);
    } finally {
      target.delete();
    }
    const scores = new Float64Array(structures.length).fill(NaN);
    Chemistry._each(structures, (mol, i) =>
      scores[i] = Chemistry.tanimoto(fp, mol.get_morgan_fp_as_uint8array(Chemistry.FINGERPRINT)));
    return scores;
  }

  static tanimoto(a: Uint8Array, b: Uint8Array): number {
    let both = 0;
    let either = 0;
    for (let i = 0; i < a.length; i++) {
      both += Chemistry._bits(a[i] & b[i]);
      either += Chemistry._bits(a[i] | b[i]);
    }
    return either === 0 ? 0 : both / either;
  }

  private static _identify(input: string, rule: MatchRule): CompoundIdentity {
    const mol = Chemistry._read(input);
    const std = crux.get_mol_copy(mol);
    try {
      std.cleanup();
      std.fragment_parent();
      std.uncharge();
      const smiles = std.get_smiles();
      const hashes = Object.fromEntries(MATCH_RULES.map((r) => [r,
        crux.registrationHash(smiles, r, true) as {hash: string, layers: Record<string, string>}]));
      const written = inchi.molToInchi(mol.get_molblock());
      return {
        structure: hashes[rule].layers.CANONICAL_SMILES, hash_mol: hashes[rule].hash, match_rule: rule,
        hash_all_layers: hashes.ALL_LAYERS.hash, hash_stereo_insensitive: hashes.STEREO_INSENSITIVE_LAYERS.hash,
        hash_tautomer_insensitive: hashes.TAUTOMER_INSENSITIVE_LAYERS.hash,
        inchi: written.inchi, inchikey: written.inchi ? inchi.inchiToInchiKey(written.inchi) : '',
        formula: crux.molHash(mol.get_smiles(), 'MolFormula', false, 0),
      };
    } finally {
      std.delete();
      mol.delete();
    }
  }

  private static _bits(x: number): number {
    let n = 0;
    for (; x; x &= x - 1)
      n++;
    return n;
  }

  private static _parse(input: string): Mol | null {
    return input ? crux.get_mol(input) : null;
  }

  private static _read(input: string): Mol {
    const mol = Chemistry._parse(input);
    if (mol === null)
      throw new Error('Cannot read the structure');
    return mol;
  }

  private static _each(structures: string[], visit: (mol: Mol, i: number) => void): void {
    for (const [i, s] of structures.entries()) {
      const mol = Chemistry._parse(s);
      if (mol === null)
        continue;
      try {
        visit(mol, i);
      } finally {
        mol.delete();
      }
    }
  }

  private static _indices(structures: string[], matches: (mol: Mol) => boolean): number[] {
    const hits: number[] = [];
    Chemistry._each(structures, (mol, i) => {
      if (matches(mol))
        hits.push(i);
    });
    return hits;
  }
}
