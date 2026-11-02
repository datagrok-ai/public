import * as DG from 'datagrok-api/dg';

import {PolymerType} from '../helm/types';
import {PolymerTypes} from '../helm/consts';
import {IMonomerLib, IMonomerLibBase, Monomer} from '../types/monomer-library';
import {MmcrTemps} from '../utils/cell-renderer-consts';

/* The effective monomer library for the macromolecule assembler (Bio on Crux, PLAN D8): the monomers Datagrok
 * resolves symbols with (merged libraries, duplicate preferences, an override library) as a Pistoia JSON array,
 * under a key that changes whenever the library may have. The assembler compiles a library once per worker per key. */

export type AsmLibrary = { key: string, json: string };

const POLYMER_TYPES: PolymerType[] = [PolymerTypes.PEPTIDE, PolymerTypes.RNA, PolymerTypes.CHEM, PolymerTypes.BLOB];

/** The library a sequence column resolves with: its override library (PolyTool, cell renderers) or the system one. */
export function getSeqColMonomerLib(seqCol: DG.Column<string>, sysMonomerLib: IMonomerLibBase): IMonomerLibBase {
  return seqCol.temp[MmcrTemps.overriddenLibrary] ?? sysMonomerLib;
}

/** A library monomer, not a placeholder of addMissingMonomer or a cached inline-SMILES monomer (#P1, #R2, ...). */
export function isLibraryMonomer(m: Monomer | null): m is Monomer {
  if (!m) return false;
  if (m.lib) return true;
  return m.author !== 'MISSING' && !/^#[A-Z]\d+$/.test(m.symbol) && !!(m.molfile || m.smiles);
}

function getPolymerTypes(lib: IMonomerLibBase): PolymerType[] {
  const libTypes = (lib as Partial<IMonomerLib>).getPolymerTypes?.() ?? [];
  return [...POLYMER_TYPES, ...libTypes.filter((pt) => !POLYMER_TYPES.includes(pt))];
}

/** The effective monomers: every monomer the library's lookup returns for a listed symbol, keyed by its own
 * polymer type and symbol, without runtime fields. Symbols the lookup remaps (RNA R/P to r/p, cached SMILES) and
 * placeholders are skipped; for an override library, its monomers shadow the base's.
 * @param {IMonomerLibBase} lib
 * @return {Monomer[]} */
export function getEffectiveMonomers(lib: IMonomerLibBase): Monomer[] {
  const res: Monomer[] = [];
  const seen = new Set<string>();
  for (const pt of getPolymerTypes(lib)) {
    for (const symbol of lib.getMonomerSymbolsByType(pt)) {
      const m = lib.getMonomer(pt, symbol);
      if (!isLibraryMonomer(m) || m.symbol !== symbol) continue;
      const key = `${m.polymerType ?? pt}\u0001${m.symbol}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const out: Monomer = {...m, lib: undefined, wem: undefined};
      if (!out.polymerType) out.polymerType = pt;
      res.push(out);
    }
  }
  return res;
}

/** Library version counters, per onChanged observable: an override library shares its base's observable. */
const versionByObservable = new WeakMap<object, { version: number }>();
const sourcesByLib = new WeakMap<IMonomerLibBase, { version: number, sources: string }>();

function getLibVersion(lib: IMonomerLibBase): number {
  const changed = lib.onChanged;
  if (!changed) return 0;
  let state = versionByObservable.get(changed);
  if (!state) {
    const newState = {version: 0};
    changed.subscribe(() => {++newState.version;});
    versionByObservable.set(changed, newState);
    state = newState;
  }
  return state.version;
}

/** The sources (library file names) of the effective monomers, computed once per library version. */
function getLibSources(lib: IMonomerLibBase, version: number): string {
  const cached = sourcesByLib.get(lib);
  if (cached && cached.version === version) return cached.sources;
  const sources = new Set<string>();
  for (const pt of getPolymerTypes(lib)) {
    for (const symbol of lib.getMonomerSymbolsByType(pt)) {
      const src = lib.getMonomer(pt, symbol)?.lib?.source;
      if (src) sources.add(src);
    }
  }
  const res = [...sources].sort().join(',');
  sourcesByLib.set(lib, {version, sources: res});
  return res;
}

/** The winners of duplicated symbols: user preferences change them without an onChanged event. */
function getDuplicatesFingerprint(lib: IMonomerLibBase): string {
  const duplicates = (lib as Partial<IMonomerLib>).duplicateMonomers;
  if (!duplicates) return '';
  const parts: string[] = [];
  for (const [pt, bySymbol] of Object.entries(duplicates)) {
    for (const symbol of Object.keys(bySymbol))
      parts.push(`${pt}:${symbol}=${lib.getMonomer(pt as PolymerType, symbol)?.lib?.source ?? ''}`);
  }
  return parts.sort().join(',');
}

/** FNV-1a, 32 bits, as 8 hex digits. */
function hashString(s: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < s.length; ++i) {
    h ^= s.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, '0');
}

/** The library's cache key: its source (an override library's is unique per override), the version counter bumped
 * on onChanged, and a hash of the effective sources and duplicate winners.
 * @param {IMonomerLibBase} lib
 * @return {string} */
export function getMonomerLibKey(lib: IMonomerLibBase): string {
  const version = getLibVersion(lib);
  const sources = getLibSources(lib, version);
  const duplicates = getDuplicatesFingerprint(lib);
  return `${lib.source}#${version}:${hashString(`${sources}|${duplicates}`)}`;
}

const JSON_CACHE_MAX = 4;
const jsonByKey = new Map<string, string>();

/** The effective library as a Pistoia JSON array with its key; the JSON is built once per key (last 4 kept).
 * @param {IMonomerLibBase} lib the library to resolve with, an override library included
 * @return {AsmLibrary} */
export function getMonomerLibJson(lib: IMonomerLibBase): AsmLibrary {
  const key = getMonomerLibKey(lib);
  let json = jsonByKey.get(key);
  if (json === undefined) {
    json = JSON.stringify(getEffectiveMonomers(lib));
    if (jsonByKey.size >= JSON_CACHE_MAX) jsonByKey.delete(jsonByKey.keys().next().value!);
  } else
    jsonByKey.delete(key); // re-inserted below as the most recent
  jsonByKey.set(key, json);
  return {key, json};
}
