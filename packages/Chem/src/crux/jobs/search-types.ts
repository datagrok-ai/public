/**
 * The search types beside Contains / Not contains (Exact, Stereo agnostic, Included in, Not included in) over a
 * column's index segments in the Crux workers (crux/crux-search-types.ts sends them, crux.worker.ts runs them), and on
 * the main thread's RDKit module for the molecules Crux cannot read. They answer as Chem's RDKit search answers
 * (RdKitService.searchSubstructureWithFps, RdKitServiceWorkerSubstructure.searchSubstructure) on the molecules its
 * workers hold, each row's molecule as getMolSafe reads it (sanitized; the pool caches them as it computes the pattern
 * fingerprints, every row of a column of up to 10,000 rows per worker), by its canonical SMILES (Exact), its SMILES
 * without stereo (Stereo agnostic), and the pattern fingerprint screen and the substructure match of the molecule in
 * the query (Included in).
 */
import type {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {SubstructureSearchType} from '../../constants';
import {getMolSafe, getQueryMolSafe} from '../../utils/mol-creation_rdkit';
import {hasNewLines} from '../../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../../utils/chem-constants';
import {WorkerJob, throwIfTrap} from '../crux-api';

/** The bytes of a pattern fingerprint (RDKit's 2048 bits), as RdKitService.filterByPatternFps compares them. */
const PATTERN_BYTES = 256;
const NO_STEREO = '{"doIsomericSmiles":false}';

/**
 * What a segment's searches read, per row, from the molecule getMolSafe reads (computed when a search first needs it,
 * each molecule parsed once per kind asked): its canonical SMILES (RDKit.js's `get_smiles()`), its SMILES without
 * stereo, its pattern fingerprint (`get_pattern_fp_as_uint8array()`); null where getMolSafe reads none (or the call
 * throws), as the RDKit pool's getFingerprints leaves them.
 */
export class SearchTypeKeys {
  private canonical: (string | null)[] | null = null;
  private flat: (string | null)[] | null = null;
  private patterns: Uint8Array | null = null;
  private hasPattern: Uint8Array | null = null;
  private unreadRows: number[] | null = null;

  constructor(readonly molecules: string[]) {}

  /** The rows with a molecule getMolSafe reads none from (an empty cell is no such row). */
  unread(module: RDModule): number[] {
    this.ensure(module, false);
    return this.unreadRows!;
  }

  canonicalSmiles(module: RDModule): (string | null)[] {
    this.ensure(module, false);
    return this.canonical!;
  }

  flatSmiles(module: RDModule): (string | null)[] {
    this.ensure(module, false);
    return this.flat!;
  }

  /** Each row's pattern fingerprint (PATTERN_BYTES from `row * PATTERN_BYTES`) and whether it has one. */
  patternFps(module: RDModule): {fps: Uint8Array, has: Uint8Array} {
    this.ensure(module, true);
    return {fps: this.patterns!, has: this.hasPattern!};
  }

  private ensure(module: RDModule, patterns: boolean): void {
    const smiles = this.canonical === null;
    if (!smiles && (!patterns || this.patterns !== null))
      return;
    const n = this.molecules.length;
    if (smiles) {
      this.canonical = new Array<string | null>(n).fill(null);
      this.flat = new Array<string | null>(n).fill(null);
      this.unreadRows = [];
    }
    if (patterns) {
      this.patterns = new Uint8Array(n * PATTERN_BYTES);
      this.hasPattern = new Uint8Array(n);
    }
    for (let i = 0; i < n; i++) {
      const molecule = this.molecules[i];
      if (!molecule)
        continue;
      const mol = getMolSafe(molecule, {}, module).mol;
      if (!mol) {
        if (smiles)
          this.unreadRows!.push(i);
        continue;
      }
      try {
        if (smiles) {
          // as getFingerprints writes the canonical SMILES (createSmiles) and the search the one without stereo
          this.canonical![i] = guarded(() => mol.get_smiles());
          this.flat![i] = guarded(() => mol.get_smiles(NO_STEREO));
        }
        if (patterns) {
          const fp = guarded(() => mol.get_pattern_fp_as_uint8array());
          if (fp && fp.length === PATTERN_BYTES) {
            this.patterns!.set(fp, i * PATTERN_BYTES);
            this.hasPattern![i] = 1;
          }
        }
      } finally {
        mol.delete();
      }
    }
  }
}

/** `call()`, or null where it throws (a trap rethrown). */
function guarded<T>(call: () => T): T | null {
  try {
    return call();
  } catch (e) {
    throwIfTrap(e);
    return null;
  }
}

/** The query's SMILES as the RDKit worker's searchSubstructure writes it for Exact (`stereo`) and Stereo agnostic: of
 * `get_mol(query)` (getMolWithSmilesCheck), '' where it reads no molecule. */
function querySmiles(query: string, module: RDModule, stereo: boolean): string {
  if (query && !hasNewLines(query) && query.length > MAX_SMILES_LENGTH)
    return '';
  let mol: RDMol | null = null;
  try {
    mol = module.get_mol(query);
    return mol ? (stereo ? mol.get_smiles() : mol.get_smiles(NO_STEREO)) : '';
  } catch (e) {
    throwIfTrap(e);
    return '';
  } finally {
    mol?.delete();
  }
}

/** Whether a row's pattern fingerprint passes RdKitService.filterByPatternFps' screen: the query's bits in the row's
 * (`superstructure`: the row's bits in the query's, Included in). */
function passesScreen(fps: Uint8Array, row: number, query: Uint8Array, superstructure: boolean): boolean {
  const start = row * PATTERN_BYTES;
  for (let j = 0; j < PATTERN_BYTES; j++) {
    const bits = fps[start + j] & query[j];
    if (bits !== (superstructure ? fps[start + j] : query[j]))
      return false;
  }
  return true;
}

/** Whether the row's molecule (getMolSafe's) matches as the RDKit worker's searchBySearchType asks (Included in and its
 * complement: the molecule is a substructure of the query; Exact, for a query that is no molecule: each a substructure
 * of the other). */
function matchesRow(molecule: string, queryMol: RDMol, module: RDModule, searchType: SubstructureSearchType): boolean {
  const mol = getMolSafe(molecule, {}, module).mol;
  if (!mol)
    return false;
  try {
    if (searchType !== SubstructureSearchType.EXACT_MATCH)
      return queryMol.get_substruct_match(mol) !== '{}';
    return mol.get_substruct_match(queryMol) !== '{}' && queryMol.get_substruct_match(mol) !== '{}';
  } catch (e) {
    throwIfTrap(e);
    return false;
  } finally {
    mol.delete();
  }
}

/**
 * The rows of `keys` that match the query of `searchType` as Chem's RDKit search answers it on getMolSafe's molecules
 * (the file's head): Exact, the rows whose canonical SMILES is the query's (`get_mol(query)`), or for a query that is
 * no molecule (a SMARTS) the rows the pattern screen passes that it matches and that match it; Stereo agnostic, the
 * rows whose SMILES without stereo is the query's (none for a query that is no molecule); Included in, the rows the
 * screen passes (the row's bits in the query's) whose molecule matches in the query getQueryMolSafe builds; Not
 * included in, the other rows with a molecule (RDKit's search scans them all; the screen leaves out no molecule that
 * matches). Throws for a query getQueryMolSafe reads none from.
 */
export function searchTypeHits(module: RDModule, keys: SearchTypeKeys, query: string, molBlockFailover: string,
  searchType: SubstructureSearchType): number[] {
  const queryMol = getQueryMolSafe(query, molBlockFailover, module);
  if (!queryMol)
    throw new Error(`Chem | Search pattern cannot be set: ${query}`);
  try {
    const hits: number[] = [];
    if (searchType === SubstructureSearchType.EXACT_MATCH || searchType === SubstructureSearchType.STEREO_AGNOSTIC) {
      const stereo = searchType === SubstructureSearchType.EXACT_MATCH;
      const smiles = querySmiles(query, module, stereo);
      if (smiles) {
        const rows = stereo ? keys.canonicalSmiles(module) : keys.flatSmiles(module);
        for (let i = 0; i < rows.length; i++) {
          if (rows[i] === smiles)
            hits.push(i);
        }
        return hits;
      }
      if (!stereo) // the RDKit worker matches nothing then (searchBySearchType knows no stereo-agnostic match)
        return hits;
    }
    // RDKit's screen: Exact, the query's bits in the row's (a SMARTS query's bits are a query's: the reverse screen,
    // the row's in the query's, would leave out rows it matches); Included in, the row's in the query's
    const superstructure = searchType !== SubstructureSearchType.EXACT_MATCH;
    const queryFp = queryMol.get_pattern_fp_as_uint8array();
    const {fps, has} = keys.patternFps(module);
    for (let i = 0; i < keys.molecules.length; i++) {
      if (has[i] && passesScreen(fps, i, queryFp, superstructure) &&
        matchesRow(keys.molecules[i], queryMol, module, searchType))
        hits.push(i);
    }
    if (searchType !== SubstructureSearchType.NOT_INCLUDED_IN)
      return hits;
    const canonical = keys.canonicalSmiles(module);
    const included = new Set(hits);
    return canonical.flatMap((smiles, i) => smiles !== null && !included.has(i) ? [i] : []);
  } finally {
    queryMol.delete();
  }
}

/** The segments' keys a Crux worker holds, by segment key (CruxService's index segments). */
const segmentKeys = new Map<string, SearchTypeKeys>();

/** The answer of a segment: its hits and the rows Crux reads no molecule from (RDKit decides for them), both
 * segment-local; null when the worker holds no keys of the segment (the caller sends its molecules again). */
export type SegmentSearchResult = {hits: Uint32Array, unread: Uint32Array} | null;

/** The jobs by name (the Crux worker's `op`). */
export const SEARCH_TYPE_JOBS: {[op: string]: WorkerJob} = {
  /** A search over one segment; `molecules` (its rows' molecules) when the worker is to (re)build its keys. */
  'searchTypes::search': ({module}, key: string, molecules: string[] | null, query: string,
    molBlockFailover: string, searchType: SubstructureSearchType): SegmentSearchResult => {
    if (molecules)
      segmentKeys.set(key, new SearchTypeKeys(molecules));
    const keys = segmentKeys.get(key);
    if (!keys)
      return null;
    const hits = searchTypeHits(module, keys, query, molBlockFailover, searchType);
    return {hits: new Uint32Array(hits), unread: new Uint32Array(keys.unread(module))};
  },
};

/** Frees the keys of dropped segments (crux.worker.ts `drop`). */
export function dropSearchTypeKeys(keys: string[]): void {
  for (const key of keys)
    segmentKeys.delete(key);
}
