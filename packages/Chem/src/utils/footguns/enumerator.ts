/**
 * The Reaction Enumerator's filters, checked before a run (catalogue TRF-05, TRF-06, TRF-07): the enumerator drops a
 * blocking SMARTS (a template row's blocking_fg) or an exclusion SMARTS it does not read without a word (enumerate.ts
 * parseTemplate, filters.ts buildExclusionQmols), so a typo turns a filter off. In Crux mode a SMARTS no engine reads
 * is a block (Enumerate disabled, the words beside it); a blocking SMARTS that matches none of the building blocks is a
 * warn (it blocks nothing). An exclusion SMARTS that matched no product is said after the run (enumerate.ts). Pure: the
 * modules come in, so the package tests run it in either mode.
 */
import type {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {flag} from './rules';
import type {Flag} from './rules';

/** What the checks read of a run's inputs (shared.ts buildInputs gives them). */
export interface EnumeratorInputs {
  templates: {blockingSmartsList: string[]}[];
  buildingBlocks: string[];
  exclusionSmarts: string[];
}

/** A query of `smarts` read by `module`, or null (deleted by the caller). */
function query(module: RDModule, smarts: string): RDMol | null {
  let q: RDMol | null = null;
  try {
    q = module.get_qmol(smarts);
    if (q && q.is_valid())
      return q;
  } catch (e) {
    onError(e); // unread
  }
  q?.delete();
  return null;
}

/** What hears an error of an engine call (checks.ts: crux-module.ts noteCruxError, which turns Crux off on a trap). */
let onError: (e: unknown) => void = () => {};

/** Sets what hears the engine calls' errors. */
export function setEnumeratorErrorHook(hook: (e: unknown) => void): void {
  onError = hook;
}

/** Whether either engine reads `smarts` as a query (the run falls back to RDKit where Crux does not read it). */
function readsQuery(crux: RDModule, rdkit: RDModule | null, smarts: string): boolean {
  for (const module of [crux, rdkit]) {
    const q = module ? query(module, smarts) : null;
    if (q) {
      q.delete();
      return true;
    }
  }
  return false;
}

/** Whether `q` matches `mol` (a match that throws is none). */
function matches(mol: RDMol, q: RDMol): boolean {
  try {
    const raw = mol.get_substruct_match(q);
    return !!raw && raw !== '{}' && raw !== '[]';
  } catch (e) {
    onError(e);
    return false;
  }
}

/** Whether a blocking SMARTS matches any building block, by SMARTS and building blocks (the enumerator's validation
 * runs on every input change, a building block file can hold thousands). */
const matchCache = new Map<string, boolean>();
const MATCH_CACHE_SIZE = 200;

function matchesAnyBuildingBlock(crux: RDModule, smarts: string, buildingBlocks: string[]): boolean {
  const key = `${smarts}\u0001${buildingBlocks.length}\u0001${buildingBlocks.join('.')}`;
  const known = matchCache.get(key);
  if (known !== undefined)
    return known;
  const q = query(crux, smarts);
  let found = q === null; // read by RDKit alone: not judged here
  if (q) {
    try {
      for (const bb of buildingBlocks) {
        let mol: RDMol | null = null;
        try {
          mol = crux.get_mol(bb);
          if (mol && matches(mol, q)) {
            found = true;
            break;
          }
        } catch (e) {
          onError(e); // an unreadable building block matches nothing
        } finally {
          mol?.delete();
        }
      }
    } finally {
      q.delete();
    }
  }
  if (matchCache.size >= MATCH_CACHE_SIZE)
    matchCache.delete(matchCache.keys().next().value!);
  matchCache.set(key, found);
  return found;
}

/**
 * A template's blocking SMARTS as written: the enumerator cuts a blocking cell at ';' and '|' (shared.ts
 * extractTemplates), so a SMARTS holding ';' (its low-precedence AND: `[N;H2]`) comes in pieces whose brackets do not
 * balance. Such pieces are joined again with ';' (`cut`: the enumerator reads them apart); a piece that balances stays
 * as it is (without the join, "[N;H2" was named as "[N").
 */
function rejoined(pieces: string[]): {text: string, cut: boolean}[] {
  const out: {text: string, cut: boolean}[] = [];
  let held: string[] = [];
  let depth = 0;
  for (const p of pieces) {
    held.push(p);
    depth += (p.match(/\[/g)?.length ?? 0) - (p.match(/\]/g)?.length ?? 0);
    if (depth <= 0) {
      out.push({text: held.join(';'), cut: held.length > 1});
      held = [];
      depth = 0;
    }
  }
  if (held.length)
    out.push({text: held.join(';'), cut: held.length > 1});
  return out;
}

/**
 * The flags of an enumeration's filters: a block for each blocking or exclusion SMARTS neither engine reads (naming it
 * and its template), a warn for each blocking SMARTS that matches none of the building blocks. `crux` is Crux's module
 * (the checks run in Crux mode only), `rdkit` RDKit's (a SMARTS it reads is run by it).
 */
export function enumeratorFlags(crux: RDModule, rdkit: RDModule | null, inputs: EnumeratorInputs): Flag[] {
  const flags: Flag[] = [];
  inputs.templates.forEach((t, i) => {
    for (const {text: s, cut} of rejoined(t.blockingSmartsList)) {
      if (cut && readsQuery(crux, rdkit, s))
        flags.push(flag('blockingCut', `${s} (template ${i + 1})`));
      else if (!readsQuery(crux, rdkit, s))
        flags.push(flag('blockingUnreadable', `${s} (template ${i + 1})`));
      else if (inputs.buildingBlocks.length && !matchesAnyBuildingBlock(crux, s, inputs.buildingBlocks))
        flags.push(flag('blockingMatchesNothing', `${s} (template ${i + 1})`));
    }
  });
  for (const s of inputs.exclusionSmarts) {
    if (!readsQuery(crux, rdkit, s))
      flags.push(flag('exclusionUnreadable', s));
  }
  return flags;
}
