/**
 * Synchronous lookup of HELM monomer symbol → natural-analog letter.
 *
 * The central Bio monomer library is fetched once during package init
 * (see `initSequenceTranslatorInt` → `_package.completeInit`) and parked on
 * `_package.bioMonomerLib`. Datagrok's `@init` decorator guarantees this
 * runs before any other package function — including cell-renderer factory
 * calls and the renderer's `render()` invocations — so by the time we look
 * up an analog, the lib is always present.
 *
 * Lookups are memoized in a small Map. No async, no subscriptions, no
 * grid invalidation. If the cache is queried before init has finished
 * (defensive — shouldn't happen in normal flow), we return `null` and the
 * caller falls back to a neutral color; the next render after init will
 * resolve fully.
 */

import {IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {_package} from '../package';

/** symbol → natural-analog letter, or null if absent / no analog. */
const _cache = new Map<string, string | null>();

/** Resolve `symbol` to its single-letter natural analog. Returns `null` for
 * unknowns, the canonical letter (uppercase) for matches. Pure sync. */
export function getNaturalAnalog(symbol: string): string | null {
  if (!symbol) return null;
  const cached = _cache.get(symbol);
  if (cached !== undefined) return cached;

  let lib: IMonomerLib | null = null;
  try { lib = _package.bioMonomerLib; } catch { /* not yet initialized */ }
  if (!lib) return null;

  const analog = lookup(lib, symbol);
  _cache.set(symbol, analog);
  return analog;
}

function lookup(lib: IMonomerLib, symbol: string): string | null {
  for (const pt of lib.getPolymerTypes()) {
    const m = lib.getMonomer(pt, symbol);
    const na = m?.naturalAnalog;
    if (na && typeof na === 'string' && na.length === 1)
      return na.toUpperCase();
  }
  return null;
}
