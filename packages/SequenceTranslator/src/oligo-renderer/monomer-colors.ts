/**
 * Synchronous color resolution against the central Bio monomer library.
 *
 * `_package.bioMonomerLib` is wired up at package init (see
 * `initSequenceTranslatorInt`), so by the time any cell renders, the library
 * is always present. The lib's `getMonomerColors(biotype, symbol)` returns
 * `{ textcolor?, backgroundcolor?, linecolor? }` (all optional) — the
 * library itself handles natural-analog fallback for custom symbols.
 *
 * Lookups are memoized in a small map keyed by `${kind}:${symbol}`.
 */

import {HelmTypes} from '@datagrok-libraries/js-draw-lite/src/types/org';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {_package} from '../package';

/** The four "kinds" we render — maps to the appropriate `HelmType`. */
export type MonomerKind = 'sugar' | 'base' | 'linker' | 'chem';

const KIND_TO_HELM_TYPE: Record<MonomerKind, HelmType> = {
  sugar: HelmTypes.SUGAR as HelmType,
  base: HelmTypes.BASE as HelmType,
  linker: HelmTypes.LINKER as HelmType,
  chem: HelmTypes.CHEM as HelmType,
};

export interface MonomerColorTriple {
  backgroundcolor: string | null;
  textcolor: string | null;
  linecolor: string | null;
}

const EMPTY: MonomerColorTriple = {backgroundcolor: null, textcolor: null, linecolor: null};
const _cache = new Map<string, MonomerColorTriple>();

/** Resolve background/text/line colors for a HELM monomer. Returns `null`s
 * when the library has no entry. Pure sync — backed by `_package.bioMonomerLib`. */
export function getMonomerColors(kind: MonomerKind, symbol: string): MonomerColorTriple {
  if (!symbol) return EMPTY;
  const key = `${kind}:${symbol}`;
  const cached = _cache.get(key);
  if (cached) return cached;

  let result: MonomerColorTriple = EMPTY;
  try {
    const lib = _package.bioMonomerLib;
    const colors = lib.getMonomerColors(KIND_TO_HELM_TYPE[kind], symbol);
    if (colors) {
      result = {
        backgroundcolor: colors.backgroundcolor ?? null,
        textcolor: colors.textcolor ?? null,
        linecolor: colors.linecolor ?? null,
      };
    }
  } catch { /* lib not yet initialized — return all-null */ }

  if (_cache.size > 512) _cache.clear();
  _cache.set(key, result);
  return result;
}
