/**
 * Shared types, constants, and cache utilities for the interactive atom
 * highlighting bridge between Chem (2D SMILES) and Molstar (3D poses).
 *
 * Extracted from molstar-viewer.ts to keep the viewer class focused on
 * Molstar lifecycle management while this module handles the cross-package
 * event protocol and atom-mapping types.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Structure} from 'molstar/lib/mol-model/structure';

import {_package} from '../../package';

// -- Constants ---------------------------------------------------------------

/** Event name for cross-package atom selection synchronization.
 *  Mirrors CHEM_INTERACTIVE_SELECTION_EVENT from Chem's constants.ts. */
export const CHEM_SELECTION_EVENT = 'chem-interactive-selection-changed';

/** Event name fired BY this viewer when the user hovers a ligand atom in
 *  the 3D pose. Chem's rdkit-cell-renderer listens and renders the
 *  corresponding 2D atom highlight (the reverse of CHEM_SELECTION_EVENT).
 *  Mirrors CHEM_MOL3D_HOVER_EVENT from Chem's constants.ts. */
export const CHEM_MOL3D_HOVER_EVENT = 'chem-mol3d-hover-changed';

/** Payload shape for CHEM_MOL3D_HOVER_EVENT. `atom3DSerial` is the Molstar
 *  1-based atom id (same numbering that `computeSerials` produces for the
 *  forward direction), or null when the cursor leaves all ligand atoms.
 *  `mode` mirrors the 2D modifier semantics. */
export interface Mol3DHoverEventArgs {
  mol3DColumnName: string;
  rowIdx: number;
  atom3DSerial: number | null;
  mode: 'preview' | 'paint' | 'erase';
}

// -- Types -------------------------------------------------------------------

/** Atom-index mapping produced by Chem's mapAtomIndices2Dto3D.
 *  Mirrors the AtomIndexMapping interface from atom-index-mapper.ts
 *  (defined here to avoid a cross-package import). */
export interface AtomMapping3D {
  mapping: number[];
  method: string;
  mappedCount: number;
  pdbSerials?: number[];
}

/** Shape of the cross-package `chem-interactive-selection-changed` event. */
export interface ChemSelectionEventArgs {
  column?: unknown;
  rowIdx: number;
  atoms: number[];
  mapping3D?: AtomMapping3D | null;
  persistent?: boolean;
  clearAll?: boolean;
  mol3DColumnName?: string;
}

/** Cached entry for atom selection events. */
export interface SelectionCacheEntry {
  atoms: number[];
  mapping3D: AtomMapping3D | null;
}

// -- Cache -------------------------------------------------------------------

/** Builds a composite cache key: dfId-dfName-columnName-rowIdx. */
export function selectionCacheKey(
  dfId: string, dfName: string, colName: string, rowIdx: number,
): string {
  return `${dfId}-${dfName}-${colName}-${rowIdx}`;
}

/**
 * Per-row cache of atom selection events. Uses composite keys
 * (dfId-dfName-colName-rowIdx) to avoid collisions across different
 * dataframes, columns, or tabs.
 */
let _selectionCache: DG.LruCache<string, SelectionCacheEntry> | null = null;

export function getSelectionCache(): DG.LruCache<string, SelectionCacheEntry> {
  if (!_selectionCache)
    _selectionCache = new DG.LruCache<string, SelectionCacheEntry>();
  return _selectionCache;
}

// Register at module load — no lazy guard needed.
grok.events.onCustomEvent(CHEM_SELECTION_EVENT)
  .subscribe((_args: unknown) => {
    const {
      rowIdx = -1, atoms = [], persistent, clearAll, mapping3D, column,
    } = (_args as ChemSelectionEventArgs) ?? {};
    const isPersistent = persistent !== false;
    const col = column as DG.Column | undefined;
    const dfId = col?.dataFrame?.id ?? '';
    const dfName = col?.dataFrame?.name ?? '';
    const colName = col?.name ?? '';

    _package.logger.debug(
      `[molstar-picker-global] caching selection event atomsLen=${atoms.length} rowIdx=${rowIdx} persistent=${isPersistent}`);

    if (isPersistent) {
      const cache = getSelectionCache();
      if (clearAll)
        _selectionCache = new DG.LruCache<string, SelectionCacheEntry>();
      else {
        const key = selectionCacheKey(dfId, dfName, colName, rowIdx);
        if (atoms.length > 0)
          cache.set(key, {atoms, mapping3D: mapping3D ?? null});
        else
          cache.set(key, {atoms: [], mapping3D: null});
      }
    }
  });

// -- Pure computations -------------------------------------------------------

/** Converts 2D atom indices to 3D PDB serial numbers using the
 *  pre-computed mapping. Falls back to heavy-atom serial order from
 *  the Molstar Structure, then to naive index+1. */
export function computeSerials(
  atomIndices: number[],
  mapping3D?: AtomMapping3D | null,
  structure?: Structure,
): number[] {
  if (mapping3D?.mapping) {
    const pdbSerials = mapping3D.pdbSerials;
    const serials: number[] = [];
    for (const i of atomIndices) {
      const mapped = mapping3D.mapping[i];
      if (mapped < 0) continue;
      serials.push(
        pdbSerials && mapped < pdbSerials.length ? pdbSerials[mapped] : mapped + 1);
    }
    _package.logger.debug(
      `[molstar-picker] _computeSerials: method=${mapping3D.method} atoms=[${atomIndices}] serials=[${serials}] hasPdbSerials=${!!pdbSerials}`);
    return serials;
  }
  if (structure) {
    const heavySerials: number[] = [];
    try {
      for (const unit of structure.units) {
        const {elements} = unit;
        const atomicNumber = unit.model.atomicHierarchy.atoms.type_symbol;
        for (let j = 0; j < elements.length; j++) {
          const eI = elements[j];
          const sym = atomicNumber.value(eI);
          if (sym !== 'H' && sym !== 'D')
            heavySerials.push(eI + 1);
        }
      }
    } catch {
      /* fall through */
    }
    if (heavySerials.length > 0) {
      const serials: number[] = [];
      for (const i of atomIndices) {
        if (i >= 0 && i < heavySerials.length)
          serials.push(heavySerials[i]);
      }
      return serials;
    }
  }
  return atomIndices.map((i) => i + 1);
}
