/* eslint-disable max-len */
/**
 * Types, constants, and pure helpers for the 2D↔3D atom-highlighting bridge.
 * Extracted from molstar-viewer.ts; the selection cache lives on MolstarViewer.selectionCache.
 */

import {Structure} from 'molstar/lib/mol-model/structure';
import {AtomMapping3D} from '@datagrok-libraries/chem-meta/src/types';

import {_package} from '../../package';

// -- Constants ---------------------------------------------------------------

/** Payload for CHEM_MOL3D_HOVER_EVENT. `atom3DSerial: null` = cursor left the atom. */

/** Cached entry for atom selection events. */
export interface SelectionCacheEntry {
  atoms: number[];
  mapping3D: AtomMapping3D | null;
}

// -- Pure computations -------------------------------------------------------

export function selectionCacheKey(
  dfId: string, dfName: string, colName: string, rowIdx: number,
): string {
  return `${dfId}-${dfName}-${colName}-${rowIdx}`;
}

/** Converts 2D atom indices to 3D PDB serial numbers. Falls back to
 *  heavy-atom serial order from the Molstar Structure, then index+1. */
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
