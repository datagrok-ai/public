import * as DG from 'datagrok-api/dg';

import {MonomerMapValue} from './types';
import {IMonomerLibBase} from '../types/monomer-library';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';


/** The monomers' atoms and bonds in their library colours; `alpha` overrides the colours' own opacity (the one monomer
 * under the pointer stands out more than a whole molecule coloured by its monomers). */
export function getMolHighlight(monomerMaps: Iterable<MonomerMapValue>, monomerLib: IMonomerLibBase, alpha?: number): ISubstruct {
  const hlAtoms: { [key: number]: number[] } = {};
  const hlBonds: { [key: number]: number[] } = {};

  for (const monomerMapValue of monomerMaps) {
    const wem = monomerLib.getWebEditorMonomer(monomerMapValue.biotype, monomerMapValue.symbol)!;
    const mColorStr = wem.backgroundcolor ?? DG.Color.toRgb(DG.Color.mouseOverRows);
    // a library colour is any CSS colour: the natural amino acids' are rgb(...), read as hex they were all red
    const mColor = DG.Color.fromHtml(mColorStr);
    const rgba = DG.Color.hexToPercentRgb(mColorStr) ??
      (mColor != null ? [DG.Color.r(mColor) / 255, DG.Color.g(mColor) / 255, DG.Color.b(mColor) / 255, 0.3] : [1.0, 0.0, 0.0, 0.7]);
    const mColorA = alpha === undefined ? rgba : [rgba[0], rgba[1], rgba[2], alpha];
    for (const mAtom of monomerMapValue.atoms)
      hlAtoms[mAtom] = mColorA;
    for (const mBond of monomerMapValue.bonds)
      hlBonds[mBond] = mColorA;
  }

  const resSubstruct: ISubstruct = {
    atoms: Object.keys(hlAtoms).map((k) => parseInt(k)),
    bonds: Object.keys(hlBonds).map((k) => parseInt(k)),
    highlightAtomColors: hlAtoms,
    highlightBondColors: hlBonds,
  };

  return resSubstruct;
}
