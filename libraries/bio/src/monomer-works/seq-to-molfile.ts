import * as DG from 'datagrok-api/dg';

import {MonomerMapValue} from './types';
import {hexToPercentRgb} from './utils';
import {IMonomerLibBase} from '../types/monomer-library';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';


export function getMolHighlight(monomerMaps: Iterable<MonomerMapValue>, monomerLib: IMonomerLibBase): ISubstruct {
  const hlAtoms: { [key: number]: number[] } = {};
  const hlBonds: { [key: number]: number[] } = {};

  for (const monomerMapValue of monomerMaps) {
    const wem = monomerLib.getWebEditorMonomer(monomerMapValue.biotype, monomerMapValue.symbol)!;
    const mColorStr = wem.backgroundcolor;
    // eslint-disable-next-line deprecation/deprecation
    const mColorA = hexToPercentRgb(mColorStr ?? DG.Color.toRgb(DG.Color.mouseOverRows)) ?? [1.0, 0.0, 0.0, 0.7];
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
