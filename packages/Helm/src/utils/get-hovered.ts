import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {HelmType, PolymerType, Atom, Mol, ISeqMonomer, HelmMol, HelmAtom} from '@datagrok-libraries/bio/src/helm/types';

import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';

import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {GAP_SYMBOL, GapOriginals, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';

import {HelmMonomerPlacer} from '../helm-monomer-placer';

export function getHoveredMonomerFromEditorMol(
  argsX: number, argsY: number, mol: HelmMol, cellHeight?: number
): HelmAtom | null {
  /** @return {[number, number]} [atom, distance] */
  function getNearest(excluded: (number | undefined)[]): [number | undefined, number | undefined] {
    let atom: number | undefined = undefined;
    let distance: number | undefined = undefined;
    for (let atomI = 0; atomI < mol.atoms.length; ++atomI) {
      if (!excluded.includes(atomI)) {
        const aX: number = mol.atoms[atomI].p.x;
        const aY: number = mol.atoms[atomI].p.y;
        const distanceToAtomI: number = Math.sqrt((argsX - aX) ** 2 + (argsY - aY) ** 2);
        if (distance === undefined || distance > distanceToAtomI) {
          atom = atomI;
          distance = distanceToAtomI;
        }
      }
    }
    return [atom, distance];
  }

  const [firstAtomI, firstDistance] = getNearest([]);
  const [secondAtomI, secondDistance] = getNearest([firstAtomI]);

  let resAtom: HelmAtom | null = null;
  if (firstAtomI !== undefined && firstDistance !== undefined) {
    const firstAtom = mol.atoms[firstAtomI];
    if (secondAtomI !== undefined && secondDistance !== undefined) {
      if (firstDistance < secondDistance * 0.45)
        resAtom = firstAtom;
    } else {
      if (cellHeight && firstDistance < 0.35 * cellHeight)
        resAtom = firstAtom;
    }
  }
  return resAtom;
}

export function getHoveredMonomerFallback(
  argsX: number, _argsY: number, gridCell: DG.GridCell, helmPlacer: HelmMonomerPlacer
): ISeqMonomer | null {
  let hoveredSeqMonomer: ISeqMonomer | null = null;
  const [allParts, lengths, sumLengths] = helmPlacer.getCellAllPartsLengths(gridCell.tableRowIndex!);
  const maxIndex = Object.values(lengths).length - 1;
  let left = 0;
  let right = maxIndex;
  let found = false;
  let iterCount: number = 0;

  let mid = 0;
  if (argsX > sumLengths[0]) {
    while (!found && iterCount < sumLengths.length) {
      mid = Math.floor((right + left) / 2);
      if (argsX >= sumLengths[mid] && argsX <= sumLengths[mid + 1]) {
        left = mid;
        found = true;
      } else if (argsX < sumLengths[mid])
        right = mid - 1;
      else if (argsX > sumLengths[mid + 1])
        left = mid + 1;

      if (left == right)
        found = true;

      iterCount++;
    }
  }
  left = (argsX >= sumLengths[left]) ? left : left - 1; // correct left to between sumLengths
  if (left >= 0)
    hoveredSeqMonomer = getSeqMonomerFromHelm(left, allParts[left], allParts[0]);
  return hoveredSeqMonomer;
}

export function getSeqMonomerFromHelmAtom(atom: HelmAtom): ISeqMonomer {
  const polymerType = helmTypeToPolymerType(atom.bio!.type);
  const canonicalSymbol = atom.elem === GapOriginals[NOTATION.HELM] ? GAP_SYMBOL : atom.elem;
  return {position: parseInt(atom.bio!.continuousId as string) - 1, symbol: atom.elem, biotype: atom.bio!.type};
}

/** Linear
 * @deprecated
 */
function getSeqMonomerFromHelm(pos: number, symbol: string, helmPrefix: string): ISeqMonomer {
  let resSeqMonomer: ISeqMonomer | undefined = undefined;
  const polymerTypeList: PolymerType[] = ['RNA', 'PEPTIDE', 'CHEM', 'BLOB', 'G'];
  for (const polymerType of polymerTypeList) {
    if (helmPrefix.startsWith(polymerType)) {
      let helmType: HelmType = HelmTypes.AA;
      if (polymerType == 'RNA')
        helmType = HelmTypes.NUCLEOTIDE;
      else if (polymerType == 'PEPTIDE')
        helmType = HelmTypes.AA;
      else if (polymerType == 'CHEM')
        helmType = HelmTypes.CHEM;
      resSeqMonomer = {position: pos, symbol: symbol, biotype: helmType};
    }
  }
  if (!resSeqMonomer)
    throw new Error(`Monomer not found for symbol = '${symbol}' and helmPrefix = '${helmPrefix}'.`);
  return resSeqMonomer;
}
