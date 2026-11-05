import {ISeqMonomer, HelmMol, HelmAtom} from '@datagrok-libraries/bio/src/helm/types';

import {SeqValueBase} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

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

export function getSeqMonomerFromHelmAtom(seqValue: SeqValueBase, atom: HelmAtom): ISeqMonomer {
  const pos: number = atom.bio!.continuousId - 1;
  let canonicalSymbol = atom.elem;
  try { // for DNA/RNA, we might add phosphates and sugars
    if (!seqValue.isDna() && !seqValue.isRna())
      canonicalSymbol = seqValue.getSplitted().getCanonical(pos);
    else
      canonicalSymbol = seqValue.getSplittedWithSugarsAndPhosphates().getCanonical(pos);
  } catch (_) {}
  return {position: atom.bio!.continuousId - 1, symbol: canonicalSymbol, biotype: atom.bio!.type};
}
