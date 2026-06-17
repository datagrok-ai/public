import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Linkage} from './pt-misc';

// hwe migration (Phase 6): PolyTool no longer builds a live JSDraw2 / HELM Web
// Editor `HelmMol` (which the standalone `@datagrok-libraries/hwe` model — being
// immutable — cannot reproduce). The "mol" assembled here is a private graph
// consumed only by `helmMolToNotation` (template round-trip) and `Chain.check`
// (count/elem/continuousId validation); both read a tiny field subset. So we
// model it with plain objects instead of the editor's mutable mol graph,
// dropping the `JSDraw2`/`org.helm.webeditor` globals entirely.

export type PtBio = {type: HelmType, i: number, j: number, continuousId: number};
export type PtAtom = {T: 'ATOM', elem: string, bio: PtBio};
export type PtBond = {a1: PtAtom, a2: PtAtom, r1: number, r2: number};
export type PtMol = {atoms: PtAtom[], bonds: PtBond[], name?: string};

export function getHelmMol(linkages: Linkage[], mainFragments: string[][], _helmHelper?: IHelmHelper): PtMol {
  const atoms: PtAtom[] = [];
  const bonds: PtBond[] = [];

  let counter = 0;
  for (let i = 0; i < mainFragments.length; i++) {
    for (let j = 0; j < mainFragments[i].length; j++) {
      if (!!mainFragments[i][j]) {
        const elem = mainFragments[i][j];
        const bio: PtBio = {type: HelmTypes.AA, i: i, j: j, continuousId: counter};
        const atom: PtAtom = {T: 'ATOM', elem, bio};
        atoms.push(atom);

        if (j !== 0) {
          const atom1 = atoms[counter - 1];
          const atom2 = atoms[counter];
          bonds.push({a1: atom1, a2: atom2, r1: 2, r2: 1});
        }

        counter++;
      }
    }
  }

  for (let i = 0; i < linkages.length; i++) {
    const atom1 = atoms[linkages[i].fMonomer - 1];
    const atom2 = atoms[linkages[i].sMonomer - 1];
    bonds.push({a1: atom1, a2: atom2, r1: linkages[i].fR, r2: linkages[i].sR});
  }
  return {atoms, bonds};
}

export function helmMolToNotation(mol: PtMol): string {
  const atoms = mol.atoms;
  const bonds = mol.bonds;
  const chains: number[] = [];
  const specialBonds: number[] = [];
  for (let i = 0; i < bonds.length; i++) {
    if (bonds[i].a1.bio.i !== bonds[i].a2.bio.i)
      specialBonds.push(i);
  }

  for (let i = 0; i < atoms.length; i++) {
    const atomChain = atoms[i].bio.i;
    if (atomChain + 1 > chains.length)
      chains.push(1);
    else
      chains[atomChain]++;
  }

  const simpleChains: string[][] = new Array(chains.length);
  let counter = 0;
  for (let i = 0; i < chains.length; i++) {
    const simpleChain: string[] = new Array(chains[i]);
    for (let j = 0; j < chains[i]; j++) {
      simpleChain[j] = atoms[counter].elem;
      counter++;
    }

    simpleChains[i] = simpleChain;
  }

  let res = '';
  for (let i = 0; i < simpleChains.length; i++) {
    let chainAdd = '';

    for (let j = 0; j < simpleChains[i].length; j++)
      chainAdd += `${j == 0 ? '' : '-'}${simpleChains[i][j]}`;

    if (i !== 0) {
      const rxp = /(\(.\d+\))/;
      const match = chainAdd.match(rxp);
      chainAdd = chainAdd.replace(match?.[0]!, '');
      const group = match ? match?.[0]! : '';
      chainAdd = `${group}{${chainAdd}}`;
    } else {
      if (simpleChains.length > 1) {
        const firstAtomLinks = bonds[specialBonds[0]].a1.bio.i == 0 && bonds[specialBonds[0]].a1.bio.j == 0;
        const secondAtomLinks = bonds[specialBonds[0]].a2.bio.i == 1 && bonds[specialBonds[0]].a1.bio.j == 0;
        if (firstAtomLinks && secondAtomLinks)
          chainAdd += '-';
      }
    }

    res += chainAdd;
  }

  return res;
}
