import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {Atom, IHelmBio, HelmMol, HelmType, JSDraw2ModuleType, OrgType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Linkage} from './pt-misc';

declare const JSDraw2: JSDraw2ModuleType;
declare const org: OrgType;


export type PtBio = IHelmBio & { i: number, j: number };
type PtAtom = Atom<HelmType, PtBio>

export function getHelmMol(linkages: Linkage[], mainFragments: string[][], helmHelper: IHelmHelper): HelmMol {
  const resHwe = helmHelper.createHelmWebEditor();
  const resMol = resHwe.editor.m;

  let counter = 0;
  const p = new JSDraw2.Point(0, 0);
  for (let i = 0; i < mainFragments.length; i++) {
    for (let j = 0; j < mainFragments[i].length; j++) {
      if (!!mainFragments[i][j]) {
        const elem = mainFragments[i][j];
        const bio: PtBio = {type: HelmTypes.AA, i: i, j: j, continuousId: counter};
        const atom = new JSDraw2.Atom<HelmType, IHelmBio>(p, elem, bio);
        resMol.addAtom(atom);

        if (j !== 0) {
          const atom1 = resMol.atoms[counter - 1];
          const atom2 = resMol.atoms[counter];
          const bond = new JSDraw2.Bond<HelmType, IHelmBio>(atom1, atom2);
          bond.r1 = 2;
          bond.r2 = 1;
          resMol.addBond(bond);
        }

        counter++;
        p.x += JSDraw2.Editor.BONDLENGTH; // Inspired by HELMWebEditor
      }
      p.y += 4 * JSDraw2.Editor.BONDLENGTH; // Inspired by HELMWebEditor
    }
  }

  for (let i = 0; i < linkages.length; i++) {
    const atom1 = resMol.atoms[linkages[i].fMonomer - 1];
    const atom2 = resMol.atoms[linkages[i].sMonomer - 1];
    const bond = new JSDraw2.Bond<HelmType, IHelmBio>(atom1, atom2);
    bond.r1 = linkages[i].fR;
    bond.r2 = linkages[i].sR;
    resMol.addBond(bond);
  }
  return resMol;
}

export function helmMolToNotation(mol: HelmMol): string {
  const atoms = mol.atoms;
  const bonds = mol.bonds;
  const chains: number[] = [];
  const specialBonds: number[] = [];
  for (let i = 0; i < bonds.length!; i++) {
    //@ts-ignore
    if (bonds[i].a1.bio.i !== bonds[i].a2.bio.i)
      specialBonds.push(i);
  }

  for (let i = 0; i < atoms.length!; i++) {
    //@ts-ignore
    const atomChain = atoms[i].bio?.i;
    if (atomChain + 1 > chains.length)
      chains.push(1);
    else
      chains[atomChain]++;
  }

  const simpleChains: string[][] = new Array(chains.length);
  let counter = 0;
  for (let i = 0; i < chains.length!; i++) {
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
        //@ts-ignore
        const firstAtomLinks = bonds[specialBonds[0]].a1.bio.i == 0 && bonds[specialBonds[0]].a1.bio.j == 0;
        //@ts-ignore
        const secondAtomLinks = bonds[specialBonds[0]].a2.bio.i == 1 && bonds[specialBonds[0]].a1.bio.j == 0;
        if (firstAtomLinks && secondAtomLinks)
          chainAdd += '-';
      }
    }

    res += chainAdd;
  }

  return res;
}
