import {Atom, IHelmBio, HelmMol, HelmType, JSDraw2ModuleType, OrgType,} from '@datagrok-libraries/bio/src/helm/types';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {cleanupHelmSymbol} from '@datagrok-libraries/bio/src/helm/utils';

import {getOuterIdx, getInnerIdx} from './pt-misc';
import {Rules, RuleLink, RuleReaction, getMonomerPairs} from './pt-rules';

declare const JSDraw2: JSDraw2ModuleType;
declare const org: OrgType;

type Linkage = {
  fChain: number,
  sChain: number,
  /** Continuous 1-based numbering */ fMonomer: number,
  /** Continuous 1-based numbering */ sMonomer: number,
  fR: number,
  sR: number
}
export type PtBio = IHelmBio & { i: number, j: number };

type PtAtom = Atom<HelmType, PtBio>

export class Chain {
  linkages: Linkage[];
  monomers: string[][];
  mol: HelmMol;

  constructor(
    monomers: string[][], linkages: Linkage[], mol: HelmMol,
    protected helmHelper: IHelmHelper,
  ) {
    this.linkages = linkages;
    this.monomers = monomers;
    this.mol = mol;
  }

  /** Parse harmonized sequence (template) from pseudo helm */
  static parseHelm(helm: string, helmHelper: IHelmHelper) {
    const ea = /(\w+\{.*\})\$(.*)\$(.*)\$(.*)\$/g.exec(helm)!;
    // const fragmentation = helm.split('$');
    const fragmentation = [ea[1], ea[2], ea[3], ea[4]];

    const rawFragments = fragmentation[0].split('|');
    const rawLinkages = fragmentation[1].split('|');

    const monomers = new Array<Array<string>>(rawFragments.length);
    const linkages: Linkage[] = [];

    const resHwe = helmHelper.createHelmWebEditor();
    const resMol = resHwe.editor.m;

    let counter = 0;
    const p = new JSDraw2.Point(0, 0);
    //HELM parsing
    for (let i = 0; i < rawFragments.length; i++) {
      const idxStart = rawFragments[i].indexOf('{');
      const idxEnd = rawFragments[i].indexOf('}');

      monomers[i] = rawFragments[i].slice(idxStart + 1, idxEnd).split('.').map((s) => cleanupHelmSymbol(s));
      for (let j = 0; j < monomers[i].length; j++) {
        const elem = monomers[i][j];
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

    //HELM parsing
    for (let i = 0; i < rawLinkages.length; i++) {
      if (rawLinkages[i] !== '' && rawLinkages[i] !== 'V2.0') {
        const rawData = rawLinkages[i].split(',');
        const fChainIdx = parseInt(rawData[0].replace('PEPTIDE', '')) - 1;
        const sChainIdx = parseInt(rawData[1].replace('PEPTIDE', '')) - 1;
        const rawDataConnections = rawData[2].split('-');
        const rawDataConnection1 = rawDataConnections[0].split(':');
        const rawDataConnection2 = rawDataConnections[1].split(':');

        linkages.push({
          fChain: fChainIdx,
          sChain: sChainIdx,
          fMonomer: getOuterIdx(parseInt(rawDataConnection1[0]), fChainIdx, monomers),
          sMonomer: getOuterIdx(parseInt(rawDataConnection2[0]), sChainIdx, monomers),
          fR: parseInt(rawDataConnection1[1].replace('R', '')),
          sR: parseInt(rawDataConnection2[1].replace('R', '')),
        });
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

    return new Chain(monomers, linkages, resMol, helmHelper);
  }

  /** Get macromolecule from harmonized sequence (template) */
  applyRules(rules: Rules): Chain {
    const newNotation = this.getNotation();
    return Chain.fromNotation(newNotation, rules, this.helmHelper);
  }

  static fromNotation(sequence: string, rules: Rules, helmHelper: IHelmHelper): Chain {
    const heterodimerCode = rules.heterodimerCode;
    const homodimerCode = rules.homodimerCode;

    const mainFragments: string[] = [];
    const linkages: Linkage[] = [];

    //NOTICE: this works only with simple single heterodimers
    const heterodimeric = heterodimerCode !== null ? sequence.split(`(${rules.heterodimerCode!})`) : '';
    if (heterodimerCode !== null && heterodimeric.length > 1) {
      linkages.push({fChain: 0, sChain: 1, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
      mainFragments.push(heterodimeric[1].replaceAll('{', '').replaceAll('}', ''));
      mainFragments.push(heterodimeric[2].replaceAll('{', '').replaceAll('}', ''));
    } else {
      mainFragments.push(sequence);
    }

    //NOTICE: this works only with simple single dimers
    for (let i = 0; i < mainFragments.length; i++) {
      if (homodimerCode !== null && mainFragments[i].includes(`(${homodimerCode!})`)) {
        const idxSequence = mainFragments.length;

        linkages.push({fChain: i, sChain: idxSequence, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
        const rawDimer = mainFragments[i].replace(`(${homodimerCode!})`, '');
        const idx = rawDimer.indexOf('{');
        const linker = rawDimer.slice(0, idx);
        const body = rawDimer.replace(linker, '').replaceAll('{', '').replaceAll('}', '');

        mainFragments[i] = linker + body;
        mainFragments.push(body);
      }
    }

    for (let i = 0; i < mainFragments.length; i++) {
      if (homodimerCode !== null && mainFragments[i].includes(`(${homodimerCode!})`)) {
        const idxSequence = mainFragments.length;

        linkages.push({fChain: i, sChain: idxSequence, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
        const rawDimer = mainFragments[i].replace(`(${homodimerCode!})`, '');
        const idx = rawDimer.indexOf('{');
        const linker = rawDimer.slice(0, idx);
        const body = rawDimer.replace(linker, '').replaceAll('{', '').replaceAll('}', '');

        mainFragments[i] = linker + body;
        mainFragments.push(body);
      }
    }

    const monomers = new Array<Array<string>>(mainFragments.length);

    for (let i = 0; i < mainFragments.length; i++) {
      const rawMonomers = mainFragments[i].split('-');
      const linkedPositions = this.getLinkedPositions(rawMonomers, rules.linkRules);
      const [monomersCycled, allPos1, allPos2, allAttaches1, allAttaches2] =
        this.getAllCycles(rules.linkRules, rawMonomers, linkedPositions);

      const monomersReady = new Array<string>(monomersCycled.length);
      // for (let j = 0; j < monomersCycled.length; j++)
      //   monomersReady[j] = `[${monomersCycled[j]}]`;

      for (let j = 0; j < allPos1.length; j++) {
        linkages.push({
          fChain: i,
          sChain: i,
          fMonomer: allPos1[j],
          sMonomer: allPos2[j],
          fR: allAttaches1[j],
          sR: allAttaches2[j],
        });
      }

      monomers[i] = monomersCycled;
    }

    const monomersAll: string[][] = [];

    const resHwe = helmHelper.createHelmWebEditor();
    const resMol = resHwe.editor.m;

    let counter = 0;
    const p = new JSDraw2.Point(0, 0);
    for (let i = 0; i < monomers.length; i++) {
      const linkedPositions = this.getLinkedPositions(monomers[i], rules.reactionRules);
      const [monomersCycled, allPos1, allPos2, ruleN] =
        this.getAllReactants(rules.reactionRules, monomers[i], linkedPositions);

      if (allPos1.length >= 1) {
        const ch1 = new Array<string>(allPos2[0] - 1);
        const ch2 = new Array<string>(monomersCycled.length - allPos2[0]);
        for (let j = 0; j < allPos2[0] - 1; j++) {
          const elem = ch1[j] = monomersCycled[j];
          const bio: PtBio = {type: HelmTypes.AA, i: i, j: j, continuousId: counter};
          const atom: PtAtom = new JSDraw2.Atom<HelmType, PtBio>(p, elem, bio);
          resMol.addAtom(atom);

          if (j > 0) {
            const atom1 = resMol.atoms[counter - 1];
            const atom2 = resMol.atoms[counter];
            const bond = new JSDraw2.Bond<HelmType, IHelmBio>(atom1, atom2);
            bond.r1 = 2;
            bond.r2 = 1;
            resMol.addBond(bond);
          }
          counter++;
        }

        for (let j = allPos2[0]; j < monomersCycled.length; j++) {
          const elem = ch2[j - allPos2[0]] = monomersCycled[j];
          const bio: PtBio = {type: HelmTypes.AA, i: i, j: j, continuousId: counter};
          const atom: PtAtom = new JSDraw2.Atom<HelmType, PtBio>(p, elem, bio);
          resMol.addAtom(atom);

          if (j > allPos2[0]) {
            const atom1 = resMol.atoms[counter - 1];
            const atom2 = resMol.atoms[counter];
            const bond = new JSDraw2.Bond<HelmType, IHelmBio>(atom1, atom2);
            bond.r1 = 2;
            bond.r2 = 1;
            resMol.addBond(bond);
          }
          counter++;
        }

        resMol.atoms[allPos1[0] - 1].elem = ch1[allPos1[0] - 1] = rules.reactionRules[ruleN[0]].name;

        for (let j = 0; j < linkages.length; j++) {
          if (linkages[j].fMonomer > allPos2[0]) {
            linkages[j].fMonomer -= allPos2[0];
            linkages[j].fChain++;
          }
          if (linkages[j].sMonomer > allPos2[0]) {
            linkages[j].sMonomer -= allPos2[0];
            linkages[j].sChain++;
          }
        }
        linkages.push({
          fChain: 0,
          sChain: 0,
          fMonomer: allPos1[0],
          sMonomer: allPos2[0] - 1,
          fR: 3,
          sR: 2,
        });

        linkages.push({
          fChain: 0,
          sChain: 1,
          fMonomer: allPos1[0],
          sMonomer: 1,
          fR: 4,
          sR: 1,
        });

        const monomersReady1 = new Array<string>(ch1.length);
        for (let j = 0; j < ch1.length; j++)
          monomersReady1[j] = `[${ch1[j]}]`;
        const monomersReady2 = new Array<string>(ch2.length);
        for (let j = 0; j < ch2.length; j++)
          monomersReady2[j] = `[${ch2[j]}]`;

        monomersAll.push(ch1);
        monomersAll.push(ch2);
      } else {
        for (let j = 0; j < monomers[i].length; j++) {
          const elem = monomers[i][j];
          const bio: PtBio = {type: HelmTypes.AA, i: i, j: j, continuousId: counter};
          const atom: PtAtom = new JSDraw2.Atom<HelmType, PtBio>(p, elem, bio);
          resMol.addAtom(atom);

          if (j > 0) {
            const atom1 = resMol.atoms[counter - 1];
            const atom2 = resMol.atoms[counter];
            const bond = new JSDraw2.Bond<HelmType, IHelmBio>(atom1, atom2);
            bond.r1 = 2;
            bond.r2 = 1;
            resMol.addBond(bond);
          }
          counter++;
        }
        monomersAll.push(monomers[i]);
      }
    }

    for (const l of linkages) {
      const atom1 = resMol.atoms[l.fMonomer - 1];
      const atom2 = resMol.atoms[l.sMonomer - 1];
      const bond = new JSDraw2.Bond<HelmType, IHelmBio>(atom1, atom2);
      bond.r1 = l.fR;
      bond.r2 = l.sR;
      resMol.addBond(bond);
    }

    const chain = new Chain(monomersAll, linkages, resMol, helmHelper);
    return chain;
  }

  /** Parse harmonized sequence notation (template)  */
  static parseNotation(sequence: string, helmHelper: IHelmHelper): Chain {
    const mainFragments: string[][] = [];

    const linkages: Linkage[] = [];

    const resHwe = helmHelper.createHelmWebEditor();
    const resMol = resHwe.editor.m;

    const rxp = /(\(.\d+\))?\{[^\}]*\}/g;
    const seqs: string [] = [];
    seqs.push(sequence.replaceAll(rxp, ''));

    //const l = (rxpRes?.length) ?? -1;

    const matches = sequence.matchAll(rxp);
    //const rxpRes = rxp.exec(sequence);
    for (const m of matches) {
      const str = m![0];
      if (str)
        seqs.push(str);
    }

    let counter = 0;
    for (let i = 0; i < seqs.length; i++) {
      const splMonomers = seqs[i].split('-');
      const monomers: string [] = new Array<string>(splMonomers.length);
      let spmCount: number = 0;
      for (let j = 0; j < splMonomers.length; j++) {
        const monomer = splMonomers[j].replace('{', '').replace('}', '');
        if (monomer !== '') {
          monomers[j] = monomer;
          counter++;
          spmCount++;
        } else {
          linkages.push({fChain: i, sChain: i + 1, fMonomer: counter, sMonomer: counter + 1, fR: 1, sR: 1});
        }
      }
      mainFragments.push(monomers.slice(0, spmCount));
    }

    counter = 0;
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

    const chain = new Chain(mainFragments, linkages, resMol, helmHelper);
    return chain;
  }

  getHelmChanged(changeNumber: number, monomer: string): string {
    //TODO: make more efficient
    let counter = 0;
    let idx1 = 0;
    let idx2 = 0;
    loop1:
    for (let i = 0; i < this.monomers.length; i++) {
      loop2:
      for (let j = 0; j < this.monomers[i].length; j++) {
        if (counter == changeNumber) {
          idx1 = i;
          idx2 = j;
          break loop1;
        }
        counter++;
      }
    }

    const previous = this.monomers[idx1][idx2];

    this.monomers[idx1][idx2] = `[${monomer}]`;
    const res = this.getHelm();
    this.monomers[idx1][idx2] = previous;

    return res;
  }

  /** Gets harmonized sequence (template) pseudo helm */
  getNotationHelm(): string {
    return this.getHelm();
  }

  /** Gets harmonized sequence (template) pseudo helm */
  getHelm(): string {
    let helm = '';
    for (let i = 0; i < this.monomers.length; i++) {
      if (i > 0)
        helm += '|';

      helm += `PEPTIDE${i + 1}{`;

      for (let j = 0; j < this.monomers[i].length; j++) {
        if (j > 0)
          helm += '.';
        const symbol = this.monomers[i][j];
        helm += symbol.length > 1 ? `[${symbol}]` : symbol;
      }
      helm += `}`;
    }

    helm += '$';

    for (let i = 0; i < this.linkages.length; i++) {
      if (i > 0)
        helm += '|';
      helm += `PEPTIDE${this.linkages[i].fChain + 1},PEPTIDE${this.linkages[i].sChain + 1},`;

      helm += `${getInnerIdx(this.linkages[i].fMonomer - 1, this.monomers)[0] + 1}:R${this.linkages[i].fR}-`;
      helm += `${getInnerIdx(this.linkages[i].sMonomer - 1, this.monomers)[0] + 1}:R${this.linkages[i].sR}`;
    }

    helm += '$$$' + 'V2.0';
    return helm;
  }

  getNotation(): string {
    const atoms = this.mol.atoms;
    const bonds = this.mol.bonds;
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

  protected static getLinkedPositions(monomers: string[], rules: RuleLink[] | RuleReaction []):
    [number, number, number][] {
    const result: [number, number, number][] = [];

    for (let i = 0; i < rules.length; i++) {
      let firstFound = false;
      let secondFound = false;
      let firstIsFirst = false;
      let firstEntryIndex = -1;
      let secondEntryIndex = -1;
      const add = `(${rules[i].code})`;

      const [firstMonomers, secondMonomers] = getMonomerPairs(rules[i]);

      if (firstMonomers.length > 0) {
        for (let j = 0; j < firstMonomers.length; j++) {
          for (let k = 0; k < monomers.length; k++) {
            if (monomers[k].includes(add)) {
              if (firstFound) {
                if (firstIsFirst && monomers[k] == secondMonomers[j] + add) {
                  secondFound = true;
                  secondEntryIndex = k;
                  break;
                } else if (!firstIsFirst && monomers[k] == firstMonomers[j] + add) {
                  secondFound = true;
                  secondEntryIndex = k;
                  break;
                } else {
                  continue;
                }
              } else {
                if (monomers[k] == firstMonomers[j] + add) {
                  firstFound = true;
                  firstIsFirst = true;
                  firstEntryIndex = k;
                } else if (monomers[k] == secondMonomers[j] + add) {
                  firstFound = true;
                  firstIsFirst = false;
                  firstEntryIndex = k;
                } else {
                  continue;
                }
              }
            }
          }

          if (!(firstFound && secondFound))
            result.push([-1, -1, -1]);
          else if (firstIsFirst)
            result.push([firstEntryIndex, secondEntryIndex, i]);
          else
            result.push([secondEntryIndex, firstEntryIndex, i]);
        }
      } else {
        for (let k = 0; k < monomers.length; k++) {
          if (monomers[k].includes(add)) {
            if (firstFound) {
              if (firstIsFirst && monomers[k]) {
                secondFound = true;
                secondEntryIndex = k;
                break;
              } else if (!firstIsFirst && monomers[k]) {
                secondFound = true;
                secondEntryIndex = k;
                break;
              } else {
                continue;
              }
            } else {
              firstFound = true;
              firstIsFirst = true;
              firstEntryIndex = k;
            }
          }
        }

        if (!(firstFound && secondFound))
          result.push([-1, -1, -1]);
        else if (firstIsFirst)
          result.push([firstEntryIndex, secondEntryIndex, i]);
        else
          result.push([secondEntryIndex, firstEntryIndex, i]);
      }
    }

    return result;
  }

  protected static getAllCycles(rules: RuleLink[], monomers: string [], positions: [number, number, number][]):
    [string [], number [], number [], number [], number []] {
    const allPos1: number [] = [];
    const allPos2: number [] = [];
    const allAttaches1: number [] = [];
    const allAttaches2: number [] = [];
    const ruleCount = rules.length;

    for (let i = 0; i < ruleCount; i++) {
      if (positions[i][0] == -1)
        continue;

      const firstMonomer = monomers[positions[i][0]];
      const secondMonomer = monomers[positions[i][1]];

      // monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstSubstitution);
      // monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondSubstitution);
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(0)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(0)', '');
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(1)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(1)', '');
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(2)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(2)', '');
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(3)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(3)', '');

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      allAttaches1.push(rules[i].firstLinkingGroup);
      allAttaches2.push(rules[i].secondLinkingGroup);
    }

    return [monomers, allPos1, allPos2, allAttaches1, allAttaches2];
  }

  protected static getAllReactants(rules: RuleReaction[], monomers: string [], positions: [number, number, number][]):
    [string [], number [], number [], number []] {
    const allPos1: number [] = [];
    const allPos2: number [] = [];
    const rule: number [] = [];
    const ruleCount = rules.length;

    for (let i = 0; i < ruleCount; i++) {
      if (positions[i][0] == -1)
        continue;

      const firstMonomer = monomers[positions[i][0]];
      const secondMonomer = monomers[positions[i][1]];
      // monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstMonomer);
      // monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondMonomer);

      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(0)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(0)', '');
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(1)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(1)', '');
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(2)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(2)', '');
      monomers[positions[i][0]] = monomers[positions[i][0]].replace('(3)', '');
      monomers[positions[i][1]] = monomers[positions[i][1]].replace('(3)', '');

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      rule.push(positions[i][2]);
    }

    return [monomers, allPos1, allPos2, rule];
  }

  public check(throwError: boolean = false): string[] {
    const errors: string[] = [];

    const chainsMonomerCount = this.monomers.map((ch) => ch.length).reduce((acc, curr) => acc + curr, 0);
    if (this.mol.atoms.length !== chainsMonomerCount) {
      errors.push(`The mol atoms count ${this.mol.atoms.length} does not match ` +
        `the total number ${chainsMonomerCount} of chains' monomers.`);
    }

    const internalBondsCount = this.monomers.map((ch) => ch.length - 1).reduce((acc, curr) => acc + curr, 0);
    const chainsBondCount = internalBondsCount + this.linkages.length;
    if (this.mol.bonds.length !== chainsBondCount) {
      errors.push(`The mol bonds count ${this.mol.bonds.length} does not match ` +
        `the total number ${chainsBondCount} in- and inter-chain linkages.`);
    }

    let counter: number = 0;
    for (let spIdx = 0; spIdx < this.monomers.length; ++spIdx) {
      const chain = this.monomers[spIdx];
      for (let mIntIdx = 0; mIntIdx < chain.length; ++mIntIdx) {
        try {
          const m = chain[mIntIdx];
          const a = this.mol.atoms[counter];
          if (a.bio!.continuousId !== counter)
            errors.push(`Atom #${counter} has incorrect .bio.continuousId: ${a.bio!.continuousId}.`);
          if (a.elem !== m)
            errors.push(`Atom #${counter} elem: '${a.elem}' does not match chain monomer: '${m}'.`);
        } finally { counter++; }
      }
    }
    if (throwError && errors.length > 0)
      throw new Error(`Chain errors:\n${errors.map((e) => `  ${e}`).join('\n')}`);
    return errors;
  }
}
