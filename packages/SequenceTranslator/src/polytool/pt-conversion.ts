import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {HelmTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {IMonomerLib, IMonomerLibBase, Monomer, MonomerLibData, RGroup} from '@datagrok-libraries/bio/src/types';
import {RDModule, RDMol, RDReaction, MolList, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {HELM_REQUIRED_FIELD as REQ, HELM_OPTIONAL_FIELDS as OPT, HELM_RGROUP_FIELDS, HELM_OPTIONAL_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {HelmMol, HelmType, JSDraw2ModuleType, OrgType} from '@datagrok-libraries/bio/src/helm/types';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {Rules, RuleLink, RuleReaction} from './pt-rules';
import {InvalidReactionError, MonomerNotFoundError} from './types';

import {_package} from '../package';

declare const JSDraw2: JSDraw2ModuleType;
declare const org: OrgType;


export class Chain {
  linkages: { fChain: number, sChain: number, fMonomer: number, sMonomer: number, fR: number, sR: number }[];
  monomers: string[][];
  mol: HelmMol | null;

  constructor(
    monomers: string[][],
    linkages: { fChain: number, sChain: number, fMonomer: number, sMonomer: number, fR: number, sR: number }[],
    mol: HelmMol | null) {
    this.linkages = linkages;
    this.monomers = monomers;
    this.mol = mol;
  }

  static fromHelm(helm: string) {
    const fragmentation = helm.split('$');
    const rawFragments = fragmentation[0].split('|');
    const rawLinkages = fragmentation[1].split('|');

    const monomers = new Array<Array<string>>(rawFragments.length);
    const linkages: {
      fChain: number,
      sChain: number,
      fMonomer: number,
      sMonomer: number,
      fR: number,
      sR: number
    }[] = [];

    //HELM parsing
    for (let i = 0; i < rawFragments.length; i++) {
      const idxStart = rawFragments[i].indexOf('{');
      const idxEnd = rawFragments[i].indexOf('}');

      monomers[i] = rawFragments[i].slice(idxStart + 1, idxEnd).split('.');
    }

    //HELM parsing
    for (let i = 0; i < rawLinkages.length; i++) {
      if (rawLinkages[i] !== '' && rawLinkages[i] !== 'V2.0') {
        const rawData = rawLinkages[i].split(',');
        const seq1 = (rawData[0].replace('PEPTIDE', '') as unknown as number) - 1;
        const seq2 = (rawData[1].replace('PEPTIDE', '') as unknown as number) - 1;
        const rawDataConnctions = rawData[2].split('-');
        const rawDataConnction1 = rawDataConnctions[0].split(':');
        const rawDataConnction2 = rawDataConnctions[1].split(':');

        linkages.push({
          fChain: seq1,
          sChain: seq2,
          fMonomer: rawDataConnction1[0] as unknown as number,
          sMonomer: rawDataConnction2[0] as unknown as number,
          fR: rawDataConnction1[1].replace('R', '') as unknown as number,
          sR: rawDataConnction2[1].replace('R', '') as unknown as number,
        });
      }
    }

    return new Chain(monomers, linkages, null);
  }

  static fromNotation(sequence: string, rules: Rules): Chain {
    const heterodimerCode = rules.heterodimerCode;
    const homodimerCode = rules.homodimerCode;
    const mainFragments: string[] = [];

    const linkages: {
      fChain: number,
      sChain: number,
      fMonomer: number,
      sMonomer: number,
      fR: number,
      sR: number
    }[] = [];

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

    for (let i = 0; i < monomers.length; i++) {
      const linkedPositions = this.getLinkedPositions(monomers[i], rules.reactionRules);
      const [monomersCycled, allPos1, allPos2, ruleN] =
        this.getAllReactants(rules.reactionRules, monomers[i], linkedPositions);

      if (allPos1.length >= 1) {
        const ch1 = new Array<string>(allPos2[0] - 1);
        const ch2 = new Array<string>(monomersCycled.length - allPos2[0]);
        for (let j = 0; j < allPos2[0] - 1; j++)
          ch1[j] = monomersCycled[j];

        for (let j = allPos2[0]; j < monomersCycled.length; j++)
          ch2[j - allPos2[0]] = monomersCycled[j];


        ch1[allPos1[0] - 1] = rules.reactionRules[ruleN[0]].name;

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
        monomersAll.push(monomers[i]);
      }
    }

    const chain = new Chain(monomersAll, linkages, null);
    return chain;
  }

  static async parseNotation(sequence: string): Promise<Chain> {
    const mainFragments: string[][] = [];

    const linkages: {
      fChain: number,
      sChain: number,
      fMonomer: number,
      sMonomer: number,
      fR: number,
      sR: number
    }[] = [];


    const hh = await getHelmHelper();
    // const sampleHwe = hh.createHelmWebEditor();
    // sampleHwe.editor.setHelm('PEPTIDE1{R.P.D.[meI]}$$$$');

    const resHwe = hh.createHelmWebEditor();
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
      for (let j = 0; j < splMonomers.length; j++) {
        const monomer = splMonomers[j].replace('{', '').replace('}', '');
        if (monomer !== '') {
          monomers[j] = monomer;
          counter++;
        } else {
          linkages.push({fChain: i, sChain: i + 1, fMonomer: counter, sMonomer: counter + 1, fR: 1, sR: 1});
        }
      }

      mainFragments.push(monomers);
    }

    counter = 0;
    const p = new JSDraw2.Point(0, 0);
    for (let i = 0; i < mainFragments.length; i++) {
      for (let j = 0; j < mainFragments[i].length; j++) {
        if (!!mainFragments[i][j]) {
          const elem = mainFragments[i][j];
          const bio = {type: HelmTypes.AA, i: i, j: j};
          const atom = new JSDraw2.Atom<HelmType>(p, elem, bio);
          resMol.addAtom(atom);

          if (j !== 0) {
            const atom1 = resMol.atoms[counter - 1];
            const atom2 = resMol.atoms[counter];
            const bond = new JSDraw2.Bond<HelmType>(atom1, atom2);
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
      const bond = new JSDraw2.Bond<HelmType>(atom1, atom2);
      bond.r1 = linkages[i].fR;
      bond.r2 = linkages[i].sR;
      resMol.addBond(bond);
    }

    const chain = new Chain(mainFragments, linkages, resMol);
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

  getNotationHelm(): string {
    return this.getHelm();
  }

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
      helm += `${this.linkages[i].fMonomer}:R${this.linkages[i].fR}-`;
      helm += `${this.linkages[i].sMonomer}:R${this.linkages[i].sR}`;
    }

    helm += '$$$';
    return helm;
  }

  getNotation(): string {
    const atoms = this.mol!.atoms;
    const bonds = this.mol!.bonds;
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
    const result: [number, number, number][] = new Array<[number, number, number]>(rules.length);

    for (let i = 0; i < rules.length; i++) {
      let firstFound = false;
      let secondFound = false;
      let firstIsFirst = false;
      let firstEntryIndex = -1;
      let secondEntryIndex = -1;
      const add = `(${rules[i].code})`;
      for (let j = 0; j < monomers.length; j++) {
        if (monomers[j].includes(add)) {
          if (firstFound) {
            if (firstIsFirst && monomers[j] == rules[i].secondMonomer + add) {
              secondFound = true;
              secondEntryIndex = j;
              break;
            } else if (!firstIsFirst && monomers[j] == rules[i].firstMonomer + add) {
              secondFound = true;
              secondEntryIndex = j;
              break;
            } else {
              continue;
            }
          } else {
            if (monomers[j] == rules[i].firstMonomer + add) {
              firstFound = true;
              firstIsFirst = true;
              firstEntryIndex = j;
            } else if (monomers[j] == rules[i].secondMonomer + add) {
              firstFound = true;
              firstIsFirst = false;
              firstEntryIndex = j;
            } else {
              continue;
            }
          }
        }
      }

      if (!(firstFound && secondFound))
        result[i] = [-1, -1, -1];
      else if (firstIsFirst)
        result[i] = [firstEntryIndex, secondEntryIndex, i];
      else
        result[i] = [secondEntryIndex, firstEntryIndex, i];
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
      monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstSubstitution);
      monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondSubstitution);

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
      monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstMonomer);
      monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondMonomer);

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      rule.push(positions[i][2]);
    }

    return [monomers, allPos1, allPos2, rule];
  }
}

/** The main PolyTool convert engine. Returns list of Helms. Covered with tests. */
export function doPolyToolConvert(sequences: string[], rules: Rules): string[] {
  const helms = new Array<string>(sequences.length);
  for (let i = 0; i < sequences.length; i++) {
    try {
      if (sequences[i] == null) { helms[i] = ''; } else {
        const chain = Chain.fromNotation(sequences[i], rules);
        helms[i] = chain.getHelm();
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
      helms[i] = '';
    }
  }
  return helms;
}

function getMonomersMolBlocks(monomer1: Monomer, monomer2: Monomer): [string, string] {
  const mb1 = monomer1.molfile;
  let mb2 = monomer2.molfile;
  const addGroups = monomer1.rgroups.length;

  //mol v2000 monomer
  const rgpIdx = mb2.indexOf('M  RGP');
  if (rgpIdx !== -1) {
    const groupsCountStr = mb2.substring(rgpIdx + 6, rgpIdx + 9);
    const groupsCount = Number(groupsCountStr);

    for (let i = 0; i < groupsCount; i++) {
      const start = rgpIdx + 9 + 4 + i * 8;
      const end = rgpIdx + 9 + 8 + i * 8;
      const rGroupSpecifier = mb2.substring(start, end);
      const groupPosition = Number(rGroupSpecifier) + addGroups;
      const digits = Math.floor(Math.log10(groupPosition) + 1);
      const newSpecifier = ' '.repeat(4 - digits) + String(groupPosition);
      mb2 = mb2.substring(0, start) + newSpecifier + mb2.substring(end, mb2.length);
    }
  }

  //TODO: same for v3000 monomer

  return [mb1, mb2];
}

function getSyntheticMolBlock(rdkit: RDModule, reaction: string,
  mb1: string, mb2: string, monomerName: string): string {
  let rxn: RDReaction | null = null;
  let mols: MolList | null = null;
  let mol1: RDMol | null = null;
  let mol2: RDMol | null = null;
  let rctns: RDReactionResult | null = null;
  let molP: RDMol | null = null;
  let molBlock = '';

  try {
    rxn = rdkit.get_rxn(reaction);
    if (!rxn) throw new InvalidReactionError(reaction);
    mols = new rdkit.MolList();
    mol1 = rdkit.get_mol(mb1!);
    mol2 = rdkit.get_mol(mb2!);
    mols.append(mol1!);
    mols.append(mol2!);

    rctns = rxn.run_reactants(mols, 1);
    //const size = rctns.size();
    const element = rctns.get(0);

    molP = element.next();
    molBlock = molP?.get_molblock();//molP?.get_v3Kmolblock();//
  } catch (err: any) {
    const [errMsg, _errStack] = errInfo(err);
    grok.shell.error(`Can not assemble monomer '${monomerName}': ${errMsg}.`);
    throw err;
  } finally {
    rxn?.delete();
    mols?.delete();
    mol1?.delete();
    mol2?.delete();
    rctns?.delete();
    molP?.delete();
  }

  return molBlock;
}

function getNewGroups(monomer1: Monomer, monomer2: Monomer): RGroup[] {
  const groups = new Array<RGroup>(monomer1?.rgroups.length! + monomer2?.rgroups.length!);
  const length1 = monomer1?.rgroups.length!;
  const length2 = monomer2?.rgroups.length!;

  for (let i = 0; i < length1; i++)
    groups[i] = monomer1?.rgroups[i]!;

  for (let i = 0; i < length2; i++) {
    const rGroupSpecifier = monomer2?.rgroups[i]!.label.replace('R', '');
    const groupPosition = Number(rGroupSpecifier) + length1;
    const group: RGroup = {
      //@ts-ignore
      [HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE]: monomer2?.rgroups[i].capGroupSMILES.replace(rGroupSpecifier, String(groupPosition)),
      [HELM_RGROUP_FIELDS.ALTERNATE_ID]: monomer2?.rgroups[i].alternateId.replace(rGroupSpecifier, String(groupPosition)),
      [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: monomer2?.rgroups[i].capGroupName,
      [HELM_RGROUP_FIELDS.LABEL]: monomer2?.rgroups[i].label.replace(rGroupSpecifier, String(groupPosition)),
    };

    groups[i + length1] = group;
  }

  return groups;
}

export function getNewMonomer(rdkit: RDModule, mLib: IMonomerLib, rule: RuleReaction): [string, Monomer] {
  const reacSmarts = rule.reaction;
  const monomerName = rule.name;

  const monomer1 = mLib.getMonomer('PEPTIDE', rule.firstMonomer);
  if (!monomer1) throw new MonomerNotFoundError('PEPTIDE', rule.firstMonomer);
  const monomer2 = mLib.getMonomer('PEPTIDE', rule.secondMonomer);
  if (!monomer2) throw new MonomerNotFoundError('PEPTIDE', rule.secondMonomer);

  const [mb1, mb2] = getMonomersMolBlocks(monomer1!, monomer2!);
  const molBlock = getSyntheticMolBlock(rdkit, reacSmarts, mb1, mb2, monomerName);
  const groups: RGroup[] = getNewGroups(monomer1!, monomer2!);

  const resMonomer: Monomer = {
    [REQ.SYMBOL]: monomerName,
    [REQ.NAME]: monomerName,
    [REQ.MOLFILE]: molBlock,
    [REQ.AUTHOR]: '',
    [REQ.ID]: 0,
    [REQ.RGROUPS]: groups,
    [REQ.SMILES]: '',
    [REQ.POLYMER_TYPE]: 'PEPTIDE',
    [REQ.MONOMER_TYPE]: 'Backbone',
    [REQ.CREATE_DATE]: null,

    // // @ts-ignore
    // lib: {source: 'Reaction'},
  };

  resMonomer[OPT.META] = Object.assign(resMonomer[OPT.META] ?? {},
    {'colors': {'default': {line: '#2083D5', text: '#2083D5', background: '#F2F2F5'}}});

  return [monomerName, resMonomer];
}

export async function getOverriddenLibrary(rules: Rules): Promise<IMonomerLibBase> {
  const monomerLibHelper = await getMonomerLibHelper();
  const systemMonomerLib = monomerLibHelper.getMonomerLib();

  const rdkit = await getRdKitModule();
  const argLib: { [symbol: string]: Monomer } = {};

  for (let i = 0; i < rules.reactionRules.length; i++) {
    const [name, monomer] = getNewMonomer(rdkit, systemMonomerLib, rules.reactionRules[i]);
    argLib[name] = monomer;
  }

  const overrideMonomerLibData: MonomerLibData = {[PolymerTypes.PEPTIDE]: argLib};
  const overriddenMonomerLib = systemMonomerLib.override(overrideMonomerLibData,
    'ST-PT-reactions.' + wu.repeat(1).map(() => Math.floor((Math.random() * 36)).toString(36)).take(4).toArray().join(''));
  return overriddenMonomerLib;
}
