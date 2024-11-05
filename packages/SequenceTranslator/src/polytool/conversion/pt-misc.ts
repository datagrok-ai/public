import * as grok from 'datagrok-api/grok';

import wu from 'wu';
import {PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {IMonomerLib, IMonomerLibBase, Monomer, MonomerLibData, RGroup} from '@datagrok-libraries/bio/src/types';
import {RDModule, RDMol, RDReaction, MolList, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {HELM_REQUIRED_FIELD as REQ,
  HELM_OPTIONAL_FIELDS as OPT, HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {Rules, RuleReaction, getMonomerPairs} from './pt-rules';
import {InvalidReactionError, MonomerNotFoundError} from '../types';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

/** Gets 0-based in-index (simple polymer) of out-index (continuous) {@link idx} */
export function getInnerIdx(outIdx: number, monomers: string[][]): [number, number] {
  // let prevSpCount = 0;
  // for (let spI = 0; spI < monomers.length && idx >= (prevSpCount + monomers[spI].length); ++spI)
  //   prevSpCount += monomers[spI].length;
  // return idx - prevSpCount;
  let inIdx = outIdx;
  let spIdx: number;
  for (spIdx = 0; spIdx < monomers.length && inIdx >= monomers[spIdx].length; ++spIdx)
    inIdx -= monomers[spIdx].length;
  return [inIdx, spIdx];
}

/** Gets 0-based out-index of 0-based in-index {@link inIdx} monomer of simple polymer {@link spIdx} */
export function getOuterIdx(inIdx: number, spIdx: number, monomers: string[][]): number {
  let outIdx = 0;
  for (let i = 0; i < spIdx; ++i)
    outIdx += monomers[i].length;
  return outIdx + inIdx;
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
      [HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE]: monomer2?.rgroups[i].capGroupSMILES
        .replace(rGroupSpecifier, String(groupPosition)),
      [HELM_RGROUP_FIELDS.ALTERNATE_ID]: monomer2?.rgroups[i].alternateId
        .replace(rGroupSpecifier, String(groupPosition)),
      [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: monomer2?.rgroups[i].capGroupName,
      [HELM_RGROUP_FIELDS.LABEL]: monomer2?.rgroups[i].label.replace(rGroupSpecifier, String(groupPosition)),
    };

    groups[i + length1] = group;
  }

  return groups;
}

export function getNewMonomers(rdkit: RDModule, mLib: IMonomerLib, rule: RuleReaction): [string[], Monomer[]] {
  const reacSmarts = rule.reaction;
  const monomerName = rule.name;

  const [firstMonomers, secondMonomers] = getMonomerPairs(rule);

  const monomerNames = new Array<string>(firstMonomers.length);
  const resMonomers = new Array<Monomer>(firstMonomers.length);

  for (let i = 0; i < firstMonomers.length; i ++) {
    const monomer1 = mLib.getMonomer('PEPTIDE', firstMonomers[i]);
    if (!monomer1) throw new MonomerNotFoundError('PEPTIDE', firstMonomers[i]);
    const monomer2 = mLib.getMonomer('PEPTIDE', secondMonomers[i]);
    if (!monomer2) throw new MonomerNotFoundError('PEPTIDE', secondMonomers[i]);

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

    monomerNames[i] = monomerName;
    resMonomers[i] = resMonomer;
  }

  return [monomerNames, resMonomers];
}

export async function getOverriddenLibrary(rules: Rules): Promise<IMonomerLibBase> {
  const monomerLibHelper = await getMonomerLibHelper();
  const systemMonomerLib = monomerLibHelper.getMonomerLib();

  const rdkit = await getRdKitModule();
  const argLib: { [symbol: string]: Monomer } = {};

  for (let i = 0; i < rules.reactionRules.length; i++) {
    const [names, monomers] = getNewMonomers(rdkit, systemMonomerLib, rules.reactionRules[i]);
    for (let j = 0; j < names.length; j ++)
      argLib[names[j]] = monomers[j];
  }

  const overrideMonomerLibData: MonomerLibData = {[PolymerTypes.PEPTIDE]: argLib};
  const overriddenMonomerLib = systemMonomerLib.override(overrideMonomerLibData,
    'ST-PT-reactions.' + wu.repeat(1).map(() => Math.floor((Math.random() * 36))
      .toString(36)).take(4).toArray().join(''));
  return overriddenMonomerLib;
}
