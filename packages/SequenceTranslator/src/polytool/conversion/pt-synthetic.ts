import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {Rules, RuleReaction, getMonomerPairs} from './pt-rules';
import {InvalidReactionError, MonomerNotFoundError} from '../types';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

import * as grok from 'datagrok-api/grok';

import wu from 'wu';
import {PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {IMonomerLib, IMonomerLibBase, Monomer, MonomerLibData, RGroup} from '@datagrok-libraries/bio/src/types';
import {RDModule, RDMol, RDReaction, MolList, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {HELM_REQUIRED_FIELD as REQ,
  HELM_OPTIONAL_FIELDS as OPT, HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';

export async function getOverriddenLibrary(rules: Rules): Promise<IMonomerLibBase> {
  const monomerLibHelper = await getMonomerLibHelper();
  const systemMonomerLib = monomerLibHelper.getMonomerLib();

  const rdkit = await getRdKitModule();
  const argLib: { [symbol: string]: Monomer } = {};

  let names: string [] = [];
  let monomers: Monomer [] = [];

  for (let i = 0; i < rules.reactionRules.length; i++) {
    try {
      [names, monomers] = getNewMonomers(rdkit, systemMonomerLib, rules.reactionRules[i]);
    } catch (e: any) {
      names = [];
      monomers = [];
      console.error(e);
      grok.shell.warning(e);
    } finally {
      for (let j = 0; j < names.length; j ++)
        argLib[names[j]] = monomers[j];
    }
  }

  const overrideMonomerLibData: MonomerLibData = {[PolymerTypes.PEPTIDE]: argLib};
  const overriddenMonomerLib = systemMonomerLib.override(overrideMonomerLibData,
    'ST-PT-reactions.' + wu.repeat(1).map(() => Math.floor((Math.random() * 36))
      .toString(36)).take(4).toArray().join(''));
  return overriddenMonomerLib;
}

export function getNewMonomers(rdkit: RDModule, mLib: IMonomerLib, rule: RuleReaction): [string[], Monomer[]] {
  const reacSmarts = rule.reaction;
  const monomerName = rule.name;
  const totalLength = rule.firstMonomers.length + rule.secondMonomers.length + 1;
  const reactionMembers = rule.reaction.split('>>');
  const reactants = reactionMembers[0].split('.');

  const monomerNames = new Array<string>(totalLength);
  const resMonomers = new Array<Monomer>(totalLength);
  const rawMonomers = new Array<string>(totalLength);

  const mainRxn = getReactionSmirks(rdkit, reacSmarts);

  const reacCutFirst = `${reactants[0]}>>[C:1]`;
  const reacCutSecond = `${reactants[1]}>>[C:2]`;

  const rxnCutFirst = getReactionSmirks(rdkit, reacCutFirst);
  const rxnCutSecond = getReactionSmirks(rdkit, reacCutSecond);

  let counter = 0;
  for (let i = 0; i < rule.firstMonomers.length; i++) {
    const monomer = mLib.getMonomer('PEPTIDE', rule.firstMonomers[i]);
    if (!monomer) throw new MonomerNotFoundError('PEPTIDE', rule.firstMonomers[i]);

    const sMolBlock = cutReactant(rdkit, monomer.molfile, rxnCutFirst, monomer.name);
    rawMonomers[counter] = sMolBlock;
    monomerNames[counter] = `${monomer.symbol}_${monomerName}`;
    counter++;
  }
  for (let i = 0; i < rule.secondMonomers.length; i++) {
    const monomer = mLib.getMonomer('PEPTIDE', rule.secondMonomers[i]);
    if (!monomer) throw new MonomerNotFoundError('PEPTIDE', rule.secondMonomers[i]);

    const sMolBlock = cutReactant(rdkit, monomer.molfile, rxnCutSecond, monomer.name);
    rawMonomers[counter] = sMolBlock;
    monomerNames[counter] = `${monomer.symbol}_${monomerName}`;
    counter++;
  }

  // const monomer1 = mLib.getMonomer('PEPTIDE', rule.firstMonomers[0]);
  // const monomer2 = mLib.getMonomer('PEPTIDE', rule.secondMonomers[0]);
  let mol: RDMol | null = null;
  let pMolblock = '';
  try {
    mol = rdkit.get_mol(reactionMembers[1]);
    pMolblock = mol?.get_molblock();//molP?.get_v3Kmolblock();//
  } catch (err: any) {
    const [errMsg, _errStack] = errInfo(err);
    grok.shell.error(`Can not assemble monomer '${monomerName}': ${errMsg}.`);
    throw err;
  } finally {
    mol?.delete();
  }
  //pMolblock = reactionMembers[1]//cutProduct(rdkit, monomer1?.molfile, monomer2?.molfile, mainRxn, monomerName);
  rawMonomers[counter] = pMolblock;
  monomerNames[counter] = monomerName;

  //after RDKit works
  for (let i = 0; i < totalLength - 1; i ++)
    rawMonomers[i] = rawMonomers[i].replace('0.0000 C   ', '0.0000 R#  ').replace('M  RGP  2', 'M  RGP  3   1   3');
  rawMonomers[totalLength - 1] = modProduct(rawMonomers[totalLength - 1]);

  for (let i = 0; i < totalLength; i ++) {
    const isProduct = i == totalLength - 1 ? true : false;
    const resMonomer: Monomer = {
      [REQ.SYMBOL]: monomerNames[i],
      [REQ.NAME]: monomerNames[i],
      [REQ.MOLFILE]: rawMonomers[i],
      [REQ.AUTHOR]: '',
      [REQ.ID]: 0,
      [REQ.RGROUPS]: getNewGroups(isProduct),
      [REQ.SMILES]: '',
      [REQ.POLYMER_TYPE]: 'PEPTIDE',
      [REQ.MONOMER_TYPE]: 'Backbone',
      [REQ.CREATE_DATE]: null,
      // // @ts-ignore
      // lib: {source: 'Reaction'},
    };

    resMonomer[OPT.META] = Object.assign(resMonomer[OPT.META] ?? {},
      {'colors': {'default': {line: '#2083D5', text: '#2083D5', background: '#F2F2F5'}}});

    monomerNames[i] = monomerNames[i];
    resMonomers[i] = resMonomer;
  }
  mainRxn.delete();
  rxnCutFirst.delete();
  rxnCutSecond.delete();
  return [monomerNames, resMonomers];
}

function modProduct(product: string): string {
  const fullMol = product.replace('M  RAD', 'M  RGP');
  const lines = fullMol.split('\n');
  const natoms = Number(lines[3].substring(0, 3));
  const nbonds = Number(lines[3].substring(3, 6));
  const rgpShift = 4 + natoms + nbonds;

  const fAtom = Number(lines[rgpShift].substring(9, 13));
  const sAtom = Number(lines[rgpShift].substring(17, 21));

  lines[rgpShift] = lines[rgpShift].substring(0, 13) + '   1' + lines[rgpShift].substring(17, 21) + '   2';

  lines[3 + fAtom] = lines[3 + fAtom].substring(0, 30) + ' R#  0  0  0  0  0  0  0  0  0  0  0  0';
  lines[3 + sAtom] = lines[3 + sAtom].substring(0, 30) + ' R#  0  0  0  0  0  0  0  0  0  0  0  0';

  let res = '';
  for (let i = 0; i < lines.length; i++)
    res += lines[i] + '\n';

  return res;
}

function cutReactant(rdkit: RDModule, reactant: string, rxn: RDReaction, monomerName: string) {
  let mols: MolList | null = null;
  let mol: RDMol | null = null;
  let rctns: RDReactionResult | null = null;
  let molP: RDMol | null = null;
  let molBlock = '';

  try {
    mols = new rdkit.MolList();
    mol = rdkit.get_mol(reactant);
    mols.append(mol);
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
    mols?.delete();
    mol?.delete();
    rctns?.delete();
    molP?.delete();
  }

  return molBlock;
}

function getReactionSmirks(rdkit: RDModule, smirks: string): RDReaction {
  let rxn: RDReaction | null = null;
  try {
    rxn = rdkit.get_rxn(smirks);
    if (!rxn) throw new InvalidReactionError(smirks);
  } catch (err: any) {
    rxn?.delete();
    const [errMsg, _errStack] = errInfo(err);
    grok.shell.error(`Can not assemble monomer '${smirks}': ${errMsg}.`);
    throw err;
  } finally {

  }
  return rxn;
}

function getNewGroups(isProduct: boolean): RGroup[] {
  if (isProduct) {
    const groups = new Array<RGroup>(2);
    const group1: RGroup = {
      //@ts-ignore
      [HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE]: '[*:1][H]',
      [HELM_RGROUP_FIELDS.ALTERNATE_ID]: 'R1-H',
      [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: 'H',
      [HELM_RGROUP_FIELDS.LABEL]: 'R1'
    };
    const group2: RGroup = {
      //@ts-ignore
      [HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE]: '[*:2][H]',
      [HELM_RGROUP_FIELDS.ALTERNATE_ID]: 'R2-H',
      [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: 'H',
      [HELM_RGROUP_FIELDS.LABEL]: 'R2'
    };

    groups[0] = group1;
    groups[1] = group2;

    return groups;
  } else {
    const group: RGroup = {
      //@ts-ignore
      [HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE]: '[*:3][H]',
      [HELM_RGROUP_FIELDS.ALTERNATE_ID]: 'R3-H',
      [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: 'H',
      [HELM_RGROUP_FIELDS.LABEL]: 'R3'
    };

    return [group];
  }
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
