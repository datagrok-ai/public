import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MolList, RDModule, RDMol, RDReaction, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';

const MAX_PRODUCTS = 10;

export function getEqualMolsIdxs(rdkit: RDModule, firstSet: string[], secondSet: string[]): [number, number] {
  let res: [number, number] = [-1, -1];
  let molsFirst: (RDMol | null) [] = [];
  let molsSecond: (RDMol | null) [] = [];

  try {
    molsFirst = firstSet.map((m) => rdkit.get_mol(m));
    molsSecond = secondSet.map((m) => rdkit.get_mol(m));

    outer:
    for (let i = 0; i < molsFirst.length; i++) {
      for (let j = 0; j < molsSecond.length; j++) {
        const match1 = molsFirst[i]?.get_substruct_match(molsSecond[j]!);
        const match2 = molsSecond[j]?.get_substruct_match(molsFirst[i]!);
        if (match1 !== '{}' && match2 !== '{}') {
          res = [i, j];
          break outer;
        }
      }
    }
  } catch (err: any) {
    console.log(err);
  } finally {
    for (let i = 0; i < molsFirst.length; i++)
      molsFirst[i]?.delete();
    for (let i = 0; i < molsFirst.length; i++)
      molsSecond[i]?.delete();
  }

  return res;
}

export function cutFragments(rdkit: RDModule, molecules: string[], fragment: string): string[][] {
  let rxn: RDReaction | null = null;
  const res: string[][] = new Array<string[]>(molecules.length);
  try {
    rxn = getFragmentCutReaction(rdkit, fragment);
    for (let i = 0; i < molecules.length; i++)
      res[i] = cutFragment(rdkit, rxn, molecules[i]);
  } catch (err: any) {
    throw err;
  } finally {
    rxn?.delete();
  }

  return res;
}

function cutFragment(rdkit: RDModule, rxn: RDReaction, molecule: string): string[] {
  let mols: MolList | null = null;
  let mol: RDMol | null = null;
  let rctns: RDReactionResult | null = null;
  const res: string [] = [];

  try {
    mols = new rdkit.MolList();
    mol = rdkit.get_mol(molecule);
    mols.append(mol);
    rctns = rxn.run_reactants(mols, MAX_PRODUCTS);
    const size = rctns.size();
    res.length = size;
    for (let i = 0; i < size; i++) {
      const element = rctns.get(i);
      res[i] = extractFirstProduct(element);
    }
  } catch (err: any) {
    throw err;
  } finally {
    mols?.delete();
    mol?.delete();
    rctns?.delete();
  }

  return res;
}

function extractFirstProduct(products: MolList): string {
  let molP: RDMol | null = null;

  let molBlock = '';

  try {
    molP = products.next();
    molBlock = molP?.get_molblock();//molP?.get_v3Kmolblock();//
  } catch (err: any) {

  } finally {
    molP?.delete();
  }

  return molBlock;
}

function getFragmentCutReaction(rdkit: RDModule, fragment: string): RDReaction {
  const reacSmarts = `${fragment}>>[*:1]`;
  return getReactionSmirks(rdkit, reacSmarts);
}

function getReactionSmirks(rdkit: RDModule, smirks: string): RDReaction {
  let rxn: RDReaction | null = null;
  try {
    rxn = rdkit.get_rxn(smirks);
    if (!rxn) throw new Error(`Invalid reaction '${smirks}'.`);
  } catch (err: any) {
    rxn?.delete();
    throw err;
  } finally {

  }
  return rxn;
}
