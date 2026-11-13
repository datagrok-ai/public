/**
 * Reactions over a column in the Crux workers (crux.worker.ts runs them): Deprotect's (analysis/deprotect.ts), each
 * molecule reacted as the main thread reacts it, a chunk per job, instead of molecule by molecule on the main thread.
 */
import type {MolList, RDModule, RDMol, RDReaction, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {hasNewLines} from '../../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../../utils/chem-constants';
import {WorkerJob, throwIfTrap} from '../crux-api';

/** The molecule's first product with `reaction` on `module` (its H atoms added, then removed) as a MOL block; the
 * molecule itself when it gives none; throws when the module reads no molecule from it. */
export function deprotectOne(module: RDModule, reaction: RDReaction, molSmiles: string): string {
  let mols: MolList | null = null;
  let mol: RDMol | null = null;
  let rctns: RDReactionResult | null = null;
  try {
    mol = module.get_mol(molSmiles);
    if (!mol)
      throw new Error(`no molecule read from ${molSmiles}`);
    mol.add_hs_in_place();
    mols = new module.MolList();
    mols.append(mol);
    rctns = reaction.run_reactants(mols, 1);
    const size = rctns.size();
    if (size === 0)
      return molSmiles;
    const products = rctns.get(0);
    try {
      const firstProduct = products.next();
      if (!firstProduct)
        return molSmiles;
      try {
        firstProduct.remove_hs_in_place();
        const deprotectedBlock = firstProduct.get_molblock();
        return deprotectedBlock ?? molSmiles;
      } finally {
        firstProduct.delete();
      }
    } finally {
      products.delete();
    }
  } finally {
    mols?.delete();
    mol?.delete();
    rctns?.delete();
  }
}

/** The jobs by name (the Crux worker's `op`). */
export const REACTION_JOBS: {[op: string]: WorkerJob} = {
  /**
   * Deprotect (analysis/deprotect.ts `addDeprotectedColumn`): each molecule's result on Crux, a SMILES over Chem's
   * length limit as it is; null for an empty cell and for a molecule Crux throws on (the main thread runs it, RDKit
   * behind Crux). A reaction Crux reads none from fails the job (the main thread does the column).
   */
  'deprotect': ({module}: {module: RDModule}, molecules: (string | null)[], smirks: string): (string | null)[] => {
    const reaction = module.get_rxn(smirks);
    if (!reaction)
      throw new Error(`Crux reads no reaction from ${smirks}`);
    try {
      return molecules.map((molSmiles) => {
        if (molSmiles == null)
          return null;
        if (molSmiles && !hasNewLines(molSmiles) && molSmiles.length > MAX_SMILES_LENGTH)
          return molSmiles;
        try {
          return deprotectOne(module, reaction, molSmiles);
        } catch (e) {
          throwIfTrap(e);
          return null;
        }
      });
    } finally {
      reaction.delete();
    }
  },
};
