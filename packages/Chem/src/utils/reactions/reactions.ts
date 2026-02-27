/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {MolList, RDModule, RDMol, RDReaction, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitService} from '../chem-common-rdkit';
import {ReactionResult, TwoComponentReactionResult} from './types';
import {SALT_AND_WATER_FRAGMENTS} from './consts';

// ---- Salt / Water removal ----

export function removeWaterAndSaltsSingle(molSmiles: string): string {
  if (!molSmiles || molSmiles.indexOf('.') === -1)
    return molSmiles ?? '';
  return molSmiles.split('.').filter((s) => !SALT_AND_WATER_FRAGMENTS.includes(s)).join('.');
}

export async function removeWaterAndSalts(mols: string[]): Promise<string[]> {
  const canonical = await (await getRdKitService()).convertMolNotation(mols, grok.chem.Notation.Smiles);
  return canonical.map((s) => removeWaterAndSaltsSingle(s));
}

// ---- Variable resolution ----

/** Substitute ${variable} placeholders in a SMARTS string with provided values. */
export function resolveReactionVariables(smarts: string, variables: Record<string, any>): string {
  let resolved = smarts;
  for (const [key, value] of Object.entries(variables))
    resolved = resolved.replace(`\${${key}}`, String(value));
  return resolved;
}

// ---- R-group correction ----

const WRONG_R_GROUPS = [1, 2, 3, 4, 5, 6].map((i) => `[*${i}]`);

function correctRGroups(smiles: string): string {
  let result = smiles;
  WRONG_R_GROUPS.forEach((rg, j) => {
    if (result.includes(rg))
      result = result.replaceAll(rg, `[*:${j + 1}]`);
  });
  return result;
}

// ---- Core: run a single RDKit reaction on one molecule ----

function runSingleReaction(
  rxn: RDReaction, rdkit: RDModule, smiles: string, maxProducts: number = 2,
): string | null {
  let mol: RDMol | null = null;
  let molList: MolList | null = null;
  let result: RDReactionResult | null = null;
  try {
    mol = rdkit.get_mol(smiles);
    molList = new rdkit.MolList();
    molList.append(mol);
    result = rxn.run_reactants(molList, maxProducts);
    if (result.size() === 0)
      return null;
    const productList = result.get(0);
    const productMol = productList.next();
    return productMol?.get_smiles() ?? null;
  } finally {
    try {result?.delete();} catch (_) {/* noop */}
    try {molList?.delete();} catch (_) {/* noop */}
    try {mol?.delete();} catch (_) {/* noop */}
  }
}

// ---- Core: run a two-component RDKit reaction on a molecule pair ----

function runTwoComponentSingle(
  rxn: RDReaction, rdkit: RDModule, smiles1: string, smiles2: string, maxProducts: number = 2,
): string | null {
  let mol1: RDMol | null = null;
  let mol2: RDMol | null = null;
  let molList: MolList | null = null;
  let result: RDReactionResult | null = null;
  try {
    mol1 = rdkit.get_mol(smiles1);
    mol2 = rdkit.get_mol(smiles2);
    molList = new rdkit.MolList();
    molList.append(mol1);
    molList.append(mol2);
    result = rxn.run_reactants(molList, maxProducts);
    if (result.size() === 0)
      return null;
    const productList = result.get(0);
    const productMol = productList.next();
    return productMol?.get_smiles() ?? null;
  } finally {
    try {result?.delete();} catch (_) {/* noop */}
    try {molList?.delete();} catch (_) {/* noop */}
    try {mol2?.delete();} catch (_) {/* noop */}
    try {mol1?.delete();} catch (_) {/* noop */}
  }
}

// ---- Batch UI yielding helper ----

async function yieldIfNeeded(i: number, batchSize: number = 50): Promise<void> {
  if (i % batchSize === 0)
    await new Promise<void>((resolve) => setTimeout(resolve, 0));
}

// ==========================================
//  TRANSFORMATION REACTION (single column)
// ==========================================

/**
 * Run a transformation reaction on a list of molecules.
 * The reaction SMARTS should have exactly one reactant.
 */
export async function runTransformationReaction(
  rdkit: RDModule,
  reactionSmarts: string,
  mols: string[],
  options: {
    removeSaltsAndWater?: boolean;
    resolvedVariables?: Record<string, any>;
    showInfo?: boolean;
  } = {},
): Promise<ReactionResult> {
  const {removeSaltsAndWater = true, resolvedVariables, showInfo = false} = options;

  let resolvedSmarts = reactionSmarts;
  if (resolvedVariables)
    resolvedSmarts = resolveReactionVariables(resolvedSmarts, resolvedVariables);

  const rxn = rdkit.get_rxn(resolvedSmarts);
  if (!rxn)
    throw new Error(`Invalid reaction SMARTS: "${resolvedSmarts}".`);

  try {
    const canonical = await (await getRdKitService()).convertMolNotation(mols, grok.chem.Notation.Smiles);
    const cleaned = removeSaltsAndWater ? canonical.map((s) => removeWaterAndSaltsSingle(s)) : canonical;
    const products: string[] = new Array(mols.length);

    let errors = 0;
    let noProducts = 0;
    let emptyMolecules = 0;

    for (let i = 0; i < cleaned.length; i++) {
      await yieldIfNeeded(i);
      const smiles = cleaned[i] ?? mols[i];
      if (!smiles) {
        emptyMolecules++;
        products[i] = '';
        continue;
      }
      try {
        const product = runSingleReaction(rxn, rdkit, smiles);
        if (product == null) {
          noProducts++;
          products[i] = smiles;
        } else
          products[i] = correctRGroups(product);
      } catch {
        errors++;
        products[i] = smiles;
      }
    }

    if (showInfo) {
      const success = mols.length - errors - noProducts - emptyMolecules;
      grok.shell.info(
        `Reaction complete: ${success} products, ${errors} errors, ${noProducts} no reaction, ${emptyMolecules} empty.`,
        {timeout: 10},
      );
    }

    return {products, errors, noProducts, emptyMolecules, reactionSmarts: resolvedSmarts};
  } finally {
    try {rxn.delete();} catch (_) {/* noop */}
  }
}

// ==========================================
//  TWO-COMPONENT REACTION (two columns)
// ==========================================

/**
 * Run a two-component reaction between two lists of molecules.
 *
 * Modes:
 *  - 'pairwise': react row i from list 1 with row i from list 2 (uses shorter length)
 *  - 'matrix': react every combination (N*M products)
 */
export async function runTwoComponentReaction(
  rdkit: RDModule,
  reactionSmarts: string,
  reactants1: string[],
  reactants2: string[],
  options: {
    mode?: 'pairwise' | 'matrix';
    removeSaltsAndWater?: boolean;
    resolvedVariables?: Record<string, any>;
    maxMatrixProducts?: number;
    showInfo?: boolean;
  } = {},
): Promise<TwoComponentReactionResult> {
  const {mode = 'pairwise', removeSaltsAndWater = true, resolvedVariables, maxMatrixProducts = 100_000, showInfo = false} = options;

  let resolvedSmarts = reactionSmarts;
  if (resolvedVariables)
    resolvedSmarts = resolveReactionVariables(resolvedSmarts, resolvedVariables);

  const rxn = rdkit.get_rxn(resolvedSmarts);
  if (!rxn)
    throw new Error(`Invalid reaction SMARTS: "${resolvedSmarts}".`);

  try {
    const service = await getRdKitService();
    const canon1 = await service.convertMolNotation(reactants1, grok.chem.Notation.Smiles);
    const canon2 = await service.convertMolNotation(reactants2, grok.chem.Notation.Smiles);
    const clean1 = removeSaltsAndWater ? canon1.map((s) => removeWaterAndSaltsSingle(s)) : canon1;
    const clean2 = removeSaltsAndWater ? canon2.map((s) => removeWaterAndSaltsSingle(s)) : canon2;

    if (mode === 'pairwise')
      return await runPairwise(rxn, rdkit, clean1, clean2, reactants1, reactants2, showInfo);
    else
      return await runMatrix(rxn, rdkit, clean1, clean2, reactants1, reactants2, maxMatrixProducts, showInfo);
  } finally {
    try {rxn.delete();} catch (_) {/* noop */}
  }
}

async function runPairwise(
  rxn: RDReaction, rdkit: RDModule,
  clean1: string[], clean2: string[],
  orig1: string[], orig2: string[],
  showInfo: boolean,
): Promise<TwoComponentReactionResult> {
  const len = Math.min(clean1.length, clean2.length);
  const r1: string[] = [];
  const r2: string[] = [];
  const prods: string[] = [];
  let errors = 0;
  let noProducts = 0;

  for (let i = 0; i < len; i++) {
    await yieldIfNeeded(i);
    const s1 = clean1[i] ?? orig1[i];
    const s2 = clean2[i] ?? orig2[i];
    if (!s1 || !s2) {
      r1.push(s1 ?? '');
      r2.push(s2 ?? '');
      prods.push('');
      noProducts++;
      continue;
    }
    try {
      const product = runTwoComponentSingle(rxn, rdkit, s1, s2);
      r1.push(s1);
      r2.push(s2);
      prods.push(product ? correctRGroups(product) : '');
      if (!product) noProducts++;
    } catch {
      errors++;
      r1.push(s1);
      r2.push(s2);
      prods.push('');
    }
  }

  if (showInfo) {
    const success = len - errors - noProducts;
    grok.shell.info(`Two-component (pairwise): ${success} products, ${errors} errors, ${noProducts} no reaction.`, {timeout: 10});
  }

  return {reactants1: r1, reactants2: r2, products: prods, errors, noProducts};
}

async function runMatrix(
  rxn: RDReaction, rdkit: RDModule,
  clean1: string[], clean2: string[],
  orig1: string[], orig2: string[],
  maxProducts: number, showInfo: boolean,
): Promise<TwoComponentReactionResult> {
  const r1: string[] = [];
  const r2: string[] = [];
  const prods: string[] = [];
  let errors = 0;
  let noProducts = 0;
  let count = 0;

  for (let i = 0; i < clean1.length && count < maxProducts; i++) {
    for (let j = 0; j < clean2.length && count < maxProducts; j++) {
      await yieldIfNeeded(count);
      count++;
      const s1 = clean1[i] ?? orig1[i];
      const s2 = clean2[j] ?? orig2[j];
      if (!s1 || !s2) {
        r1.push(s1 ?? '');
        r2.push(s2 ?? '');
        prods.push('');
        noProducts++;
        continue;
      }
      try {
        const product = runTwoComponentSingle(rxn, rdkit, s1, s2);
        r1.push(s1);
        r2.push(s2);
        prods.push(product ? correctRGroups(product) : '');
        if (!product) noProducts++;
      } catch {
        errors++;
        r1.push(s1);
        r2.push(s2);
        prods.push('');
      }
    }
  }

  if (showInfo) {
    const success = count - errors - noProducts;
    grok.shell.info(`Two-component (matrix ${clean1.length}x${clean2.length}): ${success} products, ${errors} errors, ${noProducts} no reaction.`, {timeout: 10});
  }

  return {reactants1: r1, reactants2: r2, products: prods, errors, noProducts};
}
