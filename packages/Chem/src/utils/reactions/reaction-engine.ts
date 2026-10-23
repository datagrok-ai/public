/**
 * Chem's reactions by the "Chemical Engine" setting (src/engine/README.md): in Crux mode, once the Crux module is ready
 * on the main thread, `get_rxn`, `run_reactants` and the molecules they read and write run on it (RDKit.js's calls;
 * the same products as RDKit on Chem's reactions, reaction files and enumeration demo: crux-core's decision
 * chem-wasm-ports), with RDKit behind them; on RDKit's module otherwise. Reaction drawing (the reaction renderer, the
 * dialogs' previews) stays RDKit's in both modes: crux has no `drawReaction` yet.
 */
import {RDModule, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ChemicalEngine, engineInEffect} from '../../engine/chem-engine';
import {cruxLoadFailed, cruxModuleSync, noteCruxError} from '../../engine/crux-module';
import {cruxOrRdKitAsync, warnCruxFallback} from '../../engine/fallback';

/** The engine the last reaction was set up on (tests). */
let lastEngine: ChemicalEngine | null = null;

/** The engine the last reaction (an EngineReaction, an enumeration) was set up on (tests). */
export function lastReactionEngine(): ChemicalEngine | null {
  return lastEngine;
}

/** Tests: SMILES the Crux module reads as no molecule in reactions (see {@link setCruxUnreadableReactantsForTests}). */
let unreadableForTests: Set<string> | null = null;

/** Tests: the Crux module reads these molecules as none in reactions (null: none); RDKit reads them as they are. */
export function setCruxUnreadableReactantsForTests(smiles: string[] | null): void {
  unreadableForTests = smiles ? new Set(smiles) : null;
}

/** The Crux module reactions run on in Crux mode (null: RDKit's, `rdkit`). */
function cruxReactionModule(): RDModule | null {
  const crux = engineInEffect() === ChemicalEngine.Crux ? cruxModuleSync() : null;
  lastEngine = crux ? ChemicalEngine.Crux : ChemicalEngine.RDKit;
  const unreadable = unreadableForTests;
  if (!crux || !unreadable)
    return crux;
  return {...crux, get_mol: (s: string, details?: string) => unreadable.has(s) ? null : crux.get_mol(s, details)} as
    unknown as RDModule;
}

/** `run(module)` on the engine's module: Crux's in Crux mode, and RDKit's (`rdkit`) when that fails. */
export async function onReactionEngine<T>(call: string, rdkit: RDModule, run: (module: RDModule) => Promise<T>):
  Promise<T> {
  const crux = cruxReactionModule();
  return crux ? cruxOrRdKitAsync(call, () => run(crux), () => run(rdkit)) : run(rdkit);
}

/** `module.get_rxn(smarts)`, or null where it gives none or throws. */
function getRxn(module: RDModule, smarts: string): RDReaction | null {
  try {
    return module.get_rxn(smarts) ?? null;
  } catch (e) {
    noteCruxError(e);
    return null;
  }
}

/**
 * One reaction, run molecule by molecule on the engine's module, with RDKit's behind Crux's per molecule: in Crux mode
 * `run` reacts on Crux, and on RDKit when Crux throws (a molecule it cannot read: logged once per call name) or after
 * the Crux module failed; a reaction Crux cannot read runs on RDKit. In RDKit mode, RDKit's as before.
 */
export class EngineReaction {
  private rdkitRxn: RDReaction | null = null;

  private constructor(private readonly call: string, private readonly rdkit: RDModule, private readonly smarts: string,
    private readonly crux: {module: RDModule, rxn: RDReaction} | null) {}

  /** The reaction of `smarts`, or null where RDKit reads none either (the callers' "invalid reaction"). */
  static create(call: string, rdkit: RDModule, smarts: string): EngineReaction | null {
    const module = cruxReactionModule();
    const rxn = module ? getRxn(module, smarts) : null;
    if (module && !rxn) {
      warnCruxFallback(call, `Crux reads no reaction from ${smarts}`);
      lastEngine = ChemicalEngine.RDKit;
    }
    const reaction = new EngineReaction(call, rdkit, smarts, rxn ? {module: module!, rxn} : null);
    return reaction.crux || reaction.rdkitReaction() ? reaction : null;
  }

  private rdkitReaction(): RDReaction | null {
    this.rdkitRxn ??= this.rdkit.get_rxn(this.smarts) ?? null;
    return this.rdkitRxn;
  }

  /** `react(module, rxn)` on Crux in Crux mode, on RDKit when it throws there; RDKit's error propagates. */
  run<T>(react: (module: RDModule, rxn: RDReaction) => T): T {
    if (this.crux && !cruxLoadFailed()) {
      try {
        return react(this.crux.module, this.crux.rxn);
      } catch (e) {
        noteCruxError(e);
        warnCruxFallback(this.call, e);
      }
    }
    const rxn = this.rdkitReaction();
    if (!rxn)
      throw new Error(`Invalid reaction SMARTS: "${this.smarts}".`);
    return react(this.rdkit, rxn);
  }

  delete(): void {
    try {this.crux?.rxn.delete();} catch {/* noop */}
    try {this.rdkitRxn?.delete();} catch {/* noop */}
  }
}
