/**
 * The "Chemical engine" package setting (RDKit | Crux) and the one place that reads it. See `README.md` beside this
 * file: in Crux mode, what has moved to Crux runs on Crux, with RDKit behind every Crux path; other packages get RDKit
 * objects from Chem in both modes (`getRdKitModule`, `getRdKitService`).
 */
import * as DG from 'datagrok-api/dg';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_package} from '../package';
import {cruxModuleSync, cruxLoadFailed, loadCrux} from './crux-module';
import {cruxOrRdKit, cruxOrRdKitAsync, warnCruxFallback} from './fallback';

export enum ChemicalEngine {
  RDKit = 'RDKit',
  Crux = 'Crux',
}

/** The package property (package.json `properties`). */
const ENGINE_PROPERTY = 'ChemicalEngine';
/** The setting is read from the platform at most this often: the renderer asks for it on every cell. */
const SETTING_REFRESH_MS = 1000;

let engineOverride: ChemicalEngine | null = null;
let settingValue: ChemicalEngine = ChemicalEngine.RDKit;
let settingReadAt = -Infinity;
/** Changes whenever the engine in effect changes (a switch, Crux ready, Crux failed): caches keyed by engine compare
 * it to know they must be emptied. */
let generation = 0;
let lastInEffect: ChemicalEngine | null = null;
/** A load of the Crux module was asked for and has not settled yet. */
let waitingForCrux = false;
const switchListeners: (() => void)[] = [];

/**
 * A Chem package setting as the server has it (undefined when unset or unreadable). The package's own `_package` gives
 * it; in the package-test bundle, whose `_package` is not the published package, the published package does.
 */
export function chemSetting(name: string): any {
  const read = (settings: () => {[name: string]: any} | null | undefined) => {
    try {
      return settings() ?? null;
    } catch {
      return null;
    }
  };
  const settings = read(() => _package.settings) ??
    read(() => DG.Func.find({package: 'Chem', name: 'getRdKitModule'})[0]?.package?.settings);
  return settings?.[name];
}

/** Overrides the package setting until reset with null (tests, and scripts that run a suite in Crux mode). */
export function setChemicalEngine(engine: ChemicalEngine | null): void {
  engineOverride = engine;
  settingReadAt = -Infinity;
  engineInEffect();
}

/** The engine the setting (or the override) asks for. */
export function requestedEngine(): ChemicalEngine {
  if (engineOverride)
    return engineOverride;
  const now = performance.now();
  if (now - settingReadAt >= SETTING_REFRESH_MS) {
    settingReadAt = now;
    settingValue = chemSetting(ENGINE_PROPERTY) === ChemicalEngine.Crux ? ChemicalEngine.Crux : ChemicalEngine.RDKit;
  }
  return settingValue;
}

/** Whether Crux is asked for and has not failed to load: the Crux paths (or their RDKit fallbacks) run. */
export function isCruxRequested(): boolean {
  return requestedEngine() === ChemicalEngine.Crux && !cruxLoadFailed();
}

/**
 * The engine in effect now: Crux once its module is ready in Crux mode, else RDKit. Asking for it in Crux mode
 * starts loading the module; when it is ready (or fails), the engine in effect changes and the listeners hear it.
 */
export function engineInEffect(): ChemicalEngine {
  let engine = ChemicalEngine.RDKit;
  if (isCruxRequested()) {
    if (cruxModuleSync())
      engine = ChemicalEngine.Crux;
    else if (!waitingForCrux) {
      waitingForCrux = true;
      const settled = () => {
        waitingForCrux = false;
        engineInEffect();
      };
      loadCrux('engine').then(settled, settled);
    }
  }
  if (engine !== lastInEffect) {
    const first = lastInEffect === null;
    lastInEffect = engine;
    generation++;
    if (!first) {
      for (const listener of switchListeners)
        listener();
    }
  }
  return engine;
}

/** A number that changes with every change of {@link engineInEffect}. */
export function engineGeneration(): number {
  engineInEffect();
  return generation;
}

/** Calls `listener` whenever the engine in effect changes (not on the first read); returns its removal. */
export function onEngineSwitch(listener: () => void): () => void {
  switchListeners.push(listener);
  return () => {
    const i = switchListeners.indexOf(listener);
    if (i >= 0)
      switchListeners.splice(i, 1);
  };
}

/**
 * A call by the setting's engine (a worker job, a calculation, a script's replacement): in Crux mode `crux()`, and
 * `rdkit()` (the RDKit path: the RDKit pool, a script, the container) when it fails, logged as a fallback
 * ({@link cruxOrRdKitAsync}); in RDKit mode `rdkit()`.
 */
export function byEngineAsync<T>(call: string, crux: () => Promise<T>, rdkit: () => Promise<T>): Promise<T> {
  return isCruxRequested() ? cruxOrRdKitAsync(call, crux, rdkit) : rdkit();
}

/** The module drawing is done with: Crux's (RDKit.js-shaped) in Crux mode once it is ready, else RDKit's. */
export function renderingModule(rdKitModule: RDModule): {module: RDModule, engine: ChemicalEngine} {
  if (engineInEffect() === ChemicalEngine.Crux)
    return {module: cruxModuleSync()!, engine: ChemicalEngine.Crux};
  return {module: rdKitModule, engine: ChemicalEngine.RDKit};
}

/**
 * A one-molecule job on the main thread (a context panel's, a drawing) by the engine in effect: in Crux mode, once Crux
 * is ready, `job(Crux's module, true)`, and `job(rdKitModule, false)` where it gives null (a molecule Crux cannot read)
 * or throws, logged as a fallback (`fallback`: the log's message, when RDKit gave a result); in RDKit mode on
 * `rdKitModule`. A cache `job` fills holds one map per module: objects of one module never reach the other.
 */
export function byEngineModule<T>(call: string, rdKitModule: RDModule,
  job: (module: RDModule, crux: boolean) => T | null, fallback?: string): T | null {
  const {module, engine} = renderingModule(rdKitModule);
  if (engine !== ChemicalEngine.Crux)
    return job(rdKitModule, false);
  let threw = false; // (cruxOrRdKit logs a throw)
  const crux = cruxOrRdKit(call, () => job(module, true), () => {
    threw = true;
    return null;
  });
  if (crux !== null)
    return crux;
  const rdkit = job(rdKitModule, false);
  if (!threw && rdkit !== null)
    warnCruxFallback(call, fallback ?? 'Crux read no molecule, RDKit was asked');
  return rdkit;
}

/**
 * `read(mol)` of `molecule` parsed on one engine's `module` (`get_mol`, or `get_qmol` with `asQuery`), deleted after;
 * null where it reads none. Under Crux (`crux`, as {@link byEngineModule} hands it) what the module throws is thrown,
 * so that RDKit does the molecule (the fallback notes a trap); on RDKit it is null.
 */
export function withMol<T>(module: RDModule, crux: boolean, molecule: string, read: (mol: RDMol) => T | null,
  asQuery = false): T | null {
  let mol: RDMol | null = null;
  try {
    mol = asQuery ? module.get_qmol(molecule) : module.get_mol(molecule);
    return mol ? read(mol) : null;
  } catch (e) {
    if (crux)
      throw e;
    return null;
  } finally {
    mol?.delete();
  }
}

/** In Crux mode, waits for the Crux module (true when it is ready; false in RDKit mode or when it failed to load). */
export async function cruxReady(): Promise<boolean> {
  if (!isCruxRequested())
    return false;
  try {
    await loadCrux('engine');
  } catch {
    // logged by the loader; RDKit for the session
  }
  return engineInEffect() === ChemicalEngine.Crux;
}
